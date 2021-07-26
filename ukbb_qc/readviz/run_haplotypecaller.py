import logging
import os
import pandas as pd
import tqdm

import hail as hl

from tgg.batch.batch_utils import (
    check_storage_bucket_region,
    HG38_REF_PATHS,
    localize_file,
    init_arg_parser,
    init_job,
    run_batch,
    set_gcloud_project,
)

from ukbb_qc.resources.basics import readviz_haplotype_caller_path


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("run_haplotypecaller")
logger.setLevel(logging.INFO)


GCLOUD_PROJECT = "broad-mpg-gnomad"
"""
Google cloud project.
"""

# Dockerfile for this image is in gnomad-readviz repo
DOCKER_IMAGE = "gcr.io/broad-mpg-gnomad/gnomad-readviz:ukbb"
"""
Docker image to be passed to hail batch. 

Built with GATK version 4.0.10.1 (version used to process UKBB crams).
"""

EXCLUDE_INTERVALS = (
    "gs://gnomad-bw2/exclude_intervals_with_non_ACGT_bases_in_GRCh38__150bp_window.bed"
)
"""
Variants to exclude when running HaplotypeCaller.
"""

PADDING_AROUND_VARIANT = 200
"""
Amount of padding to add around each variant when running HaplotypeCaller.
"""


def parse_args():
    """Parse command line args."""

    p = init_arg_parser(default_cpu=1, default_billing_project="gnomad-production")
    p.add_argument(
        "-p",
        "--output-dir",
        help="Where to write haplotype caller output.",
        default=f"{readviz_haplotype_caller_path}/outputs/",
    )
    p.add_argument(
        "-n",
        "--num-samples-to-process",
        help="For testing, process only the first N samples.",
        type=int,
    )
    p.add_argument(
        "-s",
        "--sample-to-process",
        help="For testing, process only the given sample id(s).",
        nargs="+",
    )
    p.add_argument(
        "cram_and_tsv_paths_table",
        help="A text file containing at least these columns: sample_id, cram_path",
        default=f"{readviz_haplotype_caller_path()}/inputs/step4_output_cram_and_tsv_paths_table.tsv",
    )
    args = p.parse_args()

    return p, args


def main():
    """
    Run HaplotypeCaller to generate bamouts.

    Step 5 of readviz pipeline.
    """
    p, args = parse_args()

    df = pd.read_table(args.cram_and_tsv_paths_table)
    if {"sample_id", "cram_path", "crai_path", "variants_tsv_bgz"} - set(df.columns):
        p.error(f"{args.tsv_path} must contain 'sample_id', 'cram_path' columns")

    # Check that all buckets are in "US-CENTRAL1" or are multi-regional to avoid egress charges to the Batch cluster
    set_gcloud_project(GCLOUD_PROJECT)
    if args.cluster:
        check_storage_bucket_region(df.cram_path)

    if not args.force:
        hl.init(log="/dev/null", quiet=True)

    # Process samples
    with run_batch(args, batch_name=f"HaplotypeCaller -bamout") as batch:
        counter = 0
        for _, row in tqdm.tqdm(df.iterrows(), unit=" rows", total=len(df)):
            if args.sample_to_process and row.sample_id not in set(
                args.sample_to_process
            ):
                continue

            input_filename = os.path.basename(row.cram_path)
            output_prefix = input_filename.replace(".bam", "").replace(".cram", "")

            output_bam_path = os.path.join(
                args.output_dir, f"{output_prefix}.bamout.bam"
            )
            output_bai_path = os.path.join(
                args.output_dir, f"{output_prefix}.bamout.bai"
            )

            if (
                not args.force
                and hl.hadoop_is_file(output_bam_path)
                and hl.hadoop_is_file(output_bai_path)
            ):
                logger.info(
                    f"Output files exist (eg. {output_bam_path}). Skipping {input_filename}..."
                )
                continue

            counter += 1
            if args.num_samples_to_process and counter > args.num_samples_to_process:
                break

            j = init_job(
                batch,
                f"readviz: {row.sample_id}",
                DOCKER_IMAGE if not args.raw else None,
                args.cpu,
                args.memory,
            )
            j.command(
                f"""gcloud -q auth activate-service-account --key-file=/gsa-key/key.json"""
            )
            local_exclude_intervals = localize_file(j, EXCLUDE_INTERVALS)
            local_fasta = localize_file(j, HG38_REF_PATHS.fasta, use_gcsfuse=True)
            local_fasta_fai = localize_file(j, HG38_REF_PATHS.fai, use_gcsfuse=True)
            localize_file(j, HG38_REF_PATHS.dict, use_gcsfuse=True)
            local_tsv_bgz = localize_file(j, row.variants_tsv_bgz)
            local_cram_path = localize_file(j, row.cram_path)
            local_crai_path = localize_file(j, row.crai_path)

            j.command(
                f"""echo --------------

echo "Start - time: $(date)"
df -kh


# 1) Convert variants_tsv_bgz to sorted interval list

gunzip -c "{local_tsv_bgz}" | awk '{{ OFS="\t" }} {{ print( "chr"$1, $2, $2 ) }}' | bedtools slop -b {PADDING_AROUND_VARIANT} -g {local_fasta_fai} > variant_windows.bed

# Sort the .bed file so that chromosomes are in the same order as in the input_cram file.
# Without this, if the input_cram has a different chromosome ordering (eg. chr1, chr10, .. vs. chr1, chr2, ..)
# than the interval list passed to GATK tools' -L arg, then GATK may silently skip some of regions in the -L intervals.
# The sort is done by first retrieving the input_cram header and passing it to GATK BedToIntervalList.

java -Xms2g -jar /gatk/gatk.jar PrintReadsHeader \
    --gcs-project-for-requester-pays {GCLOUD_PROJECT} \
    -R {local_fasta} \
    -I "{local_cram_path}" \
    -O header.bam

java -Xms2g -jar /gatk/gatk.jar BedToIntervalList \
    --SORT true \
    --SEQUENCE_DICTIONARY header.bam \
    --INPUT variant_windows.bed \
    --OUTPUT variant_windows.interval_list

# 2) Get reads from the input_cram for the intervals in variant_windows.interval_list

time java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+DisableAttachMechanism -XX:MaxHeapSize=2000m -Xmx30000m \
    -jar /gatk/GATK35.jar \
    -T HaplotypeCaller \
    -R {local_fasta} \
    -I "{local_cram_path}" \
    -L variant_windows.interval_list \
    -XL {local_exclude_intervals} \
    --disable_auto_index_creation_and_locking_when_reading_rods \
    -ERC GVCF \
    --max_alternate_alleles 3 \
    -variant_index_parameter 128000 \
    -variant_index_type LINEAR \
    --read_filter OverclippedRead \
    -bamout "{output_prefix}.bamout.bam" \
    -o "{output_prefix}.gvcf"  |& grep -v "^DEBUG"

bgzip "{output_prefix}.gvcf"
tabix "{output_prefix}.gvcf.gz"

gsutil -m cp "{output_prefix}.bamout.bam" {args.output_dir}
gsutil -m cp "{output_prefix}.bamout.bai" {args.output_dir}
gsutil -m cp "{output_prefix}.gvcf.gz" {args.output_dir}
gsutil -m cp "{output_prefix}.gvcf.gz.tbi" {args.output_dir}

ls -lh
echo --------------; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------

"""
            )


if __name__ == "__main__":
    main()
