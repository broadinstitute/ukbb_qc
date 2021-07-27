import asyncio
import logging
from tqdm import tqdm

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import parallel_file_exists

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
        "--gcloud-project",
        help="Google cloud project. Default is 'broad-mpg-gnomad'.",
        default="broad-mpg-gnomad",
    )
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
        "--cram-and-tsv_paths-table",
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
    _, args = parse_args()
    hl.init(log="/dev/null", quiet=True)
    project = args.gcloud_project
    set_gcloud_project(project)

    logger.info("Making sure input cram_and_tsv_paths_table arg is valid...")
    bams = {}
    samples = {}
    with hl.hadoop_open(args.cram_and_tsv_paths_table) as c:
        # Confirm header has all required columns
        header = c.readline().strip().split("\t")
        if {"sample_id", "cram_path", "crai_path", "variants_tsv_bgz"} - set(header):
            raise DataException(
                "%s must contain 'sample_id', 'cram_path', 'crai_path', variants_tsv_bgz' columns!"
            )

        for line in c:
            sample, cram, crai, variants_tsv_bgz = line.strip().split("\t")

            # Store output BAM path
            bam = f"{args.output_prefix}/{sample}.bamout.bam"
            bai = f"{args.output_prefix}/{sample}.bamout.bai"
            bams[sample] = bam

            # Store sample information
            samples[sample] = [cram, crai, variants_tsv_bgz, bam, bai]

            logger.info(
                "Checking that all input crams are 'US-CENTRAL1' or multi-regional buckets..."
            )
            # Check that all buckets are in "US-CENTRAL1" or are multi-regional to avoid egress charges to the Batch cluster
            check_storage_bucket_region(cram)

    logger.info("Checking if any output bams already exist...")
    bam_exists = asyncio.get_event_loop().run_until_complete(
        parallel_file_exists(list(bams.values()))
    )
    samples_without_bams = []
    for sample in bams:
        if not bam_exists[bams[sample]]:
            samples_without_bams.append(sample)

    # Process samples
    with run_batch(args, batch_name=f"HaplotypeCaller -bamout") as batch:
        with tqdm(
            range(len(samples)), desc="submit HC batch job for samples",
        ) as pbar:
            for sample in samples:
                cram, crai, variants_tsv_bgz, bam, bai = samples[sample]

                j = init_job(
                    batch,
                    f"readviz: {sample}",
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
                local_tsv_bgz = localize_file(j, variants_tsv_bgz)
                local_cram_path = localize_file(j, cram)

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
    --gcs-project-for-requester-pays {project} \
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
    -bamout "{sample}.bamout.bam" \
    -o "{sample}.gvcf"  |& grep -v "^DEBUG"

bgzip "{sample}.gvcf"
tabix "{sample}.gvcf.gz"

gsutil -m cp "{sample}.bamout.bam" {args.output_dir}
gsutil -m cp "{sample}.bamout.bai" {args.output_dir}
gsutil -m cp "{sample}.gvcf.gz" {args.output_dir}
gsutil -m cp "{sample}.gvcf.gz.tbi" {args.output_dir}

ls -lh
echo --------------; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------

"""
                )


if __name__ == "__main__":
    main()
