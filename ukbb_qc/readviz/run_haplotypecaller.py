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
        help="Google cloud project. Default is 'gnomad-production'.",
        default="gnomad-production",
    )
    p.add_argument(
        "-p",
        "--output-dir",
        help="Where to write haplotype caller output.",
        default=f"{readviz_haplotype_caller_path}/outputs/",
    )
    p.add_argument(
        "--cram-and-tsv-paths-table",
        help="A text file containing at least these columns: sample_id, cram_path, variants_tsv_bgz",
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
    output_dir = args.output_dir

    logger.info("Making sure input cram_and_tsv_paths_table arg is valid...")
    bamouts = {}
    samples = {}
    with hl.hadoop_open(args.cram_and_tsv_paths_table) as c:
        # Confirm header has all required columns
        header = c.readline().strip().split("\t")
        if {"sample_id", "cram", "crai", "variants_tsv_bgz"} - set(header):
            raise DataException(
                "%s must contain 'sample_id', 'cram', 'crai', 'variants_tsv_bgz' columns!"
            )

        for line in c:
            values = dict(zip(header, line.strip().split("\t")))

            # Store output BAM path
            bamout = f"{args.output_dir}/{values['sample_id']}.bamout.bam"
            bamouts[values["sample_id"]] = bamout

            # Store sample information
            samples[values["sample_id"]] = [
                values["cram"],
                values["crai"],
                values["variants_tsv_bgz"],
            ]

            logger.info(
                "Checking that all input crams are 'US-CENTRAL1' or multi-regional buckets..."
            )
            # Check that all buckets are in "US-CENTRAL1" or are multi-regional to avoid egress charges to the Batch cluster
            check_storage_bucket_region(values["cram"])

    logger.info("Checking if any output bams already exist...")
    bamout_exists = parallel_file_exists(list(bamouts.values()))
    samples_without_bamouts = []
    for sample in bamouts:
        if not bamout_exists[bamouts[sample]]:
            samples_without_bamouts.append(sample)

    # Process samples
    with run_batch(args, batch_name="HaplotypeCaller -bamout") as batch:
        for sample in tqdm(samples_without_bamouts, unit="samples"):
            cram, crai, variants_tsv_bgz = samples[sample]

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
            localize_file(j, crai)

            # NOTE: Currently, we don't store the gVCFs produced by this step because it isn't used
            # TODO: Revisit if we should store gVCFs for downstream comparison --
            # we should keep gVCFs if there is a cost-effective way to compare this gVCF call with
            # variant call in full callset
            j.command(
                f"""echo --------------

echo "Start - time: $(date)"
df -kh


# 1) Convert variants_tsv_bgz to sorted interval list

gunzip -c "{local_tsv_bgz}" | grep -v chrom| awk '{{ OFS="\t" }} {{ print( "chr"$6, $7, $7 ) }}' | bedtools slop -b {PADDING_AROUND_VARIANT} -g {local_fasta_fai} > variant_windows.bed

# Sort the .bed file so that chromosomes are in the same order as in the input_cram file.
# Without this, if the input_cram has a different chromosome ordering (eg. chr1, chr10, .. vs. chr1, chr2, ..)
# than the interval list passed to GATK tools' -L arg, then GATK may silently skip some of regions in the -L intervals.
# The sort is done by first retrieving the input_cram header and passing it to GATK BedToIntervalList.

samtools view -H {local_cram_path} -T {local_fasta} > header.bam

java -Xms2g -jar /gatk/gatk.jar BedToIntervalList \
    --SORT true \
    --SEQUENCE_DICTIONARY header.bam \
    --INPUT variant_windows.bed \
    --OUTPUT variant_windows.interval_list

# 2) Get reads from the input_cram for the intervals in variant_windows.interval_list

time java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+DisableAttachMechanism -XX:MaxHeapSize=2000m -Xmx30000m \
    -jar /gatk/gatk.jar \
    HaplotypeCaller \
    -R {local_fasta} \
    -I "{local_cram_path}" \
    -L variant_windows.interval_list \
    -XL {local_exclude_intervals} \
    -ERC GVCF \
    --gvcf-gq-bands 10 \
    --gvcf-gq-bands 20 \
    --gvcf-gq-bands 30 \
    --gvcf-gq-bands 40 \
    --gvcf-gq-bands 50 \
    --gvcf-gq-bands 60 \
    --gvcf-gq-bands 70 \
    --gvcf-gq-bands 80 \
    --gvcf-gq-bands 90 \
    --use-new-qual-calculator true \
    --contamination-fraction-to-filter 0.0069722000000000004 \
    --annotation-group StandardAnnotation \
    --annotation-group StandardHCAnnotation \
    --annotation-group AS_StandardAnnotation \
    --indel-size-to-eliminate-in-ref-model 10 \
    --use-alleles-trigger false \
    --disable-optimizations false \
    --just-determine-active-regions false \
    --dont-genotype false \
    --max-mnp-distance 0 \
    --max-disc-ar-extension 25 \
    --max-gga-ar-extension 300 \
    --kmer-size 10 \
    --kmer-size 25 \
    --dont-increase-kmer-sizes-for-cycles false \
    --allow-non-unique-kmers-in-ref false \
    --num-pruning-samples 1 \
    --recover-dangling-heads false \
    --do-not-recover-dangling-branches false \
    --min-dangling-branch-length 4 \
    --consensus false \
    --max-num-haplotypes-in-population 128 \
    --error-correct-kmers false \
    --min-pruning 2 \
    --debug-graph-transformations false \
    --kmer-length-for-read-error-correction 25 \
    --min-observations-for-kmer-to-be-solid 20 \
    --likelihood-calculation-engine PairHMM \
    --base-quality-score-threshold 18 \
    --pair-hmm-gap-continuation-penalty 10 \
    --pair-hmm-implementation FASTEST_AVAILABLE \
    --pcr-indel-model CONSERVATIVE \
    --phred-scaled-global-read-mismapping-rate 45 \
    --native-pair-hmm-threads 4 \
    --native-pair-hmm-use-double-precision false \
    --debug false \
    --use-filtered-reads-for-annotations false \
    --bam-writer-type CALLED_HAPLOTYPES \
    --dont-use-soft-clipped-bases false \
    --capture-assembly-failure-bam false \
    --error-correct-reads false \
    --do-not-run-physical-phasing false \
    --min-base-quality-score 10 \
    --smith-waterman JAVA \
    --annotate-with-num-discovered-alleles false \
    --heterozygosity 0.001 \
    --indel-heterozygosity 1.25E-4 \
    --heterozygosity-stdev 0.01 \
    --standard-min-confidence-threshold-for-calling 10.0 \
    --max-alternate-alleles 6 \
    --max-genotype-count 1024 \
    --sample-ploidy 2 \
    --num-reference-samples-if-no-call 0 \
    --genotyping-mode DISCOVERY \
    --genotype-filtered-alleles false \
    --output-mode EMIT_VARIANTS_ONLY \
    --all-site-pls false \
    --min-assembly-region-size 50 \
    --max-assembly-region-size 300 \
    --assembly-region-padding 100 \
    --max-reads-per-alignment-start 50 \
    --active-probability-threshold 0.002 \
    --max-prob-propagation-distance 50 \
    --interval-set-rule UNION \
    --interval-padding 0 \
    --interval-exclusion-padding 0 \
    --interval-merging-rule ALL \
    --read-validation-stringency SILENT \
    --seconds-between-progress-updates 10.0 \
    --disable-sequence-dictionary-validation false \
    --create-output-bam-index true \
    --create-output-bam-md5 false \
    --create-output-variant-index true \
    --create-output-variant-md5 false \
    --lenient false \
    --add-output-sam-program-record true \
    --add-output-vcf-command-line true \
    --cloud-prefetch-buffer 40 \
    --cloud-index-prefetch-buffer -1 \
    --disable-bam-index-caching false \
    --sites-only-vcf-output false \
    --help false \
    --version false \
    --showHidden false \
    --verbosity INFO \
    --QUIET false \
    --use-jdk-deflater false \
    --use-jdk-inflater false\
    --gcs-max-retries 20 \
    --gcs-project-for-requester-pays {project} \
    --disable-tool-default-read-filters false \
    --minimum-mapping-quality 20 \
    --disable-tool-default-annotations false \
    --enable-all-annotations false \
    --dont-trim-active-regions false \
    -bamout "{sample}.bamout.bam" \
    -O "{sample}.gvcf"  |& grep -v "^DEBUG"

gsutil -m cp "{sample}.bamout.bam" {output_dir}
gsutil -m cp "{sample}.bamout.bai" {output_dir}

ls -lh
echo --------------; free -h; df -kh; uptime; set +xe; echo "Done - time: $(date)"; echo --------------

"""
            )


if __name__ == "__main__":
    main()
