import logging
from tqdm import tqdm

import hail as hl

from gnomad.utils.file_utils import parallel_file_exists

from tgg.batch.batch_utils import (
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
logger = logging.getLogger("check_haplotypecaller_bamouts")
logger.setLevel(logging.INFO)


# Dockerfile for this image is in gnomad-readviz repo
DOCKER_IMAGE = "gcr.io/broad-mpg-gnomad/gnomad-readviz:ukbb"
"""
Docker image to be passed to hail batch. 

Built with GATK version 4.0.10.1 (version used to process UKBB crams).
Must contain samtools package.
"""


def main(args):
    """
    Check HaplotypeCaller (HC) bamouts to make sure they are not malformed.

    Necessary only because batch cluster had to be frozen during initial HC job, and ~4k jobs errored out.
    """
    hl.init(log="/dev/null", quiet=True)

    logger.info("Getting sample IDs to check and their corresponding bamouts...")
    bamouts = {}
    with open(args.ids_file) as i:
        for line in i:
            sample = line.strip()
            bamouts[
                sample
            ] = f"{readviz_haplotype_caller_path}/outputs/{sample}.bamout.bam"

    logger.info("Running file existence check on bamouts...")
    bamout_files_exist = parallel_file_exists(list(bamouts.values()))
    samples_without_bamout_files = []
    samples_with_bamout_files = []
    for sample in bamouts:
        if not bamout_files_exist[bamouts[sample]]:
            samples_without_bamout_files.append(sample)
        else:
            samples_with_bamout_files.append(sample)

    logger.info("Samples with output bam files: %i", len(samples_with_bamout_files))
    logger.info(
        "Samples without output bam files: %i", len(samples_without_bamout_files)
    )

    # Process samples
    with run_batch(args, batch_name="HaplotypeCaller -bamout checks") as batch:
        for sample, bamout in tqdm(bamouts.items(), unit="samples"):

            j = init_job(
                batch,
                f"{sample} bamout check",
                DOCKER_IMAGE if not args.raw else None,
                args.cpu,
                args.memory,
            )
            j.command(
                f"""gcloud -q auth activate-service-account --key-file=/gsa-key/key.json"""
            )
            local_bamout_path = localize_file(j, bamout, use_gcsfuse=True)
            j.command(f"""samtools quickcheck {local_bamout_path}""")


if __name__ == "__main__":
    p = init_arg_parser(default_cpu=1, default_billing_project="gnomad-production")
    p.add_argument(
        "--gcloud-project",
        help="Google cloud project. Default is 'gnomad-production'.",
        default="gnomad-production",
    )
    p.add_argument(
        "--ids-file", help="A text file containing sample IDs to check",
    )
    args = p.parse_args()
    main(args)
