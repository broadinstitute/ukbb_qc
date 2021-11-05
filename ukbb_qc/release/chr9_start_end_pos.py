import argparse
import logging
import os
import subprocess

import hail as hl
import hailtop.batch as hb

from gnomad.utils.file_utils import parallel_file_exists


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("get_chr9_shard_positions")
logger.setLevel(logging.INFO)


def get_start_end_pos(shard_path: str, pos_path: str) -> None:
    """
    Import VCF shard as MT, filter to rows, and save start and end position.

    Write start and end position to output temporary file.

    :param str shard_path: Path to input VCF shard.
    :param str coordinates_path: Path to output temporary file to store start/end positions.
    :return: None
    """
    vcf = hl.import_vcf(shard_path, reference_genome="GRCh38").rows()
    vcf = vcf.select()
    head = vcf.head(1).locus.position.take(1)[0]
    tail = vcf.tail(1).locus.position.take(1)[0]
    logger.info("Start %i, end %i", head, tail)
    shard_name = os.path.basename(shard_path)

    with hl.hadoop_open(pos_path, "w") as o:
        o.write(f"{shard_name}\t{head}\t{tail}\n")


def main(args):

    backend = hb.ServiceBackend("gnomad-production", "gnomad-kc")
    b = hb.Batch(
        backend=backend,
        default_python_image=args.image,
        default_memory=args.memory,
        default_storage=args.storage,
        default_cpu=args.cpu,
    )

    logger.info("Getting VCF shard paths...")
    path = "gs://broad-ukbb/broad.freeze_7/release/vcf/broad.freeze_7.chr9.vcf.bgz/*bgz"
    shards = (
        subprocess.check_output(["gsutil", "ls", path])
        .decode("utf8")
        .strip()
        .split("\n")
    )
    logger.info("Found %i shards", len(shards))

    logger.info("Checking if any output files already exist...")
    temp_path = "gs://broad-ukbb/broad.freeze_7/temp/"
    pos_file_map = {}
    for shard in shards:
        # Skip chr9 part-00000.bgz shard (it only has a header)
        if (
            shard
            == "gs://broad-ukbb/broad.freeze_7/release/vcf/broad.freeze_7.chr9.vcf.bgz/part-00000.bgz"
        ):
            continue
        output_file_name = f"{os.path.basename(shard).split('.')[0]}_pos.txt"
        pos_file_map[shard] = f"{temp_path}/{output_file_name}"
    output_pos_files_exist = parallel_file_exists(list(pos_file_map.values()))

    shards_to_check = []
    for shard in pos_file_map:
        if not output_pos_files_exist[pos_file_map[shard]]:
            shards_to_check.append(shard)
    logger.info("Found %i shards to check", len(shards_to_check))

    logger.info("Preparing to start batch jobs...")
    for shard in shards_to_check:
        logger.info("Working on %s", shard)
        j = b.new_python_job(name=shard)
        j.call(
            get_start_end_pos, shard, pos_file_map[shard],
        )
    b.run(wait=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--memory", help="Memory for batch job. Default is 16 GiB.", default="16Gi",
    )
    parser.add_argument(
        "--storage", help="Disk space for batch job. Default is 10GiB.", default="10Gi",
    )
    parser.add_argument(
        "--cpu",
        help="Number of CPUs for batch job. Default is 2.",
        default=2,
        type=int,
    )
    parser.add_argument(
        "--image",
        help="Docker image. Must have python, hail, dill installed.",
        default="gcr.io/broad-mpg-gnomad/tgg-methods-vm:20210915",
    )
    args = parser.parse_args()
    main(args)
