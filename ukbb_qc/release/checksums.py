import argparse
import logging
import os
import subprocess
from typing import List, Tuple

from gnomad.utils.file_utils import get_file_stats
from ukbb_qc.resources.basics import release_vcf_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("checksums")
logger.setLevel(logging.INFO)


def checksums(
    tranche_data: Tuple[str, str],
    out_file: str,
    shards: List[str],
    indices: List[str] = None,
) -> None:
    """
    Generate file stats (size and md5) for VCF shards or tabix indices.

    :param Tuple[str, str] tranche_data: Tranche data. Include data source str ('broad') and data freeze int.
    :param str out_file: Path to output file.
    :param List[str] shards: List of VCF shard file names.
    :param List[str] indices: List of VCF tabix index names. Default is None.
    :return: Nothing (writes to output file)
    :rtype: None
    """
    with open(out_file, "w") as o:
        # Write header
        o.write("shard\tsize\tmd5\n")

        # Write size and md5 for each shard
        logger.info("Writing file size and md5 for VCF shards...")
        for shard in shards:
            vcf_url = f"{release_vcf_path(*tranche_data, contig=None)}/{shard}"
            size, _, md5 = get_file_stats(vcf_url)
            o.write(f"{shard}\t{size}\t{md5}\n")

        if indices:
            logger.info("Writing file size and md5 for tabix indices...")
            for index in indices:
                index_url = f"{release_vcf_path(*tranche_data, contig=None)}/{index}"
                size, _, md5 = get_file_stats(index_url)
                o.write(f"{index}\t{size}\t{md5}\n")


def main(args):

    # NOTE: Run this script locally
    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)
    out_file = args.out

    # Get names of vcf shards
    path = release_vcf_path(*tranche_data, contig=None)
    files = (
        subprocess.check_output(["gsutil", "ls", path])
        .decode("utf8")
        .strip()
        .split("\n")
    )
    shards = [os.path.split(f)[-1] for f in files if f.endswith(".bgz")]

    # Get names of tabix indices if args.indices is set
    indices = None
    if args.indices:
        indices = [os.path.split(f)[-1] for f in files if f.endswith(".tbi")]

    logger.info("Generating checksums...")
    checksums(tranche_data, out_file, shards, indices)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "-o", "--out", help="Output file", required=True,
    )
    parser.add_argument(
        "--indices",
        help="Generate checksums for tabix indices in addition to VCF shards",
        action="store_true",
    )

    args = parser.parse_args()
    main(args)
