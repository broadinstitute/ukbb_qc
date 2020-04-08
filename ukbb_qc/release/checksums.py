import argparse
import os
import subprocess

from gnomad.utils.file_utils import get_file_stats
from ukbb_qc.resources.basics import release_vcf_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE


def main(args):

    # NOTE: Run this script locally
    data_source = "broad"
    freeze = args.freeze
    out_file = args.out

    # Get names of vcf shards
    path = release_vcf_path(data_source, freeze, contig=None)
    shards = (
        subprocess.check_output(["gsutil", "ls", path])
        .decode("utf8")
        .strip()
        .split("\n")
    )
    shards = [os.path.split(f)[-1] for f in shards if f.endswith(".bgz")]

    with open(out_file, "w") as o:
        # Write header
        o.write("shard\tsize\tmd5\n")

        # Write size and md5 for each shard
        for shard in shards:
            vcf_url = f"{release_vcf_path(data_source, freeze, contig=None)}/{shard}"
            size, int_size, md5 = get_file_stats(vcf_url)
            o.write(f"{shard}\t{size}\t{md5}\n")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "-o", "--out", help="Output file", required=True,
    )

    args = parser.parse_args()
    main(args)
