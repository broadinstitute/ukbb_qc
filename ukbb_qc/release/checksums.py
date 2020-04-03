import argparse
import base64
import json
import subprocess
from typing import Tuple
from gnomad.utils.file_utils import get_file_stats
from ukbb_qc.resources.basics import release_var_hist_path, release_vcf_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE


def main(args):

    # NOTE: Run this script locally
    data_source = args.data_source
    freeze = args.freeze
    shards = [str(c) for c in range(0, args.n_shards)]
    out_file = args.out

    with open(out_file, "w") as o:

        # write header
        o.write("shard\tsize\tmd5\n")

        for shard in shards:

            num_zeroes = (len(str(args.n_shards)) + 1) - len(shard)
            shard_name = f"part-{'0' * num_zeroes}{shard}.bgz"
            vcf_url = (
                f"{release_vcf_path(data_source, freeze, contig=None)}/{shard_name}"
            )

            size, int_size, md5 = get_file_stats(vcf_url)
            o.write(f"{shard_name}\t{size}\t{md5}\n")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s",
        "--data_source",
        help="Data source",
        choices=["regeneron", "broad"],
        default="broad",
    )
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "-n", "--n_shards", help="Number of VCF shards", default=5000, type=int,
    )
    parser.add_argument(
        "-o", "--out", help="Output file", required=True,
    )

    args = parser.parse_args()
    main(args)
