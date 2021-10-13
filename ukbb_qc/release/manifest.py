import argparse
import logging
import os
import subprocess
from typing import Dict, List, Tuple

import hail as hl

from gnomad.utils.file_utils import get_file_stats

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("checksums")
logger.setLevel(logging.INFO)


PATH = "gs://broad-ukbb-requester-pays/broad.freeze_7/sharded_vcf/"
"""
Path to requester pays bucket containing repackaged VCF shards and indices.
"""


def get_urls(in_file_path: str) -> Dict[str, str]:
    """
    Open file containing presigned URLs (created using `gsutil signurl`).

    Coordinates are required by DNAnexus.

    :param str in_file_path: Input file path. Should be stored locally.
    :return: Dictionary containing file names (key) and presigned URLs (value).
    :rtype: Dict[str, str]
    """
    file_urls = {}
    with open(in_file_path) as i:
        _ = i.readline()
        for line in i:
            # Format:
            # URL   HTTP Method Expiration  Signed URL
            file, _, _, url = line.strip().split("\t")
            file_urls[file] = url
    return file_urls


def get_coordinates(in_cloud_file_path: str) -> Tuple[int, int]:
    """
    Open temporary file containing VCF header coordinates (printed from `ukbb_header_reformat.sh`).

    Coordinates are required by DNAnexus.

    :param str in_cloud_file_path: Input file path. Should be in Google cloud storage.
    :return: Low and high header coordinates.
    :rtype: Tuple[int, int]
    """
    with hl.hadoop_open(in_cloud_file_path) as i:
        low, high = i.readline().strip().split()
    return (low, high)


def generate_manifest(
    out_file: str,
    shards: List[str],
    indices: List[str],
    url_file_path: str,
    project_id: str = "maclab-ukbb",
) -> None:
    """
    Generate manifest file for VCF shards.

    Read in signed URLs, header coordinates (required for DNAnexus), and 
    generate file stats (size and md5) for VCF shards or tabix indices.

    :param str out_file: Path to output file.
    :param List[str] shards: List of VCF shard file names.
    :param List[str] indices: List of VCF tabix index names.
    :param str url_file_path: Path to file containing presigned URLs.
    :param str project_id: Google cloud project. Required for file stats function since VCF shards and indices
        are in requester-pays bucket.
    :return: Nothing (writes to output file)
    :rtype: None
    """
    file_urls = get_urls(url_file_path)

    coordinate_temp_path = "gs://broad-ukbb/broad.freeze_7/temp"
    with open(out_file, "w") as o:
        # Write header
        o.write("file_name\tlink\tmd5sum\theader_lo\theader_hi\n")

        # Write size and md5 for each shard
        logger.info("Writing file details for VCF shards...")
        for shard in shards:
            vcf_url = f"{PATH}/{shard}"

            # Get presigned url
            url = file_urls[vcf_url]

            # Get coordinates
            low, high = get_coordinates(
                f"{coordinate_temp_path}/{shard}_coordinates.tsv"
            )

            # Get md5 checksum
            _, _, md5 = get_file_stats(vcf_url, project_id)

            # Write to output
            o.write(f"{shard}\t{url}\t{md5}\t{low}\t{high}\n")

        if indices:
            logger.info("Writing file details for tabix indices...")
            for index in indices:
                index_url = f"{PATH}/{index}"

                # Get presigned url
                url = file_urls[vcf_url]

                # Get md5 checksum
                _, _, md5 = get_file_stats(index_url, project_id)

                # Write to output
                o.write(f"{index}\t{url}\t{md5}\tNA\tNA\n")


def main(args):
    """
    Generate manifest file required to transfer data to DNAnexus.
    """

    # NOTE: Run this script locally
    out_file = args.out

    # Get names of vcf shards
    files = (
        subprocess.check_output(["gsutil", "-u", "maclab-ukbb", "ls", PATH])
        .decode("utf8")
        .strip()
        .split("\n")
    )
    shards = [os.path.split(f)[-1] for f in files if f.endswith(".gz")]

    # Get names of tabix indices if args.indices is set
    indices = [os.path.split(f)[-1] for f in files if f.endswith(".tbi")]

    logger.info("Generating manifest file...")
    generate_manifest(out_file, shards, indices, args.url)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-u", "--url", help="Path to local TSV containing presigned URLs",
    )
    parser.add_argument(
        "-o", "--out", help="Path to output manifest file", required=True,
    )

    args = parser.parse_args()
    main(args)