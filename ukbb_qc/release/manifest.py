import argparse
import base64
import logging
import os
import subprocess
from typing import Dict, List, Tuple

import hail as hl

from google.cloud import storage

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("checksums")
logger.setLevel(logging.INFO)

PATH = "gs://broad-ukbb-requester-pays/broad.freeze_7/sharded_vcf"
"""
Path to requester pays bucket containing repackaged VCF shards and indices.
"""


def get_urls(in_file_path: str) -> Dict[str, str]:
    """
    Open file containing presigned URLs (created using `gsutil signurl`).

    Pre-signed URLs are required by DNAnexus for ROR.
    Pre-signed URls are created so that a user without a Google account can download data. For more information, see:
    https://cloud.google.com/storage/docs/access-control/signed-urls.

    This is an example command for how to create a pre-signed URL that lasts 1 second for a file in a requester pays bucket:
    ```
    gsutil -u maclab-ukbb signurl -r us-east1 -d 1s -b maclab-ukbb /Users/kchao/Desktop/maclab-ukbb-b14de5ad82e8.json gs://broad-ukbb-requester-pays/broad.freeze_6/sharded_vcf/*.tbi
    ```
    For more information on how to create a pre-signed URL, see https://cloud.google.com/storage/docs/gsutil/commands/signurl.

    :param str in_file_path: Input file path. Should be stored locally.
    :return: Dictionary containing file names (key) and presigned URLs (value).
    :rtype: Dict[str, str]
    """
    file_urls = {}
    with open(in_file_path) as i:
        _ = i.readline()
        for line in i:
            # Format:
            # URL\tHTTP Method\tExpiration\tSigned URL
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
        return i.readline().strip().split()


def generate_manifest(
    out_file: str,
    shards: List[str],
    indices: List[str],
    url_file_path: str,
    md5sums: Dict[str, str],
) -> None:
    """
    Generate manifest file for VCF shards.

    Read in signed URLs, header coordinates (required for DNAnexus), and 
    generate md5 for VCF shards or tabix indices.

    :param str out_file: Path to output file.
    :param List[str] shards: List of VCF shard file names.
    :param List[str] indices: List of VCF tabix index names.
    :param str url_file_path: Path to file containing presigned URLs.
    :param Dict[str, str] md5sums: Dictionary containing file name (key) and md5sum (value).
    :return: Nothing (writes to output file)
    :rtype: None
    """
    file_urls = get_urls(url_file_path)

    coordinate_temp_path = "gs://broad-ukbb/broad.freeze_7/temp/coordinates"
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
            md5 = md5sums[shard]

            # Write to output
            o.write(f"{shard}\t{url}\t{md5}\t{low}\t{high}\n")

        if indices:
            logger.info("Writing file details for tabix indices...")
            for index in indices:
                index_url = f"{PATH}/{index}"

                # Get presigned url
                url = file_urls[index_url]

                # Get md5 checksum
                md5 = md5sums[index]

                # Write to output
                o.write(f"{index}\t{url}\t{md5}\tNA\tNA\n")


def get_md5sums(bucket: str, prefix: str, gcloud_project: str) -> Dict[str, str]:
    """
    Get md5sums.

    Uses Google's Python API to get md5 checksums:
    https://googleapis.dev/python/storage/latest/client.html#google.cloud.storage.client.Client.list_blobs
    https://googleapis.dev/python/storage/latest/blobs.html

    Using `list_blobs` pulls md5 checksums within a few minutes for the entire bucket.

    :param str bucket: Top-level bucket that contains all VCF shards/indices.
    :param str prefix: Name of buckets that contain VCF shards/indices.
    :param str gcloud_project: Name of Google cloud project.
    :return: Dictionary of VCF/VCF index file names and their md5 checksums.
    :rtype: Dict[str, str]
    """
    md5sums = {}

    logger.info("Listing blobs (faster way to get md5sums)...")
    client = storage.Client()
    bucket = storage.Bucket(client, bucket, user_project=gcloud_project)
    # List all blobs ONLY in the freeze_7 VCF bucket
    all_blobs = list(client.list_blobs(bucket, prefix=prefix))
    for blob in all_blobs:
        file_name = blob.name
        # Skip anb Blobs that do not contain VCF shard or index information
        if "repackaged" not in file_name:
            continue
        md5sums[file_name] = base64.b64decode(blob.md5_hash).hex()
    return md5sums


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
    md5sums = get_md5sums(args.bucket, args.prefix, args.gcloud_project)
    generate_manifest(out_file, shards, indices, args.url, md5sums)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-b",
        "--bucket",
        help="Name of top-level bucket containing VCF shards and indices",
        default="broad-ukbb-requester-pays",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        help="Name of buckets within top-level bucket that contain VCF shards and indices",
        default="broad.freeze_7/sharded_vcf/",
    )
    parser.add_argument(
        "-g",
        "--gcloud-project",
        help="Name of Google project associated with VCF shards/indices.",
        default="maclab-ukbb",
    )
    parser.add_argument(
        "-u", "--url", help="Path to local TSV containing presigned URLs",
    )
    parser.add_argument(
        "-o", "--out", help="Path to output manifest file", required=True,
    )

    args = parser.parse_args()
    main(args)
