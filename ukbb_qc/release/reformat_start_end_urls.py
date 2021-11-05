import argparse
import logging

from gnomad.utils.file_utils import get_file_stats

from ukbb_qc.release.manifest import get_urls

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("reformat_vcf_shard_start_end_urls")
logger.setLevel(logging.INFO)


def main(args):

    logger.info("Getting input URLs...")
    file_urls = get_urls(args.input_url_path)
    with open(args.output_tsv_path, "w") as o:
        o.write("file_name\tlink\tmd5sum\n")
        for tsv in file_urls:
            # Get md5 checksum
            _, _, md5 = get_file_stats(tsv, args.project_id)
            o.write(f"{tsv}\t{file_urls[tsv]}\t{md5}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-url-path",
        help="Path to input local file with pre-signed URLs for TSVs containing VCF shard start/end positions.",
    )
    parser.add_argument(
        "--output-tsv-path",
        help="Path to local output file. Required for file stats function because TSVs with shard start/end positions are in a requester-pays bucket.",
    )
    parser.add_argument("--project-id", help="Project ID.")
    args = parser.parse_args()
    main(args)
