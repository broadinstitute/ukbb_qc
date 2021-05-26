import argparse
import logging
from typing import List

import hailtop.batch as hb
import hail as hl

from gnomad.utils.slack import slack_notifications

from ukbb_qc.resources.basics import (
    readviz_ht_exploded_path,
    readviz_per_sample_tsv_path,
)
from ukbb_qc.slack_creds import slack_token


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("run_per_sample_tsv_export")
logger.setLevel(logging.INFO)


def get_sample_ids(ids_file: str, header: bool = False) -> List[str]:
    """
    Open file with all sample IDs and stores IDs in a list.

    :param str ids_file: Path to file containing sample IDs with variants to extract for readviz.
    :param bool header: Whether IDs file has a header line. Default is False.
    :return: List of sample IDs
    :rtype: List[str]
    """
    sample_ids = []
    with open(ids_file) as i:
        if header:
            header = i.readline()
        for line in i:
            sample_ids.append(line.strip())
    return sample_ids


def export_tsv(ht_path: str, sample_id: str, tsv_path: str) -> None:
    """
    Read in hail Table, filter to specified sample ID, and export TSV.

    :param str ht_path: Path to input hail Table.
    :param str sample_id: Sample for which to export a TSV.
    :param str tsv_path: Path to output TSV.
    :return: None
    """
    ht = hl.read_table(ht_path)
    ht = ht.filter(ht.S == sample_id)
    ht = ht.naive_coalesce(1)
    ht.export(tsv_path)


def main(args):

    backend = hb.ServiceBackend("gnomad-production", "gs://gnomad-kc/temp/")
    b = hb.Batch(
        backend=backend,
        default_cpu=1,
        default_python_image="hailgenetics/hail:0.2.67",
        project="broad-mpg-gnomad",
    )

    logger.info("Extracting sample IDs...")
    sample_ids = get_sample_ids(args.ids_file, args.header)

    logger.info("Preparing to start batch job...")
    for sample in sample_ids:
        logger.info("Working on %s", sample)
        j = b.new_python_job(name=sample)
        j.call(
            export_tsv,
            readviz_ht_exploded_path(),
            sample,
            f"{readviz_per_sample_tsv_path()}/{sample}.tsv",
        )

    b.run(wait=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ids-file", help="File with sample IDs")
    parser.add_argument(
        "--header",
        help="Whether file with sample IDs has a header line",
        action="store_true",
    )
    parser.add_argument("--slack-channel", help="Send message to Slack channel/user")
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
