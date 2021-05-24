import argparse
import logging
from typing import List

import hailtop.batch as hb

from gnomad.utils.slack import slack_notifications

from ukbb_qc.resources.basics import readviz_ht_path, readviz_per_sample_tsv_path
from ukbb_qc.slack_creds import slack_token


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("run_per_sample_tsv_export")
logger.setLevel(logging.INFO)


def get_sample_ids(ids_file: str, header: bool = False) -> List[str]:
    """
    Open file with all sample IDs and stores IDs in a list.

    :param str ids_file: List of sample IDs with variants to extract for readviz.
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


def main(args):

    backend = hb.ServiceBackend("kchao-trial", "gs://gnomad-kc/temp/")
    b = hb.Batch(backend=backend)

    logger.info("Extracting sample IDs...")
    sample_ids = get_sample_ids(args.ids_file, args.header)

    logger.info("Preparing to start batch job...")
    input = b.read_input(readviz_ht_path())

    for sample in sample_ids:
        j = b.new_job(name=sample)
        j.image("gcr.io/broad-mpg-gnomad/tgg-methods-vm")
        j.storage("50Gi")
        # NOTE: the number of cpus specified with `local` must match this number
        j.cpu(1)
        j.command(
            f"""
            import hail as hl \

            hl.init(local='local[1]') \ 
            ht = hl.read_table({input}) \
            ht = ht.filter(ht.s == {sample}) \
            ht = ht.naive_coalesce(1) \
            ht.export({j.ofile}) \
            """
        )

        b.write_output(j.ofile, f"{readviz_per_sample_tsv_path()}/{sample}.tsv")

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
