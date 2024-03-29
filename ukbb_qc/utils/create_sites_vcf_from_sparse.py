import argparse
import logging
import hail as hl

from gnomad.utils.vcf import ht_to_vcf_mt
from gnomad.utils.slack import slack_notifications
from gnomad.utils.sparse_mt import default_compute_info
from ukbb_qc.resources.basics import get_ukbb_data, logging_path, vqsr_sites_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.slack_creds import slack_token


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_sites_vcf")
logger.setLevel(logging.INFO)


def main(args):

    # NOTE: This script does NOT need a densify
    # NOTE: This script was not used for freeze 7/450K. We instead used compute_info_ht.py to allow for exclusion
    # of samples filtered by Sample QC
    hl.init(log="/create_sites_vcf.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze

    try:
        logger.info("Reading in input MT (raw sparse MT)...")
        mt = get_ukbb_data(
            data_source,
            freeze,
            key_by_locus_and_alleles=True,
            split=False,
            raw=True,
            repartition=args.repartition,
            n_partitions=args.raw_partitions,
        )
        logger.info("Filtering out ref blocks")
        mt = mt.filter_rows(hl.len(mt.alleles) > 1)

        if args.test:
            logger.info("Trimming to chr20...")
            mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chr20")])
            logger.info(f"MT count after filtration: {mt.count()}")

        logger.info("Getting info HT...")
        info_ht = default_compute_info(
            mt, site_annotations=True, n_partitions=args.partitions
        )

        info_ht = info_ht.checkpoint(
            f"{vqsr_sites_path(data_source, freeze)}", overwrite=args.overwrite,
        )
        vcf_mt = ht_to_vcf_mt(info_ht)
        hl.export_vcf(
            vcf_mt,
            f"{vqsr_sites_path(data_source, freeze, ht=False)}",
            parallel="separate_header",
        )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(data_source, freeze))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--overwrite", help="Overwrite data", action="store_true",
    )
    parser.add_argument(
        "--repartition",
        help="Repartition raw MT on read. Required for 300K",
        action="store_true",
    )
    parser.add_argument(
        "--raw_partitions",
        help="Number of desired partitions for the raw MT. Necessary only for 300K. Used only if --repartition is also specified",
        default=30000,
        type=int,
    )
    parser.add_argument(
        "--partitions", help="Number of partitions for info HT", default=5000, type=int
    )
    parser.add_argument("--test", help="Filter to chr20 to test", action="store_true")
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
