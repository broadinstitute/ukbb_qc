import argparse
import logging

import hail as hl

from gnomad.utils.slack import slack_notifications

from ukbb_qc.resources.basics import (
    get_checkpoint_path,
    logging_path,
    release_ht_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.slack_creds import slack_token


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("data_release")
logger.setLevel(logging.INFO)


def main(args):
    """
    Overwrite release HT with frequencies recalculated after fixing bugs associated with homalt hotfix.

    Run after `prepare_vcf_data_release.py` and frequencies are recalculated.
    """

    hl.init(log="/update_release_ht.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)

    try:
        logger.info("Reading in old release HT...")
        # Overwrote HT at release HT path with frequencies recalculated from first patch release,
        # which is why this HT is in temp
        ht = hl.read_table(
            "gs://broad-ukbb/broad.freeze_7/temp/broad.freeze_7.release.sites.ht"
        )

        logger.info("Reading in frequencies recalculated on homalt hotfix sites...")
        freq_ht = hl.read_table(
            get_checkpoint_path(*tranche_data, name="ukb_freq_release_patch")
        )

        logger.info(
            "Annotating new frequencies at homalt hotfix sites onto release HT..."
        )
        ht = ht.annotate(new_freq_info=hl.struct(**freq_ht[ht.key]))
        ht = ht.transmute(
            freq=hl.coalesce(ht.new_freq_info.freq, ht.freq),
            faf=hl.coalesce(ht.new_freq_info.faf, ht.faf),
            popmax=hl.coalesce(ht.new_freq_info.popmax, ht.popmax),
        )

        ht = ht.checkpoint(release_ht_path(*tranche_data), args.overwrite)
        logger.info(f"Final HT count: {ht.count()}")

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(*tranche_data))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
