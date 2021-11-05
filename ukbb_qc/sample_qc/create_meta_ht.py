import argparse
import logging

import hail as hl

from gnomad.utils.slack import slack_notifications
from ukbb_qc.resources.basics import (
    get_checkpoint_path,
    get_ukbb_data,
    logging_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import meta_ht_path
from ukbb_qc.slack_creds import slack_token


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_meta_ht")
logger.setLevel(logging.INFO)


def main(args):

    hl.init(log="/create_meta.log", default_reference="GRCh38")

    logger.info("Setting hail flag (to try to speed up computations)...")
    hl._set_flags(no_whole_stage_codegen="1")

    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)

    try:
        logger.info("Reading in raw MT and filtering to UKBB samples only...")
        ht = get_ukbb_data(
            *tranche_data,
            split=False,
            raw=True,
            ukbb_samples_only=True,
            meta_root="meta",
        ).cols()
        ht = ht.naive_coalesce(args.n_partitions)
        logger.info("Number of samples: %i", ht.count())

        logger.info("Un-nesting sample meta information...")
        ht = ht.annotate(**ht.meta)
        ht = ht.drop("meta")

        logger.info("Keying by UKBB application 26041 ID...")
        ht = ht.transmute(ukb_meta=ht.ukbb_meta)
        ht = ht.key_by(s=ht.ukb_meta.ukbb_app_26041_id)
        ht = ht.annotate(ukb_meta=ht.ukb_meta.drop("ukbb_app_26041_id"))

        logger.info("Removing extra ancestry information...")
        ht = ht.transmute(gnomad_pc_project_pop=ht.gnomad_pc_project_pop_data.pop)
        ht = ht.transmute(hybrid_pop=ht.hybrid_pop_data.pop)
        ht = ht.transmute(pan_ancestry_pop=ht.pan_ancestry_meta.pop)

        logger.info("Getting globals from original sample meta HT...")
        meta_ht = hl.read_table(meta_ht_path(*tranche_data)).select()
        logger.info("Dropping unnecessary globals from meta HT...")
        meta_ht = meta_ht.drop("platform_inference_hdbscan_parameters")
        ht = ht.annotate_globals(**meta_ht.index_globals())

        logger.info(
            "Writing new meta HT to temp (to avoid overwriting meta HT used for QC (written to `meta_ht_path`)..."
        )
        ht = ht.checkpoint(
            get_checkpoint_path(*tranche_data, "ror_meta", mt=False),
            overwrite=args.overwrite,
        )

        logger.info("Writing HT to requester pays bucket...")
        ht.write(f"{args.requester_pays_path}/meta.ht", overwrite=args.overwrite)
        ht.export(f"{args.requester_pays_path}/meta.tsv.bgz")
        for metric in ht.row:
            ht[f"{metric}"].summarize()
        logger.info(f"Final count: {ht.count()}")
        logger.info("Complete")

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(*tranche_data))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="This script creates the sample meta table for UKBB data"
    )
    parser.add_argument(
        "--n_partitions",
        help="Desired number of partitions for output HT",
        default=5000,
        type=int,
    )
    parser.add_argument(
        "-f", "--freeze", help="Current freeze", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--requester_pays_path",
        help="Path to requester pays bucket",
        default="gs://broad-ukbb-requester-pays/broad.freeze_7/sample_metadata",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
