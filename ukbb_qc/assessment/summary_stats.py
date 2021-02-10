import argparse
import logging

import hail as hl

from gnomad.assessment.summary_stats import (
    default_generate_gene_lof_summary,
    default_generate_gene_lof_matrix,
    get_summary_counts,
)
from gnomad.utils.slack import slack_notifications

from ukbb_qc.resources.basics import (
    get_ukbb_data,
    logging_path,
    release_ht_path,
    release_summary_ht_path,
    release_lof_ht_path,
    release_lof_mt_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import meta_ht_path
from ukbb_qc.slack_creds import slack_token


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("summary_stats")
logger.setLevel(logging.INFO)


def main(args):

    hl.init(log="/summary_stats.log", default_reference="GRCh38")

    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)

    try:
        if args.get_summary_counts:
            ht = hl.read_table(release_ht_path(*tranche_data)).select_globals()
            if args.interval_qc_pass_only:
                logger.info("Removing regions that fail interval QC...")
                ht = ht.filter(~ht.region_flag.fail_interval_qc)

            logger.info("Getting summary counts per variant category...")
            ht = get_summary_counts(ht)
            meta_ht = hl.read_table(meta_ht_path(*tranche_data))
            meta_ht = meta_ht.filter(
                meta_ht.sample_filters.high_quality
                & hl.is_defined(meta_ht.ukbb_meta.batch)
                & ~(meta_ht.sample_filters.related)
            )
            logger.info(f"Number of high quality samples: {meta_ht.count()}")
            ht = ht.annotate_globals(num_high_qual_samples=meta_ht.count())
            ht.write(
                release_summary_ht_path(
                    *tranche_data, intervals=args.interval_qc_pass_only
                ),
                args.overwrite,
            )

        if args.generate_gene_lof_matrix:
            # Konrad released metrics on adj GT only (gnomAD v2.1)
            mt = get_ukbb_data(*tranche_data, adj=True, split=True, meta_root="meta")
            release_ht = hl.read_table(release_ht_path(*tranche_data))
            mt = mt.annotate_rows(
                freq=release_ht[mt.row_key].freq,
                vep=release_ht[mt.row_key].vep,
                filters=release_ht[mt.row_key].filters,
                fail_interval_qc=release_ht[mt.row_key].region_flag.fail_interval_qc,
            )

            if args.interval_qc_pass_only:
                logger.info("Removing regions that fail interval QC...")
                mt = mt.filter_rows(~mt.fail_interval_qc)

            mt = default_generate_gene_lof_matrix(mt=mt, tx_ht=None)
            mt.write(
                release_lof_mt_path(
                    *tranche_data, intervals=args.interval_qc_pass_only
                ),
                args.overwrite,
            )

        if args.summarize_gene_lof_matrix:
            mt = hl.read_matrix_table(
                release_lof_mt_path(*tranche_data, intervals=args.interval_qc_pass_only)
            )

            # Reformatting pop annotations here to call `default_generate_gene_lof_summary`
            mt = mt.annotate_cols(
                meta=mt.meta.annotate(
                    pop=mt.meta.hybrid_pop_data.pop
                    if args.use_hybrid_pop
                    else mt.meta.pan_ancestry_meta.pop
                )
            )

            ht = default_generate_gene_lof_summary(mt)
            ht.write(
                release_lof_ht_path(
                    *tranche_data, intervals=args.interval_qc_pass_only
                ),
                args.overwrite,
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(data_source, freeze))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite data (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--get_summary_counts",
        help="Get summary counts per variant category",
        action="store_true",
    )
    parser.add_argument(
        "--generate_gene_lof_matrix",
        help="Generate gene LoF matrix",
        action="store_true",
    )
    parser.add_argument(
        "--summarize_gene_lof_matrix",
        help="Creates gene LoF matrix summary Table",
        action="store_true",
    )
    parser.add_argument(
        "--use_hybrid_pop",
        help="Use hybrid ancestry assignments when generating LoF matrix summary Table. Will use pan-ancestry assignment if not set",
        action="store_true",
    )
    parser.add_argument(
        "--interval_qc_pass_only",
        help="Generate summary stats on interval QC pass regions only",
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
