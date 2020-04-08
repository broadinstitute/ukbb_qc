import argparse
import logging

import hail as hl
import hdbscan

from gnomad.sample_qc.platform import (
    assign_platform_from_pcs,
    run_platform_pca,
)
from gnomad.utils.slack import try_slack
from ukbb_qc.resources.basics import logging_path
from ukbb_qc.resources.sample_qc import (
    callrate_mt_path,
    platform_pca_loadings_ht_path,
    platform_pca_assignments_ht_path,
    platform_pca_scores_ht_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.utils.utils import remove_hard_filter_samples


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("platform_pca")
logger.setLevel(logging.INFO)


def main(args):
    hl.init(log="/platform_pca.log", default_reference="GRCh38")

    data_source = "broad"
    freeze = args.freeze

    try:
        if args.run_platform_pca:
            logger.info("Running platform PCA...")
            callrate_mt = hl.read_matrix_table(
                callrate_mt_path(
                    data_source, freeze, interval_filtered=args.apply_interval_qc_filter
                )
            )

            logger.info("Removing hard filtered samples...")
            callrate_mt = remove_hard_filter_samples(
                data_source, freeze, callrate_mt, filter_rows=False
            )
            logger.info(
                f"Count after removing hard filtered samples: {callrate_mt.count_cols()}"
            )

            # NOTE: added None binarization_threshold parameter to make sure we things the same way as before parameter existed
            eigenvalues, scores_ht, loadings_ht = run_platform_pca(
                callrate_mt, binarization_threshold=None
            )
            scores_ht.write(
                platform_pca_scores_ht_path(
                    data_source, freeze, interval_filtered=args.apply_interval_qc_filter
                ),
                overwrite=args.overwrite,
            )
            logger.info(f"Scores Table count: {scores_ht.count()}")
            loadings_ht.write(
                platform_pca_loadings_ht_path(
                    data_source, freeze, interval_filtered=args.apply_interval_qc_filter
                ),
                overwrite=args.overwrite,
            )

        if args.assign_platforms:
            logger.info("Assigning platforms based on platform PCA clustering")
            scores_ht = hl.read_table(
                platform_pca_scores_ht_path(
                    data_source, freeze, interval_filtered=args.apply_interval_qc_filter
                )
            )
            platform_ht = assign_platform_from_pcs(
                scores_ht,
                hdbscan_min_cluster_size=args.hdbscan_min_cluster_size,
                hdbscan_min_samples=args.hdbscan_min_samples,
            )
            platform_ht = platform_ht.annotate_globals(
                hdbscan_min_cluster_size=args.hdbscan_min_cluster_size,
                hdbscan_min_samples=args.hdbscan_min_samples,
                interval_filtered=args.apply_interval_qc_filter,
            )
            platform_ht = platform_ht.checkpoint(
                platform_pca_assignments_ht_path(
                    data_source, freeze, interval_filtered=args.apply_interval_qc_filter
                ),
                overwrite=args.overwrite,
            )
            logger.info(f"Platform PCA Table count: {platform_ht.count()}")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(logging_path(data_source, freeze))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--apply_interval_qc_filter",  # NOTE not super sure this is necessary
        help="Filter to good intervals from interval QC",
        action="store_true",
    )
    parser.add_argument(
        "--run_platform_pca",
        help="Runs platform PCA (assumes callrate MT was computed)",
        action="store_true",
    )
    parser.add_argument(
        "--assign_platforms",
        help="Assigns platforms based on callrate PCA results using HDBSCAN",
        action="store_true",
    )
    parser.add_argument(
        "--hdbscan_min_samples",
        help="Minimum samples parameter for HDBSCAN. If not specified, --hdbscan_min_cluster_size is used.",
        type=int,
        required=False,
    )
    parser.add_argument(
        "--hdbscan_min_cluster_size",
        help="Minimum cluster size parameter for HDBSCAN.",
        type=int,
        default=100,
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
