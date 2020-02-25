import argparse
import hail as hl
import hdbscan
import logging
import numpy as np
from gnomad_hail.utils.generic import filter_to_autosomes
from gnomad_hail.utils.sample_qc import (
    assign_platform_from_pcs,
    compute_callrate_mt,
    run_platform_pca,
)
from gnomad_hail.utils.slack import try_slack
from ukbb_qc.resources.basics import capture_ht_path, get_ukbb_data
from ukbb_qc.resources.sample_qc import (
    callrate_mt_path,
    platform_pca_loadings_ht_path,
    platform_pca_results_ht_path,
    platform_pca_scores_ht_path,
)
from ukbb_qc.utils.utils import remove_hard_filter_samples


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("platform_pca")
logger.setLevel(logging.INFO)


def main(args):
    hl.init(log="/platform_pca.log", default_reference="GRCh38")

    data_source = args.data_source
    freeze = args.freeze

    if args.compute_callrate_mt:
        logger.info("Preparing data for platform PCA...")

        # Densify if not using interval qc filtered mt
        if not args.apply_interval_qc_filter:
            intervals = hl.read_table(capture_ht_path(data_source, freeze))
            mt = get_ukbb_data(data_source, freeze, adj=True)
            logger.info(f"Input MatrixTable (hardcalls) count: {mt.count()}")

            logger.info("Densifying sites...")
            mt = densify_sites(
                mt,
                intervals,
                hl.read_table(last_END_positions_ht_path(data_source, freeze)),
            )

        else:
            logger.info(
                "Reading in densified MT (created when running hard filters)..."
            )
            intervals = hl.read_table(
                interval_qc_path(data_source, freeze, "autosomes")
            )
            mt = hl.read_matrix_table(
                densified_interval_qc_path(data_source, freeze, repartitioned=True)
            )
            logger.info(f"Input MatrixTable count: {mt.count()}")

        logger.info("Removing hard filtered samples...")
        mt = remove_hard_filter_samples(data_source, freeze, mt)
        logger.info(f"Count after removing hard filtered samples: {mt.count()}")
        callrate_mt = compute_callrate_mt(mt, intervals)
        callrate_mt = callrate_mt.checkpoint(
            callrate_mt_path(
                data_source, freeze, interval_filtered=args.apply_interval_qc_filter
            ),
            overwrite=args.overwrite,
        )
        logger.info(f"Callrate MatrixTable count: {callrate_mt.count}")

    if args.run_platform_pca:
        logger.info("Running platform PCA...")
        callrate_mt = hl.read_matrix_table(
            callrate_mt_path(
                data_source, freeze, interval_filtered=args.apply_interval_qc_filter
            )
        )
        # NOTE: added None binarization_threshold parameter to make sure we things the same way as before parameter existed
        eigenvalues, scores_ht, loadings_ht = run_platform_pca(callrate_mt, None)
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
        # Regeneron freeze 4 Eigenvalues: [26489244.935849957, 2039950.6985898241, 1407875.3058482022, 1082106.1507608977, 373810.0800184624, 361301.2291929654, 324435.7483132424, 205912.4810229146, 196196.71017912056, 159808.25367132248]

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
        platform_ht = platform_ht.checkpoint(
            platform_pca_results_ht_path(
                data_source, freeze, interval_filtered=args.apply_interval_qc_filter
            ),
            overwrite=args.overwrite,
        )
        logger.info(f"Platform PCA Table count: {platform_ht.count()}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-s",
        "--data_source",
        help="Data source",
        choices=["regeneron", "broad"],
        default="broad",
    )
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )

    parser.add_argument(
        "--compute_callrate_mt",
        help="Computes an interval by sample mt of callrate that will be used to compute platform PCs",
        action="store_true",
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
