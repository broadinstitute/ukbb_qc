import argparse
import logging

import hail as hl

from gnomad.sample_qc.filtering import compute_stratified_metrics_filter
from gnomad.utils.annotations import bi_allelic_expr
from gnomad.utils.filtering import (
    filter_low_conf_regions,
    filter_to_autosomes,
)
from gnomad.utils.slack import try_slack
from ukbb_qc.resources.basics import (
    array_sample_map_ht_path,
    get_ukbb_data,
    logging_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import (
    ancestry_hybrid_ht_path,
    platform_pca_assignments_ht_path,
    platform_pop_outlier_ht_path,
    qc_temp_data_prefix,
)
from ukbb_qc.utils.utils import annotate_interval_qc_filter, remove_hard_filter_samples


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("outlier_filter")
logger.setLevel(logging.INFO)


def run_sample_qc(mt: hl.MatrixTable) -> hl.Table:
    """
    Filter input MatrixTable to bi-allelic sites, removes problematic intervals, and computes sample QC metrics

    :param MatrixTable mt: Input MatrixTable on which sample QC metrics need to be computed
    :return: MatrixTable filtered to autosomes and high-confidence regions, with computed sample QC column annotations
    :rtype: hl.MatrixTable
    """
    # Filtering to bi-allelic sites outside of low confidence regions
    mt = filter_to_autosomes(mt)
    mt = filter_low_conf_regions(
        mt, filter_lcr=True, filter_decoy=False, filter_segdup=False
    )
    mt = mt.filter_rows(bi_allelic_expr(mt))

    logger.info("Starting sample QC...")
    mt = hl.sample_qc(mt)
    return mt.cols().select("sample_qc")


def main(args):

    hl.init(log="/outlier_detection.log", default_reference="GRCh38")

    data_source = "broad"
    freeze = args.freeze
    pop_assignment_method = args.pop_assignment_method

    try:
        if args.run_mini_qc:
            # NOTE: we run outlier detection without adj filtration to get better separation between high and low quality samples
            # this is per Julia's discussion with Konrad in #ukbb_qc
            # Need all workers for the mini qc
            mt = get_ukbb_data(data_source, freeze, adj=False)
            mt = annotate_interval_qc_filter(
                data_source,
                freeze,
                mt,
                cov_filter_field=args.cov_filter_field,
                pct_samples=args.pct_samples,
                autosomes_only=True,
            )
            mt = mt.filter_rows(mt.interval_qc_pass)

            logger.info("Filtering samples that fail hard filters...")
            mt = remove_hard_filter_samples(data_source, freeze, mt, gt_field="GT")
            logger.info(f"Sample count post filtering: {mt.count_cols()}")

            logger.info(
                "Running mini sample QC for platform- and population-specific filtering..."
            )
            sample_qc_ht = run_sample_qc(mt)
            sample_qc_ht = sample_qc_ht.repartition(args.n_partitions)
            sample_qc_ht.write(
                qc_temp_data_prefix(data_source, freeze)
                + "outlier_sample_qc_intervals.ht",
                args.overwrite,
            )

        sample_qc_ht = hl.read_table(
            qc_temp_data_prefix(data_source, freeze) + "outlier_sample_qc_intervals.ht"
        )
        strata = {}

        logger.info("Annotating population assignments...")
        pop_ht = hl.read_table(ancestry_hybrid_ht_path(data_source, freeze))
        sample_qc_ht = sample_qc_ht.annotate(
            qc_pop=pop_ht[sample_qc_ht.key][pop_assignment_method]
        )
        logger.info("Annotating inferred platform assignments...")
        platform_ht = hl.read_table(
            platform_pca_assignments_ht_path(data_source, freeze)
        )
        sample_qc_ht = sample_qc_ht.annotate(
            qc_platform=platform_ht[sample_qc_ht.key].qc_platform
        )
        logger.info("Annotating with batch (tranche) as a proxy for platform...")
        sample_map_ht = hl.read_table(array_sample_map_ht_path(freeze))
        sample_qc_ht = sample_qc_ht.annotate(
            batch=sample_map_ht[sample_qc_ht.key].batch
        )

        if args.population_filter:
            strata["qc_pop"] = sample_qc_ht.qc_pop
        else:
            pop_assignment_method = "no_pop"

        if args.platform_filter or args.batch_filter:
            if args.platform_filter:
                strata["qc_platform"] = sample_qc_ht.qc_platform
                platform_assignment_method = "qc_platform"

            if args.batch_filter:
                # NOTE: used tranche as a proxy for platform in tranche 2/freeze 5/200K
                strata["qc_platform"] = sample_qc_ht.batch
                platform_assignment_method = "batch"
        else:
            platform_assignment_method = "no_platform"

        # Make qc_metrics a dict (needs to be dict for compute_stratified_metrics_filter)
        metrics = args.filtering_qc_metrics.split(",")
        qc_metrics = {}
        for m in metrics:
            qc_metrics[m] = sample_qc_ht.sample_qc[f"{m}"]

        # For each platform and population, aggregate sample QC metrics and calculate the MAD/mean/stdev
        logger.info(
            "Flagging samples failing pop/platform-specific sample qc thresholds..."
        )
        pop_platform_filter_ht = compute_stratified_metrics_filter(
            sample_qc_ht, qc_metrics, strata
        )
        pop_platform_filter_ht = pop_platform_filter_ht.repartition(args.n_partitions)
        pop_platform_filter_ht.write(
            platform_pop_outlier_ht_path(
                data_source, freeze, pop_assignment_method, platform_assignment_method
            ),
            overwrite=args.overwrite,
        )
        num_pass = pop_platform_filter_ht.aggregate(
            hl.agg.count_where(hl.len(pop_platform_filter_ht.qc_metrics_filters) == 0)
        )
        logger.info(
            f"{num_pass} exome samples found passing pop/platform-specific filtering"
        )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(data_source, freeze))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--n_partitions",
        help="Desired number of partitions. Used for ALL output HTs",
        default=5000,
        type=int,
    )

    parser.add_argument(
        "--run_mini_qc",
        help="Run mini sample qc needed for outlier filtering",
        action="store_true",
    )
    parser.add_argument(
        "--cov_filter_field",
        help="Coverage field to use to filter high coverage intervals",
        default="pct_samples_20x",
    )
    parser.add_argument(
        "--pct_samples",
        help="Percent samples at specified coverage to filter intervals",
        default=0.85,
        type=float,
    )
    platform_group = parser.add_mutually_exclusive_group()
    platform_group.add_argument(
        "--platform_filter",
        help="Include inferred platforms in outlier filtering",
        action="store_true",
    )
    platform_group.add_argument(
        "--batch_filter",
        help="Include batch (tranche) as proxy for platform in outlier filtering",
        action="store_true",
    )
    parser.add_argument(
        "--population_filter",
        help="Included inferred population in outlier filtering",
        action="store_true",
    )
    parser.add_argument(
        "--pop_assignment_method",
        help="Population assignment method to use for stratification",
        default="hybrid_pop",
        choices=["gnomad_pc_project_pop", "HDBSCAN_pop_cluster", "hybrid_pop"],
    )
    parser.add_argument(
        "--filtering_qc_metrics",
        help="List of QC metrics for filtering.",
        default=",".join(
            [
                "n_snp",
                "r_ti_tv",
                "r_insertion_deletion",
                "n_insertion",
                "n_deletion",
                "r_het_hom_var",
            ]
        ),
    )

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
