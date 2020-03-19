import argparse
import logging
import hail as hl
from gnomad.utils.slack import try_slack
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE, DataException
from ukbb_qc.resources.basics import capture_ht_path, get_checkpoint_path, get_ukbb_data
from ukbb_qc.resources.sample_qc import (
    callrate_mt_path,
    hard_filters_ht_path,
    interval_qc_path,
    sex_ht_path,
)
from ukbb_qc.utils.sparse_utils import compute_interval_callrate_dp_mt


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("apply_hard_filters")
logger.setLevel(logging.INFO)


def apply_hard_filters_expr(
    callrate_expr: hl.expr.Float64Expression,
    dp_expr: hl.expr.Float64Expression,
    sex_expr: hl.expr.StringExpression,
    min_callrate: float = 0.99,
    min_depth: float = 20.0,
) -> hl.expr.StructExpression:
    """
    Creates hard filters expression.

    :param hl.expr.Float64Expression callrate_expr: Field that contains callrate expression
    :param hl.expr.Float64Expression dp_expr: Field that contains mean depth expression
    :param hl.expr.StringExpression sex_expr: Field that contains inferred sex expression
    :param float min_callrate: Callrate threshold to be used to filter samples; default is 0.99
    :param float min_depth: Mean depth threshold to be used to filter samples; default is 20.0
    :return: Hard filters expression
    :rtype: hl.expr.StructExpression
    """

    # The default coverage/depth cutoffs were set (tranche 2/freeze 5) visually using plots:
    # p = hl.plot.histogram(mt.sample_qc.dp_stats.mean, range=(10,120), legend="Mean Sample DP")
    # p = hl.plot.histogram(mt.sample_qc.call_rate, range=(0.991, 0.997), legend="Mean Sample Callrate")

    logger.info(f"Callrate cutoff for hard filters: {min_callrate}")
    logger.info(f"Depth cutoff for hard filters: {min_depth}")

    hard_filters = hl.struct(
        # we don't have contamination/chimera for broad mt
        # "contamination": ht.freemix > 0.05,
        # "chimera": ht.pct_chimeras > 0.05,
        low_callrate=callrate_expr < min_callrate,
        ambiguous_sex=sex_expr == "Ambiguous",
        sex_aneuploidy=(
            (sex_expr != "Ambiguous") & (sex_expr != "XX") & (sex_expr != "XY")
        ),
        low_coverage=dp_expr < min_depth,
    )
    return hard_filters


def hard_filter_samples(
    data_source: str,
    freeze: int,
    mt: hl.MatrixTable,
    interval_qc_ht: hl.Table,
    sex_ht: hl.Table,
    min_callrate: float = 0.99,
    min_depth: float = 20.0,
) -> hl.Table:
    """
    Applys hard filters to samples and returns Table with samples and their hard filter status.

    This function expects the input MT to be annotated with the fields n_defined, total, and dp_sum.
    These are calculated using compute_interval_callrate_dp_mt.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param MatrixTable mt: Input MatrixTable with samples to be filtered
    :param Table interval_qc_ht: Table with high coverage intervals
    :param Table sex_ht: Table with samples and their inferred sex
    :param float min_callrate: Callrate threshold to be used to filter samples; default is 0.99
    :param float min_depth: Mean depth threshold to be used to filter samples; default is 20.0
    :return: Table with samples and their hard filter status
    :rtype: hl.Table
    """
    logger.info("Computing callrate and mean DP over high coverage intervals...")
    mt = mt.filter_rows(hl.is_defined(interval_qc_ht[mt.row_key]))
    mt = mt.checkpoint(callrate_mt_path(data_source, freeze, interval_filtered=True))
    ht = mt.annotate_cols(
        call_rate=hl.agg.sum(mt.n_defined) / hl.agg.sum(mt.total),
        mean_dp=hl.agg.sum(mt.dp_sum) / hl.agg.sum(mt.total),
    ).cols()

    logger.info("Adding sex imputation annotations...")
    ht = ht.annotate(sex=sex_ht[ht.key].sex_karyotype)

    ht = ht.checkpoint(
        get_checkpoint_path(data_source, freeze, name="interval_qc_sample_qc"),
        overwrite=True,
    )

    logger.info("Applying hard filters and writing out hard filters HT...")
    ht = ht.annotate(
        hard_filters=apply_hard_filters_expr(
            ht.call_rate, ht.mean_dp, ht.sex, min_callrate, min_depth
        )
    )
    ht = ht.annotate(
        ht.hard_filters.annotate(
            hard_filtered=(
                (ht.hard_filters.low_callrate)
                | (ht.hard_filters.ambiguous_sex)
                | (ht.hard_filters.sex_aneuploidy)
                | (ht.hard_filters.low_coverage)
            )
        )
    )
    return ht.drop("sex")


def main(args):

    hl.init(log="/apply_hard_filters.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze

    # NOTE: this function will densify
    if args.compute_interval_callrate_mt:
        logger.info("Reading in raw MT...")
        mt = get_ukbb_data(
            data_source, freeze, split=False, raw=True, key_by_locus_and_alleles=True,
        )
        logger.info(
            f"Total number of variants in raw unsplit matrix table: {mt.count_rows()}"
        )
        capture_ht = hl.read_table(capture_ht_path(data_source))
        compute_interval_callrate_dp_mt(data_source, freeze, mt, capture_ht)

    if args.apply_hard_filters:
        if not file_exists(
            callrate_mt_path(data_source, freeze, interval_filtered=False)
        ):
            raise DataException("Need to compute interval callrate MT!")

        logger.info("Reading in callrate MT, sex ht, interval qc HT...")
        callrate_mt = hl.read_matrix_table(
            callrate_mt_path(data_source, freeze, interval_filtered=False)
        )
        sex_ht = hl.read_table(sex_ht_path(data_source, freeze))
        interval_qc_ht = hl.read_table(
            interval_qc_path(data_source, freeze, "autosomes")
        )
        interval_qc_ht = interval_qc_ht.filter(
            interval_qc_ht[args.cov_filter_field] > args.pct_samples
        )

        logger.info("Hard filtering samples...")
        hard_filters_ht = hard_filter_samples(
            data_source,
            freeze,
            callrate_mt,
            interval_qc_ht,
            sex_ht,
            args.min_callrate,
            args.min_dp,
        )
        ht = ht.naive_coalesce(args.n_partitions)
        ht = ht.checkpoint(
            hard_filters_ht_path(data_source, freeze), overwrite=args.overwrite,
        )
        logger.info("Checking number of samples flagged with hard filters...")
        ht = ht.explode(ht.hard_filters)
        filters = ht.aggregate(hl.agg.counter(ht.hard_filters))
        for filt in filters:
            logger.info(f"Samples flagged due to {filt}: {filters[filt]}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int,
    )
    parser.add_argument(
        "--compute_interval_callrate_mt",
        help="Computes an interval by sample mt of callrate and depth that will be for hard filtering samples",
        action="store_true",
    )
    parser.add_argument(
        "--min_callrate", help="Minimum variant callrate", default=0.99, type=float,
    )
    parser.add_argument(
        "--min_dp", help="Minimum depth", default=20.0, type=float,
    )
    parser.add_argument(
        "--apply_hard_filters",
        help="Applies hard filters to samples in MatrixTable",
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
    parser.add_argument(
        "--n_partitions",
        help="Desired number of partitions for output Table/MatrixTable",
        default=5000,
        type=int,
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to.",
    )
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
