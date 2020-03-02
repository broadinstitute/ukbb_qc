import argparse
import hail as hl
import logging
from gnomad_hail.utils.sample_qc import add_filters_expr
from gnomad_hail.utils.sparse_mt import densify_sites


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
    min_callrate: float, 
    min_depth: float
) -> hl.expr.SetExpression:
    """
    Creates hard filters expression.

    :param hl.expr.Float64Expression callrate_expr: Field that contains callrate expression
    :param hl.expr.Float64Expression dp_expr: Field that contains mean depth expression
    :param hl.expr.StringExpression sex_expr: Field that contains inferred sex expression
    :param float min_callrate: Callrate threshold to be used to filter samples; default is 0.99
    :param float min_depth: Mean depth threshold to be used to filter samples; default is 20.0
    :return: Hard filters expression
    :rtype: hl.expr.SetExpression
    """

    # the default coverage/depth cutoffs were set visually using plots:
    # p = hl.plot.histogram(mt.sample_qc.dp_stats.mean, range=(10,120), legend="Mean Sample DP")
    # p = hl.plot.histogram(mt.sample_qc.call_rate, range=(0.991, 0.997), legend="Mean Sample Callrate")

    logger.info("Callrate cutoff for hard filters: {}".format(min_callrate))
    logger.info("Depth cutoff for hard filters: {}".format(min_depth))

    hard_filters = {
        # we don"t have contamination/chimera for regeneron vcf
        # "contamination": ht.freemix > 0.05,
        # "chimera": ht.pct_chimeras > 0.05,
        "low_callrate": callrate_expr < min_callrate,
        "ambiguous_sex": sex_expr == "Ambiguous",
        "sex_aneuploidy": ((sex_expr != "Ambiguous") & (sex_expr != "XX") & (sex_expr != "XY")),
        "low_coverage": dp_expr < min_depth,
    }
    return hard_filters


def hard_filter_samples(
    data_source: str,
    freeze: int,
    mt: hl.MatrixTable,
    interval_qc_ht: hl.Table,
    sex_ht: hl.Table,
    min_callrate: float,
    min_depth: float
) -> hl.Table:
    """
    Applys hard filters to samples and returns Table with samples and their hard filter status.

    :param data_source: One of 'regeneron' or 'broad'
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
    ht = mt.annotate_cols(
        call_rate=hl.agg.sum(mt.n_defined) / hl.agg.sum(mt.total),
        mean_dp=hl.agg.sum(mt.dp_sum) / hl.agg.sum(mt.total)
    ).cols()

    logger.info("Adding sex imputation annotations...")
    ht = ht.annotate(sex=sex_ht[ht.key].sex_karyotype)

    ht = ht.checkpoint(
        get_checkpoint_path(data_source, freeze, name="interval_qc_sample_qc"),
        overwrite=True,
    )

    logger.info("Applying hard filters and writing out hard filters HT...")
    ht = ht.annotate(apply_hard_filters_expr(ht, min_callrate, min_depth))
    ht = ht.annotate(
        hard_filters=add_filters_expr(
            apply_hard_filters_expr(
                    ht.call_rate,
                    ht.mean_dp,
                    ht.sex,
                    min_callrate,
                    min_depth
            ),
            None
        )
    )
    return ht
