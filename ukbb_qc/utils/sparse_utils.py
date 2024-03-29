import logging
from typing import List

import hail as hl

from gnomad.utils.annotations import bi_allelic_expr
from gnomad.utils.filtering import filter_to_autosomes
from ukbb_qc.resources.sample_qc import callrate_mt_path


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sparse_utils")
logger.setLevel(logging.INFO)


def compute_interval_callrate_dp_mt(
    data_source: str,
    freeze: int,
    mt: hl.MatrixTable,
    intervals_ht: hl.Table,
    bi_allelic_only: bool = True,
    autosomes_only: bool = True,
    target_pct_gt_cov: List = [10, 20],
    overwrite: bool = False,
) -> None:
    """
    Computes sample metrics (n_defined, total, mean_dp, pct_gt_20x, pct_dp_defined) per interval. 
    
    Mean depth, pct_gt_{target_pct_gt_cov}x, and pct_dp_defined annotations are used during interval QC.
    Mean depth and callrate annotations (mean_dp, n_defined, total) are used during hard filtering.
    Callrate annotations (n_defined, total) are also used during platform PCA.
    Writes call rate mt (aggregated MatrixTable) keyed by intervals row-wise and samples column-wise.
    Assumes that all overlapping intervals in intervals_ht are merged.
    NOTE: This function requires a densify! Please use an autoscaling cluster.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param MatrixTable mt: Input MatrixTable.
    :param Table intervals_ht: Table with capture intervals relevant to input MatrixTable.
    :param bool bi_allelic_only: If set, only bi-allelic sites are used for the computation.
    :param bool autosomes_only: If set, only autosomal intervals are used.
    :param List target_pct_gt_cov: Coverage levels to check for each target. Default is [10, 20].
    :param bool overwrite: Whether to overwrite existing data. Default is False.
    :return: None
    """
    logger.warning(
        "This function will densify! Make sure you have an autoscaling cluster."
    )
    logger.info("Computing call rate and mean DP MatrixTable...")
    if len(intervals_ht.key) != 1 or not isinstance(
        intervals_ht.key[0], hl.expr.IntervalExpression
    ):
        logger.warning(
            f"Call rate matrix computation expects `intervals_ht` with a key of type Interval. Found: {intervals_ht.key}"
        )
    logger.info("Densifying...")
    mt = mt.drop("gvcf_info")
    mt = hl.experimental.densify(mt)

    logger.info(
        "Filtering out lines that are only reference or not covered in capture intervals"
    )
    logger.warning(
        "Assumes that overlapping intervals in capture intervals are merged!"
    )
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)
    intervals_ht = intervals_ht.annotate(interval_label=intervals_ht.interval)
    mt = mt.annotate_rows(interval=intervals_ht.index(mt.locus).interval_label)
    mt = mt.filter_rows(hl.is_defined(mt.interval))

    if autosomes_only:
        mt = filter_to_autosomes(mt)

    if bi_allelic_only:
        mt = mt.filter_rows(bi_allelic_expr(mt))

    logger.info(
        "Grouping MT by interval and calculating n_defined, total, and mean_dp..."
    )
    mt = mt.annotate_entries(GT=hl.experimental.lgt_to_gt(mt.LGT, mt.LA))
    mt = mt.select_entries(
        GT=hl.or_missing(hl.is_defined(mt.GT), hl.struct()),
        DP=hl.if_else(hl.is_defined(mt.DP), mt.DP, 0),
    )

    # Lines 87-88 should eliminate the shuffle
    mt = mt.key_rows_by(mt.interval)
    mt = mt._filter_partitions(
        [], keep=False
    )  # hack to prevent optimization back to the original execution
    mt = (
        mt.group_rows_by(mt.interval)
        .aggregate_rows(n_var=hl.agg.count())
        .aggregate(
            n_defined=hl.agg.count_where(hl.is_defined(mt.GT)),
            dp_sum=hl.agg.sum(mt.DP),
            mean_dp=hl.agg.mean(mt.DP),
            **{
                f"pct_gt_{cov}x": hl.agg.fraction(mt.DP >= cov)
                for cov in target_pct_gt_cov
            },
            pct_dp_defined=hl.agg.count_where(mt.DP > 0) / hl.agg.count(),
        )
    )
    mt.write(
        callrate_mt_path(data_source, freeze, interval_filtered=False),
        overwrite=overwrite,
    )
