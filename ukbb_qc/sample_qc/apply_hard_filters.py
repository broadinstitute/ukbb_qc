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
    ht: hl.Table, min_callrate: float, min_depth: float
) -> hl.Table:
    """
    Creates hard filters expression and annotates Table with expression (creates hard_filters column).

    :param Table ht: Table to be annotated 
    :param float min_callrate: Callrate threshold to be used to filter samples; default is 0.99
    :param float min_depth: Mean depth threshold to be used to filter samples; default is 20.0
    :return: Table with hard_filters column
    :rtype: Table
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
        "low_callrate": ht.raw_sample_qc.call_rate < min_callrate,
        "ambiguous_sex": ht.sex == "ambiguous_sex",
        "sex_aneuploidy": ht.sex == "sex_aneuploidy",
        "low_coverage": ht.raw_sample_qc.dp_stats.mean < min_depth,
    }

    ht = ht.annotate(hard_filters=add_filters_expr(hard_filters, None))
    return ht


def hard_filter_samples(
    data_source: str,
    freeze: int,
    mt: hl.MatrixTable,
    interval_qc_ht: hl.Table,
    last_END_ht: hl.Table,
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
    :param Table last_END_ht: Table with most upstream reference block annotations
    :param Table sex_ht: Table with samples and their inferred sex
    :param float min_callrate: Callrate threshold to be used to filter samples; default is 0.99
    :param float min_depth: Mean depth threshold to be used to filter samples; default is 20.0
    :return: Table with samples and their hard filter status
    :rtype: hl.Table
    """

    logger.info("Densifying sites...")
    mt = densify_sites(mt, interval_qc_ht, last_END_ht)

    logger.info("Adding sex imputation annotations...")
    mt = mt.annotate_cols(**sex_ht[mt.col_key])

    logger.info("Computing raw sample QC metrics...")
    mt = mt.key_rows_by("locus", "alleles")
    mt = hl.sample_qc(mt)
    ht = mt.transmute_cols(
        raw_sample_qc=mt.sample_qc.select("call_rate", "gq_stats", "dp_stats")
    ).cols()
    ht = ht.checkpoint(
        get_checkpoint_path(data_source, freeze, name="interval_qc_sample_qc"),
        overwrite=True,
    )

    logger.info("Applying hard filters and writing out hard filters HT...")
    ht = apply_hard_filters_expr(ht, min_callrate, min_depth)
    return ht
