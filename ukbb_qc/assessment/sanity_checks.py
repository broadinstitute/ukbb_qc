import logging
import hail as hl
from typing import Tuple


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("sanity_checks")
logger.setLevel(logging.INFO)


def summarize_mt(mt: hl.MatrixTable) -> hl.Struct:
    """
    Gets a summary of variants in a MatrixTable.
    Prints number of variants to stdout, and checks that each chromosome has variant calls

    :param MatrixTable mt: Raw MatrixTable to be checked
    :return: Struct of MatrixTable variant summary
    :rtype: Struct
    """

    var_summary = hl.summarize_variants(mt, show=False)
    logger.info(f"Dataset has {var_summary.n_variants} variants")

    # check that all contigs have variant calls
    for contig in var_summary.contigs:
        if var_summary.contigs[contig] == 0:
            logger.warning(f"{contig} has no variants called")

    return var_summary


def check_adj(
    mt: hl.MatrixTable, mt_adj: hl.MatrixTable, gt_expr: hl.expr.CallExpression
) -> bool:
    """
    Checks if MatrixTable has been filtered using adj criteria by checking allele counts pre and post adj filtration

    :param MatrixTable mt: MatrixTable to be checked
    :param MatrixTable mt_adj: MatrixTable filtered using adj criteria
    :param hl.expr.CallExpression gt_expr: Field containing genotype information 
    :return: Bool of whether MatrixTable has been adj filtered
    :rtype: bool
    """

    pre = mt.aggregate_entries(hl.agg.counter(mt[gt_expr].n_alt_alleles()))
    logger.info(f"\nAllele distribution pre adj filtration: {pre}")
    post = mt_adj.aggregate_entries(hl.agg.counter(mt_adj[gt_expr].n_alt_alleles()))
    logger.info(f"\nAllele distribution post adj filtration: {post}")

    adj = False
    if sum(pre.values()) != sum(post.values()):
        adj = True

    return adj


def sample_check(
    ht: hl.Table, exp_ht: hl.Table, show_mismatch: bool = True,
) -> Tuple[bool, bool]:
    """
    Checks for sample mismatch between samples in two Tables.
    If there is a sample mismatch, writes unique samples to output txt file
    Assumes the keys of the two tables match (uses anti_join).

    :param Table ht: Table containing samples to be checked
    :param Table exp_ht: Table with one column containing expected samples
    :param bool show_mismatch: Boolean whether to print sample mismatches to stdout. Default is True
    :return: Tuple of bools [whether there were missing samples, whether there were extra samples]
    :rtype: Tuple[bool, bool]
    """
    # bool to store whether there are samples missing from ht
    missing = False
    # bool to store whether ht contains samples not in exp_ht
    extra = False

    missing_samples = exp_ht.anti_join(ht).select()
    n_missing_samples = missing_samples.count()
    extra_samples = ht.anti_join(exp_ht).select()
    n_extra_samples = extra_samples.count()

    if n_missing_samples > 0:
        missing = True
        logger.info(
            f"Total number of IDs that are not in the sample HT: {n_missing_samples}..."
        )
        if show_mismatch:
            missing_samples.show(n_missing_samples)

    if n_extra_samples > 0:
        extra = True
        logger.info(f"Total number of extra IDs in the sample HT: {n_extra_samples}...")
        if show_mismatch:
            extra_samples.show(n_extra_samples)

    return (missing, extra)