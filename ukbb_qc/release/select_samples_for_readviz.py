import argparse
import logging

import hail as hl

from gnomad.resources.grch37.gnomad import liftover
from gnomad.resources.grch38.gnomad import public_release
from gnomad.utils.slack import slack_notifications

from ukbb_qc.resources.basics import (
    get_ukbb_data,
    readviz_ht_path,
    unique_variants_ht_path,
)
from ukbb_qc.slack_creds import slack_token


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("get_samples_for_readviz")
logger.setLevel(logging.INFO)


# Not sure if these should be put into some common repo?
def het_hom_hemi_take_expr(mt: hl.MatrixTable) -> hl.expr.StructExpression:
    """
    Extracts sample ID and genotype quality from input MatrixTable.

    :param hl.MatrixTable mt: Input MatrixTable.
    :return: StructExpression of sample ID and GQ.
    :rtype: hl.expr.StructExpression
    """
    return hl.struct(S=mt.s, GQ=mt.GQ)


def het_expr(mt: hl.MatrixTable) -> hl.expr.BooleanExpression:
    """
    Returns whether genotypes are heterozygous.

    :param hl.MatrixTable mt: Input MatrixTable.
    :return: BooleanExpression indicating whether genotypes are het.
    :rtype: hl.expr.BooleanExpression
    """
    return mt.GT.is_het()


def hom_expr(mt: hl.MatrixTable) -> hl.expr.BooleanExpression:
    """
    Returns whether genotypes are diploid and homozygous variant.

    :param hl.MatrixTable mt: Input MatrixTable.
    :return: BooleanExpression indicating whether genotypes are diploid and homvar.
    :rtype: hl.expr.BooleanExpression
    """
    return mt.GT.is_diploid() & mt.GT.is_hom_var()


def hemi_expr(mt: hl.MatrixTable) -> hl.expr.BooleanExpression:
    """
    Returns whether genotypes are hemizygous.

    Returns missing expression if locus is not in chrX/chrY non-PAR regions.

    :param hl.MatrixTable mt: Input MatrixTable.
    :return: BooleanExpression indicating whether genotypes are hemizygous.
    :rtype: hl.expr.BooleanExpression
    """
    return hl.or_missing(
        mt.locus.in_x_nonpar() | mt.locus.in_y_nonpar(),
        mt.GT.is_haploid() & (mt.meta.sex == "male") & (mt.GT[0] == 1),
    )


def get_gnomad_variants() -> hl.Table:
    """
    Joins gnomAD v2 exomes liftover and v3.1 genomes variants.

    :return: Table containing all gnomAD variants.
    :rtype: hl.Table
    """
    exomes_ht = liftover("exomes").ht().select().select_globals()
    genomes_ht = public_release("genomes").ht().select().select_globals()
    return exomes_ht.join(genomes_ht, how="outer")


def main(args):

    hl.init(log="/select_samples_for_readviz.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = 7
    tranche_data = (data_source, freeze)

    logger.info("Reading in hardcalls MT (filtered to adj)...")
    # NOTE: sex ploidies have already been adjusted in hardcalls
    # Keeping only necessary entries
    mt = get_ukbb_data(*tranche_data, adj=True, meta_root="meta").select_entries(
        "GT", "GQ", "DP", "AD"
    )

    # NOTE: Code is assuming here that we have access to the crams for all release samples
    logger.info("Filtering to unrelated high quality samples and their variants...")
    # Need to add filter for related here because
    # `release` field just means sample has True for both `high_quality` bool
    # and has a defined UKBB batch number
    mt = mt.filter_cols(
        mt.meta.sample_filters.release & ~mt.meta.sample_filters.related
    )
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    if args.get_variants:

        ht = mt.rows()
        logger.info("Getting gnomAD variants...")
        gnomad_ht = get_gnomad_variants()

        logger.info("Filtering to variants NOT present in gnomAD...")
        ht = ht.anti_join_rows(gnomad_ht)
        ht.write(unique_variants_ht_path(*tranche_data), overwrite=args.overwrite)

    if args.get_samples:
        logger.info(
            f"Taking up to {args.num_samples} samples per site where samples are het, hom_var, or hemi"
        )
        mt = mt.annotate_rows(
            samples_w_het_var=hl.agg.filter(
                het_expr(mt),
                hl.agg.take(
                    het_hom_hemi_take_expr(mt), args.num_samples, ordering=-mt.GQ
                ),
            ),
            samples_w_hom_var=hl.agg.filter(
                hom_expr(mt),
                hl.agg.take(
                    het_hom_hemi_take_expr(mt), args.num_samples, ordering=-mt.GQ
                ),
            ),
            samples_w_hemi_var=hl.agg.filter(
                hemi_expr(mt),
                hl.agg.take(
                    het_hom_hemi_take_expr(mt), args.num_samples, ordering=-mt.GQ
                ),
            ),
        )

        ht = mt.rows()
        ht = ht.select(
            ht.samples_w_het_var, ht.samples_w_hom_var, ht.samples_w_hemi_var
        )
        ht.write(readviz_ht_path(*tranche_data), overwrite=args.overwrite)
        """
        Output schema looks like:
        ----------------------------------------
        Global fields:
            None
        ----------------------------------------
        Row fields:
            'locus': locus<GRCh38>
            'alleles': array<str>
            'samples_w_het_var': array<struct {
                S: str,
                GQ: int32
            }>
            'samples_w_hom_var': array<struct {
                S: str,
                GQ: int32
            }>
            'samples_w_hemi_var': array<struct {
                S: str,
                GQ: int32
            }>
        ----------------------------------------
        Key: ['locus', 'alleles']
        ----------------------------------------
        """


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--get_variants",
        help="Get variants unique to UKBB (not in gnomAD)",
        action="store_true",
    )
    parser.add_argument(
        "--get_samples", help="Extract samples for readviz", action="store_true",
    )
    parser.add_argument(
        "--overwrite", help="Overwrite output data", action="store_true",
    )
    parser.add_argument(
        "--num-samples",
        type=int,
        help="Number of samples to take from each genotype category at each site",
        default=10,
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
