import argparse
import logging

import hail as hl

from gnomad.resources.grch37.gnomad import liftover
from gnomad.resources.grch38.gnomad import public_release
from gnomad.utils.annotation import hemi_expr
from gnomad.utils.slack import slack_notifications

from ukbb_qc.resources.basics import (
    cram_map_path,
    get_ukbb_data,
    non_gnomad_var_ht_path,
    readviz_ht_path,
)
from ukbb_qc.slack_creds import slack_token


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("get_samples_for_readviz")
logger.setLevel(logging.INFO)


def get_gnomad_variants() -> hl.Table:
    """
    Joins gnomAD v2 exomes liftover and v3.1 genomes variants.

    :return: Table containing all gnomAD variants.
    :rtype: hl.Table
    """
    exomes_ht = liftover("exomes").ht().select().select_globals()
    # NOTE: Specifying 3.1 for gnomAD genomes version
    genomes_ht = (
        public_release(data_type="genomes")
        .versions["3.1"]
        .ht()
        .select()
        .select_globals()
    )
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

    logger.info("Filtering to samples with crams...")
    cram_ht = hl.import_table(cram_map_path(*tranche_data), no_header=True).key_by("f0")
    mt = mt.filter_cols(hl.is_defined(cram_ht[mt.s]))
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    if args.get_variants:

        ht = mt.rows()
        logger.info("Getting gnomAD variants...")
        gnomad_ht = get_gnomad_variants()

        logger.info("Extracting variants NOT present in gnomAD...")
        ht = ht.anti_join_rows(gnomad_ht)
        ht.write(non_gnomad_var_ht_path(*tranche_data), overwrite=args.overwrite)

    if args.get_samples:
        logger.info("Filtering to variants not present in gnomAD...")
        non_gnomad_var_ht = hl.read_table(non_gnomad_var_ht_path(*tranche_data))
        mt = mt.filter_rows(hl.is_defined(non_gnomad_var_ht[mt.row_key]))

        logger.info(
            f"Taking up to {args.num_samples} samples per site where samples are het, hom_var, or hemi"
        )
        mt = mt.annotate_rows(
            samples_w_het_var=hl.agg.filter(
                mt.GT.is_het(),
                hl.agg.take(
                    hl.struct(s=mt.s, GQ=mt.GQ), args.num_samples, ordering=-mt.GQ
                ),
            ),
            samples_w_hom_var=hl.agg.filter(
                mt.GT.is_hom_var() & mt.GT.is_diploid(),
                hl.agg.take(
                    hl.struct(s=mt.s, GQ=mt.GQ), args.num_samples, ordering=-mt.GQ
                ),
            ),
            samples_w_hemi_var=hl.agg.filter(
                hemi_expr(mt.locus, mt.meta.sex_imputation.sex_karyotype, mt.GT),
                hl.agg.take(
                    hl.struct(s=mt.s, GQ=mt.GQ), args.num_samples, ordering=-mt.GQ
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
                s: str,
                GQ: int32
            }>
            'samples_w_hom_var': array<struct {
                s: str,
                GQ: int32
            }>
            'samples_w_hemi_var': array<struct {
                s: str,
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
