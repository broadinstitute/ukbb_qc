import argparse
import logging
from typing import Tuple

import hail as hl

from gnomad.resources.grch37.gnomad import liftover
from gnomad.resources.grch38.gnomad import public_release
from gnomad.utils.annotations import hemi_expr
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


def get_all_gnomad_variants() -> hl.Table:
    """
    Join gnomAD v2 exomes liftover and v3.1 genomes variants.

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


def get_additional_gnomad_variants(data_type: str, input_tsv_path: str) -> hl.Table:
    """
    Get additional gnomAD variants with fewer than the maximum number of samples displayed.

    If getting gnomAD v2.1 exome variants, lift variants from GRCh37 to GRCh38.
    """

    def _liftover_gnomad_exomes_variants(ht: hl.Table) -> hl.Table:
        """
        Lift gnomAD v2.1 exomes variants to GRCh38 using public liftover table annotations.

        :param hl.Table ht: Input table of gnomAD v2.1 exomes variants.
        :return: gnomAD v2.1 exomes variants lifted over from GRCh37 to GRCh38.
        :rtype: hl.Table
        """
        # Not sure where this resource is/will be stored in gnomad methods/qc
        # Liftover table re-created with original alleles
        # This table has the following fields:
        """
        'locus': locus<grch37>
        'alleles': array<str>
        'liftover_proc_id': struct {
            'part_idx': int32
            'block_idx': int32
        }
        'liftover_locus': locus<grch38>
        'liftover_alleles': array<str>
        Key: ['locus', 'alleles']
        """
        ht = ht.annotate(locus=hl.locus(ht.f0, ht.f1, reference_genome="GRCh37"))
        ht = ht.key_by("locus", "alleles")
        exomes_ht = hl.read_table(
            "gs://gnomad-browser/gnomad-liftover/output.ht"
        ).select("liftover_locus", "liftover_alleles")
        ht = ht.annotate(**exomes_ht[ht.key])
        ht = ht.key_by()
        return ht.key_by(locus=ht.liftover_locus, alleles=ht.liftover_alleles)

    ht = hl.import_table(input_tsv_path, no_header=True)
    if data_type == "exomes":
        ht = _liftover_gnomad_exomes_variants(ht)
        ht = ht.transmute(
            alleles=[ht.f2, ht.f3], het_or_hom_or_hemi=ht.f4, n_available_samples=ht.f5,
        )
        # Drop extra field in TSV (extra count taken with sqlite query)
        return ht.drop("f6")
    else:
        ht = ht.annotate(locus=hl.locus(ht.f0, ht.f1))
        return ht.transmute(
            alleles=[ht.f2, ht.f3], zygosity=ht.f4, n_available_samples=ht.f5,
        )


def main(args):

    hl.init(log="/select_samples_for_readviz.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = 7
    tranche_data = (data_source, freeze)

    def _sample_ordering_expr(
        mt: hl.MatrixTable,
    ) -> Tuple[hl.expr.Int32Expression, hl.expr.Float64Expression]:
        """
            Generate a random number to be used along with genotype quality for ordering sample IDs.

            It can be problematic for downstream steps in the readviz pipeline when 
            several samples have many times more variants selected than in other samples. 
            To avoid this, and distribute variants more evenly across samples, add a random number as the secondary sort order. 
            This way, when many samples have an identically high GQ (as often happens for common variants), 
            the same few samples don't get selected repeatedly for all common variants.

            :param hl.MatrixTable mt: Input MatrixTable.
            :return: Tuple of genotype quality (made into a negative integer) and random float.
            :rtype: Tuple[hl.expr.Int32Expression, hl.expr.Float64Expression]
            """
        return -mt.GQ, hl.rand_unif(0, 1, seed=1)

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
        gnomad_ht = get_all_gnomad_variants()

        logger.info("Extracting variants NOT present in gnomAD...")
        ht = ht.anti_join(gnomad_ht)
        ht.write(non_gnomad_var_ht_path(*tranche_data), overwrite=args.overwrite)

    if args.get_samples:
        logger.info("Filtering to variants not present in gnomAD...")
        non_gnomad_var_ht = hl.read_table(non_gnomad_var_ht_path(*tranche_data))
        mt = mt.filter_rows(hl.is_defined(non_gnomad_var_ht[mt.row_key]))

        logger.info(
            "Getting gnomAD variants that need additional samples for readviz..."
        )
        gnomad_exomes_ht = get_additional_gnomad_variants(
            data_type="exomes", input_tsv_path=args.exomes_tsv_path
        )
        gnomad_genomes_ht = get_additional_gnomad_variants(
            data_type="genomes", input_tsv_path=args.genomes_tsv_path
        )
        mt = mt.annotate_rows(
            gnomad_exomes=hl.struct(**gnomad_exomes_ht[mt.row_key]),
            gnomad_genomes=hl.struct(**gnomad_genomes_ht[mt.row_key]),
        )
        mt = mt.filter_rows(
            hl.is_missing(mt.gnomad_exomes) & hl.is_missing(mt.gnomad_genomes),
            keep=False,
        )

        logger.info(
            f"Taking up to {args.num_samples} samples per site where samples are het, hom_var, or hemi"
        )
        mt = mt.annotate_rows(
            samples_w_het_var=hl.agg.filter(
                mt.GT.is_het(),
                hl.agg.take(
                    hl.struct(s=mt.s, GQ=mt.GQ),
                    args.num_samples,
                    ordering=_sample_ordering_expr(mt),
                ),
            ),
            samples_w_hom_var=hl.agg.filter(
                mt.GT.is_hom_var() & mt.GT.is_diploid(),
                hl.agg.take(
                    hl.struct(s=mt.s, GQ=mt.GQ),
                    args.num_samples,
                    ordering=_sample_ordering_expr(mt),
                ),
            ),
            samples_w_hemi_var=hl.agg.filter(
                hemi_expr(mt.locus, mt.meta.sex_imputation.sex_karyotype, mt.GT),
                hl.agg.take(
                    hl.struct(s=mt.s, GQ=mt.GQ),
                    args.num_samples,
                    ordering=_sample_ordering_expr(mt),
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
    parser.add_argument("--slack_channel", help="Send message to Slack channel/user")
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
