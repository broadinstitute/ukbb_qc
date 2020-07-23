import argparse
import logging
from typing import List, Optional

import hail as hl

from gnomad.resources.grch37.gnomad import liftover
from gnomad.resources.grch38.gnomad import public_release
from gnomad.utils.annotations import (
    annotate_freq,
    bi_allelic_site_inbreeding_expr,
    faf_expr,
    pop_max_expr,
)
from gnomad.utils.slack import slack_notifications
from ukbb_qc.resources.basics import get_ukbb_data, logging_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.variant_qc import var_annotations_ht_path
from ukbb_qc.slack_creds import slack_token
from ukbb_qc.utils.utils import get_hists


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("generate_frequency_data")
logger.setLevel(logging.INFO)


def generate_cohort_frequency_data(
    mt: hl.MatrixTable,
    pops: List[str],
    pop_expr: hl.expr.StringExpression,
    sex_expr: hl.expr.StringExpression,
) -> hl.Table:
    """
    Generates frequency struct annotation containing AC, AF, AN, and homozygote count for all individuals of specified ancestry from input dataset.

    :param MatrixTable mt: Input MatrixTable annotated with sample metadata information.
    :param list pops: Populations to include in frequency calculation.
    :param StringExpression pop_expr: Expression containing population information.
    :param StringExpression sex_expr: Expression containing sex karyotype information.
    :return: Table with frequency annotations in struct named `cohort_freq` and metadata in globals named `freq_meta`.
    :rtype: hl.Table
    """
    logger.info(
        f"Filtering to {pops} to calculate cohort frequency (includes related samples)..."
    )
    pops = hl.literal(pops)
    mt = mt.annotate_cols(_pop=pop_expr, _sex=sex_expr)
    mt = mt.filter_cols(pops.contains(mt._pop))
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    logger.info("Generating frequency data...")
    mt = annotate_freq(mt, sex_expr=mt._sex)
    mt = mt.drop("_pop", "_sex")
    mt = mt.select_rows(cohort_freq=mt.freq)
    cohort_freq_meta = hl.eval(mt.freq_meta)
    for i in range(len(cohort_freq_meta)):
        if i == 1:
            update_dict = {"group": "cohort_raw"}
        else:
            update_dict = {"group": "cohort"}
        cohort_freq_meta[i].update(update_dict)
    mt = mt.annotate_globals(cohort_freq_meta=cohort_freq_meta)
    return mt.select_globals("cohort_freq_meta").rows()


def generate_frequency_data(
    mt: hl.MatrixTable,
    freeze: int,
    pops_to_remove_for_popmax: List[str],
    cohort_frequency_pops: List[str],
    platform_strata: bool = False,
    tranche: bool = False,
    platform_expr: Optional[hl.expr.StringExpression] = None,
) -> hl.Table:
    """
    Generates frequency struct annotation containing AC, AF, AN, and homozygote count for input dataset stratified by population.

    Option to stratify by platform (either inferred platform or tranche as proxy for platform).

    .. note::

        This function expects dense format data.
        Frequency data is generated using gnomAD PC project population labels as the main population stratification.
        Hybrid population labels are used as the subpopulation stratification.
        However, cohort frequency (frequency including all related samples) data is generated using hybrid population labels.

    :param MatrixTable mt: Input MatrixTable.
    :param int freeze: One of the data freezes.
    :param list pops_to_remove_for_popmax: List of populations to exclude from popmax calculations. 
    :param list cohort_frequency_pops: List of populations to include in cohort frequency calculations.
    :param bool platform_strata: Whether to calculate frequencies per platform or tranche as proxy for platform.
    :param bool tranche: Whether to use tranche (as a proxy for platform) instead of inferred platform.
    :param bool calculate_by_tranche: Whether to calculate frequencies per tranche.
    :param hl.expr.StringExpression platform_expr: Expression containing platform or tranche information. Required if platform_strata is set.
    :return: Table with frequency annotations in struct named `freq` and metadata in globals named `freq_meta`.
    :rtype: Table
    """
    if platform_strata:
        if tranche:
            name = "tranche"
        else:
            name = "platform"
        additional_strata_expr = {f"{name}": platform_expr}
    else:
        additional_strata_expr = None

    logger.info("Generating cohort frequency data...")
    cohort_ht = generate_cohort_frequency_data(
        mt,
        cohort_frequency_pops,
        mt.meta.gnomad_pc_project_pop_data.pop,
        mt.meta.sex_imputation.sex_karyotype,
    )
    mt = mt.annotate_rows(cohort_freq=cohort_ht[mt.row_key].cohort_freq)
    mt = mt.annotate_globals(**cohort_ht.index_globals())

    logger.info("Filtering related samples...")
    mt = mt.filter_cols(~mt.meta.sample_filters.related)

    logger.info("Calculating InbreedingCoefficient...")
    # NOTE: This is not the ideal location to calculate this, but added here to avoid another densify
    mt = mt.annotate_rows(InbreedingCoeff=bi_allelic_site_inbreeding_expr(mt.GT))
    mt.describe()

    logger.info("Generating frequency data...")
    mt = annotate_freq(
        mt,
        sex_expr=mt.meta.sex_imputation.sex_karyotype,
        pop_expr=mt.meta.gnomad_pc_project_pop_data.pop,
        subpop_expr=mt.meta.hybrid_pop_data.pop,
        additional_strata_expr=additional_strata_expr,
    )

    # NOTE: faf and popmax are calculated on freq, not freq_expr
    # We do not want to include relateds when calculating faf and popmax
    faf, faf_meta = faf_expr(mt.freq, mt.freq_meta, mt.locus, pops_to_remove_for_popmax)
    popmax = pop_max_expr(mt.freq, mt.freq_meta, pops_to_remove_for_popmax)
    logger.info("Calculating faf and popmax...")
    mt = mt.annotate_rows(faf=faf, popmax=popmax)

    # Add cohort frequency to last index of freq struct
    # Also add cohort frequency meta to freq_meta globals
    mt = mt.select_rows(
        "InbreedingCoeff", "faf", "popmax", freq=mt.freq.extend(mt.cohort_freq),
    )
    mt = mt.select_globals(
        faf_meta=faf_meta, freq_meta=mt.freq_meta.extend(mt.cohort_freq_meta)
    )

    logger.info("Getting quality histograms...")
    mt = get_hists(mt, freeze)
    return mt.rows()


def join_gnomad(ht: hl.Table, data_type: str) -> hl.Table:
    """
    Joins UKBB frequency Table to gnomAD Tables (v2 exomes lifted Table and v3 genomes Table).

    Adds gnomAD freq, popmax, and faf structs to UKBB Table.

    :param Table ht: Input UKBB frequency Table.
    :param str data_type: One of 'exomes' or 'genomes'.
    :return: UKBB frequency Table with gnomAD frequency information added as annotation.
    :rtype: Table
    """
    if data_type == "exomes":
        gnomad_ht = (
            liftover(data_type)
            .ht()
            .select("freq", "popmax", "faf")
            .select_globals("freq_meta", "popmax_index_dict", "faf_index_dict")
        )
        ht = ht.join(gnomad_ht, how="left")
        ht = ht.rename(
            {
                "freq_1": f"gnomad_{data_type}_freq",
                "popmax_1": f"gnomad_{data_type}_popmax",
                "faf_1": f"gnomad_{data_type}_faf",
                "freq_meta_1": f"gnomad_{data_type}_freq_meta",
                "popmax_index_dict": f"gnomad_{data_type}_popmax_index_dict",
                "faf_index_dict": f"gnomad_{data_type}_faf_index_dict",
            }
        )

    else:
        gnomad_ht = (
            public_release(data_type)
            .ht()
            .select("freq", "popmax", "faf")
            .select_globals("freq_meta", "faf_index_dict")
        )
        ht = ht.join(gnomad_ht, how="left")
        ht = ht.rename(
            {
                "freq_1": f"gnomad_{data_type}_freq",
                "popmax_1": f"gnomad_{data_type}_popmax",
                "faf_1": f"gnomad_{data_type}_faf",
                "freq_meta_1": f"gnomad_{data_type}_freq_meta",
                "faf_index_dict": f"gnomad_{data_type}_faf_index_dict",
            }
        )
    return ht


def main(args):

    hl.init(log="/frequency_generation.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze

    try:
        pops_to_remove_for_popmax = set(args.pops_to_remove_for_popmax.split(","))
        cohort_frequency_pops = set(args.pops_to_include_for_cohort.split(","))
        logger.info(
            f"Excluding {pops_to_remove_for_popmax} from popmax and faf calculations"
            f"Including {cohort_frequency_pops} in cohort frequency calculations"
        )

        # Read in hardcalls
        mt = get_ukbb_data(data_source, freeze, meta_root="meta")

        if args.compute_frequency:
            logger.info("Densifying...")
            mt = hl.experimental.densify(mt)
            mt = mt.filter_rows(hl.len(mt.alleles) > 1)

            # Filter to high quality UKBB samples only
            mt = mt.filter_cols(mt.meta.sample_filters.high_quality)
            mt = mt.filter_cols(hl.is_defined(mt.meta.ukbb_meta.batch))
            mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

            if args.by_tranche:
                platform_expr = mt.meta.ukbb_meta.batch
            elif args.by_platform:
                platform_expr = mt.meta.platform_inference.qc_platform
            else:
                platform_expr = None

            logger.info("Calculating frequencies")
            ht = generate_frequency_data(
                mt,
                freeze,
                pops_to_remove_for_popmax,
                cohort_frequency_pops,
                args.by_platform,
                args.by_tranche,
                platform_expr,
            )

            ht = ht.naive_coalesce(args.n_partitions)
            ht = ht.checkpoint(
                var_annotations_ht_path("ukb_freq", data_source, freeze), args.overwrite
            )

        if args.join_gnomad:
            ht = hl.read_table(var_annotations_ht_path("ukb_freq", data_source, freeze))

            logger.info("Joining UKBB ht to gnomAD exomes and genomes liftover hts")
            ht = join_gnomad(ht, "exomes")
            ht = join_gnomad(ht, "genomes")
            ht.write(
                var_annotations_ht_path("join_freq", data_source, freeze),
                args.overwrite,
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(data_source, freeze))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--n_partitions",
        help="Desired number of partitions for output HTs",
        default=5000,
        type=int,
    )
    parser.add_argument(
        "--pops_to_remove_for_popmax",
        help="Pops to exclude for popmax",
        default="asj,fin,oth",
    )
    parser.add_argument(
        "--pops_to_include_for_cohort",
        help="Pops to include in cohort frequency calculations",
    )
    parser.add_argument(
        "--compute_frequency", help="Compute frequency data data", action="store_true",
    )
    platform_args = parser.add_mutually_exclusive_group()
    platform_args.add_argument(
        "--by_platform",
        help="Calculate frequencies stratified by inferred platform",
        action="store_true",
    )
    platform_args.add_argument(
        "--by_tranche", help="Use tranche as a proxy for platform", action="store_true",
    )
    parser.add_argument(
        "-j",
        "--join_gnomad",
        help="Join table with gnomAD tables to get gnomAD frequencies",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
