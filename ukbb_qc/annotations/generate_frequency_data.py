import argparse
import logging
from typing import List

import hail as hl

from gnomad.resources.grch37.gnomad import liftover
from gnomad.resources.grch38.gnomad import public_release
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.slack_creds import slack_token
from gnomad.utils.annotations import (
    age_hists_expr,
    annotate_freq,
    bi_allelic_site_inbreeding_expr,
    get_adj_expr,
    faf_expr,
    pop_max_expr,
    qual_hist_expr,
)
from gnomad.utils.file_utils import write_temp_gcs
from gnomad.utils.filtering import filter_to_adj
from gnomad.utils.slack import slack_notifications
from ukbb_qc.resources.basics import (
    array_sample_map_ht_path,
    capture_ht_path,
    get_ukbb_data,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.variant_qc import var_annotations_ht_path
from ukbb_qc.utils.utils import get_age_ht


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("generate_frequency_data")
logger.setLevel(logging.INFO)


def get_hists(mt: hl.MatrixTable, freeze: int) -> hl.MatrixTable:
    """
    Gets age (at recruitment; field 21022) and qual hists for UKBB

    :param MatrixTable mt: Input MatrixTable
    :param int freeze: One of the data freezes
    :return: MatrixTable with age and qual hists
    :rtype: MatrixTable
    """
    logger.info("Getting age hists")
    age_ht = get_age_ht(freeze)
    mt = mt.annotate_cols(age=age_ht[mt.col_key].age)
    mt = mt.annotate_rows(**age_hists_expr(mt.adj, mt.GT, mt.age))

    logger.info("Annotating with qual hists (on raw data)")
    mt = mt.annotate_rows(**qual_hist_expr(mt.GT, mt.GQ, mt.DP, mt.AD))

    logger.info("Annotating with qual hists (adj)")
    mt_filt = filter_to_adj(mt)
    mt_filt = mt_filt.select_rows()
    mt_filt = mt_filt.annotate_rows(
        **qual_hist_expr(mt_filt.GT, mt_filt.GQ, mt_filt.DP, mt_filt.AD)
    )
    ht = mt_filt.rows()
    mt = mt.annotate_rows(adj_qual_hists=ht[mt.row_key])
    return mt


def calculate_inbreeding_coefficient(
    mt: hl.MatrixTable, data_source: str, freeze: int, overwrite: bool
) -> None:
    """
    Calculate inbreeding coefficient annotation and writes Table.

    :param MatrixTable mt: Input MatrixTable with related samples removed.
    :param str data_source: One of "regeneron" or "broad".
    :param int freeze: One of the data freezes.
    :param bool overwrite: Whether to overwrite data.
    :return: None
    :rtype: None
    """
    ht = mt.annotate_rows(InbreedingCoeff=bi_allelic_site_inbreeding_expr(mt.GT)).rows()
    ht = ht.select("InbreedingCoeff")
    ht.write(
        var_annotations_ht_path("inbreeding_coefficient", data_source, freeze),
        overwrite=overwrite,
    )


def generate_cohort_frequency_data(
    mt: hl.MatrixTable, pops: List[str]
) -> hl.MatrixTable:
    """
    Generates frequency struct annotation containing AC, AF, AN, and homozygote count for all NFE individuals from input dataset.

    Option to stratify by tranche as proxy for platform.

    :param MatrixTable mt: Input MatrixTable annotated with sample metadata information.
    :param list pops: Populations to include in frequency calculation.
    :return: MatrixTable with frequency annotations in struct named `cohort_freq` and metadata in globals named `freq_meta`.
    :rtype: hl.MatrixTable
    """
    logger.info(
        f"Filtering to {pops} to calculate cohort frequency (includes related samples)..."
    )
    pops = hl.literal(pops)
    mt = mt.filter_cols(pops.contains(mt.meta.hybrid_pop_data.pop))
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    logger.info("Generating frequency data...")
    mt = annotate_freq(mt, sex_expr=mt.meta.sex_imputation.sex_karyotype,)
    mt = mt.transmute_rows(cohort_freq=mt.freq).select_rows("cohort_freq")
    cohort_freq_meta = hl.eval(mt.freq_meta)
    for i in range(len(cohort_freq_meta)):
        if i == 0:
            update_dict = {"group": "cohort_adj"}
        elif i == 1:
            update_dict = {"group": "cohort_raw"}
        else:
            update_dict = {"group": "cohort"}
        cohort_freq_meta[i].update(update_dict)
    mt = mt.annotate_globals(cohort_freq_meta=cohort_freq_meta)
    return mt.select_globals("cohort_freq_meta")


def generate_frequency_data(
    mt: hl.MatrixTable,
    data_source: str,
    freeze: int,
    pops_to_remove_for_popmax: List[str],
    cohort_frequency_pops: List[str],
    overwrite: bool,
    calculate_by_tranche: bool,
) -> hl.Table:
    """
    Generates frequency struct annotation containing AC, AF, AN, and homozygote count for input dataset stratified by population.

    Option to stratify by tranche as proxy for platform.

    :param MatrixTable mt: Input MatrixTable.
    :param str data_source: One of "regeneron" or "broad".
    :param int freeze: One of the data freezes.
    :param list pops_to_remove_for_popmax: List of populations to exclude from popmax calculations. 
    :param list cohort_frequency_pops: List of populations to include in cohort frequency calculations.
    :param bool overwrite: Whether to overwrite data.
    :param bool calculate_by_tranche: Whether to calculate frequencies per tranche.
    :return: Table with frequency annotations in struct named `freq` and metadata in globals named `freq_meta`.
    :rtype: Table
    """
    if calculate_by_tranche:
        additional_strata_expr = {"tranche": mt.tranche}
    else:
        additional_strata_expr = None

    logger.info("Generating cohort frequency data...")
    mt = generate_cohort_frequency_data(mt, cohort_frequency_pops)

    logger.info("Filtering related samples...")
    mt = mt.filter_cols(~mt.meta.sample_filters.related)

    logger.info("Calculating InbreedingCoefficient...")
    # NOTE: This is not the ideal location to calculate this, but added here to avoid another densify
    calculate_inbreeding_coefficient(mt, data_source, freeze, overwrite)

    logger.info("Generating frequency data...")
    mt = annotate_freq(
        mt,
        sex_expr=mt.meta.sex_imputation.sex_karyotype,
        pop_expr=mt.meta.gnomad_pc_project_data.pop,
        subpop_expr=mt.meta.hybrid_pop_data.pop,
        additional_strata_expr=additional_strata_expr,
    )

    # Add cohort frequency to last index of freq struct
    # Also add cohort frequency meta to freq_meta globals
    freq_expr = mt.freq.extend(mt.cohort_freq)
    freq_meta = mt.freq_meta.extend(mt.cohort_freq_meta)
    mt = mt.select_globals()

    logger.info("Calculating faf...")
    faf, faf_meta = faf_expr(mt.freq, mt.freq_meta, mt.locus, pops_to_remove_for_popmax)
    mt = mt.select_rows(
        freq=freq_expr,
        faf=faf,
        popmax=pop_max_expr(mt.freq, mt.freq_meta, pops_to_remove_for_popmax),
    )
    mt = mt.annotate_globals(faf_meta=faf_meta, freq_meta=freq_meta)

    logger.info("Getting hists")
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
    pops_to_remove_for_popmax = set(args.pops_to_remove_for_popmax.split(","))
    cohort_frequency_pops = set(args.pops_to_include_for_cohort.split(","))
    logger.info(
        f"Excluding {pops_to_remove_for_popmax} from popmax and faf calculations"
        f"Including {cohort_frequency_pops} in cohort frequency calculations"
    )

    # Read in hardcalls and filter to high quality UKBB samples only
    mt = get_ukbb_data(data_source, freeze, ukbb_samples_only=True, meta_root="meta")
    mt = mt.filter_cols(mt.meta.sample_filters.high_quality)

    if not args.skip_densify:
        logger.info("Reading in capture HT a")
        capture_ht = hl.read_table(capture_ht_path(data_source))

        logger.info("Computing adj and sex adjusted genotypes")
        mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)
        mt = mt.select_entries(
            GT=adjusted_sex_ploidy_expr(mt.locus, mt.GT, mt.meta.sex_karyotype),
            adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD),
        )

        logger.info("Densifying...")
        mt = hl.experimental.densify(mt)
        mt = mt.filter_rows(hl.len(mt.alleles) > 1)
        mt = mt.filter_rows(hl.is_defined(capture_ht[mt.locus]))
        mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    if args.by_platform:
        logger.info(
            "Annotating sample with tranche information (from array sample map ht)"
        )
        sample_map_ht = hl.read_table(array_sample_map_ht_path(freeze))
        mt = mt.annotate_cols(tranche=sample_map_ht[mt.col_key].batch)

        logger.info("Calculating frequencies")
        ht = generate_frequency_data(
            mt,
            data_source,
            freeze,
            pops_to_remove_for_popmax,
            cohort_frequency_pops,
            args.overwrite,
            args.by_platform,
        )

        ht = ht.repartition(args.n_partitions)
        ht = ht.checkpoint(
            var_annotations_ht_path("ukb_freq", data_source, freeze,), args.overwrite,
        )

    if args.join_gnomad:
        ht = hl.read_table(var_annotations_ht_path("ukb_freq", data_source, freeze,))

        logger.info("Joining UKBB ht to gnomAD exomes and genomes liftover hts")
        ht = join_gnomad(ht, "exomes")
        ht = join_gnomad(ht, "genomes")
        write_temp_gcs(
            ht,
            var_annotations_ht_path("join_freq", data_source, freeze,),
            args.overwrite,
        )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--n_partitions",
        help="Desired number of partitions for output HT",
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
        "--skip_densify", help="Skip densifying data", action="store_true",
    )
    parser.add_argument(
        "--calculate_frequencies",
        help="Calculate most frequency data",
        action="store_true",
    )
    parser.add_argument(
        "--by_tranche",
        help="Also calculate frequencies by UKBB data tranche",
        action="store_true",
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
