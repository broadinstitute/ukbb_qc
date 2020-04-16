import argparse
import logging
from typing import List

import hail as hl

from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.resources.grch37.gnomad import liftover
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
from gnomad.utils.slack import try_slack
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


def generate_frequency_data(
    mt: hl.MatrixTable,
    freeze: int,
    POPS_TO_REMOVE_FOR_POPMAX: List[str],
    calculate_by_platform: bool = False,
) -> hl.Table:
    """
    Generates frequency struct annotation containing AC, AF, AN, and homozygote count for dataset stratified by population.

    Option to stratify by tranche as proxy for platform.

    :param MatrixTable mt: Input MatrixTable.
    :param int freeze: One of the data freezes
    :param list POPS_TO_REMOVE_FOR_POPMAX: List of populations to exclude from popmax calculations.
    :param bool calculate_by_platform: Whether to calculate frequencies per tranche.
    :return: Table with frequency annotations in struct named `freq` and metadata in globals named `freq_meta`
    :rtype: Table
    """
    if calculate_by_platform:
        additional_strata_expr = {"tranche": mt.tranche}
    else:
        additional_strata_expr = None

    logger.info("Generating frequency data...")
    mt = annotate_freq(
        mt,
        sex_expr=mt.meta.sex_imputation.sex_karyotype,
        pop_expr=mt.meta.hybrid_pop_data.hybrid_pop,
        additional_strata_expr=additional_strata_expr,
    )

    # Select freq, FAF and popmax
    faf, faf_meta = faf_expr(mt.freq, mt.freq_meta, mt.locus, POPS_TO_REMOVE_FOR_POPMAX)
    mt = mt.select_rows(
        "freq",
        faf=faf,
        popmax=pop_max_expr(mt.freq, mt.freq_meta, POPS_TO_REMOVE_FOR_POPMAX),
    )
    mt = mt.annotate_globals(faf_meta=faf_meta)

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
            hl.read_table(
                "gs://gnomad-public/release/3.0/ht/genomes/gnomad.genomes.r3.0.sites.ht"
            )
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
    POPS_TO_REMOVE_FOR_POPMAX = set(args.pops.split(","))
    logger.info(
        f"Excluding {POPS_TO_REMOVE_FOR_POPMAX} from popmax and faf calculations"
    )

    # Read in hardcalls and filter to high quality UKBB samples only
    mt = get_ukbb_data(data_source, freeze, ukbb_samples_only=True, meta_root="meta")
    mt = mt.filter_cols(mt.meta.sample_filters.high_quality)
    mt = mt.filter_rows(hl.agg.any(mt.LGT.is_non_ref() | hl.is_defined(mt.END)))

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

    if args.filter_related:
        logger.info("Filtering out related samples and their variants...")
        mt = mt.filter_cols(~mt.meta.sample_filters.related)
        mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    if args.by_platform:
        logger.info(
            "Annotating sample with tranche information (from array sample map ht)"
        )
        sample_map_ht = hl.read_table(array_sample_map_ht_path(freeze))
        mt = mt.annotate_cols(tranche=sample_map_ht[mt.col_key].batch)

    if args.calculate_frequencies:
        logger.info("Calculating InbreedingCoeff to avoid another densify")
        # Note: this is not the ideal location for this, but adding here to avoid another densify
        # TODO: Still need to decide how to handle variants where there are no hets, this returns NA
        if args.filter_related:
            mt = mt.annotate_rows(
                InbreedingCoeff=bi_allelic_site_inbreeding_expr(mt.GT)
            )

        logger.info("Calculating frequencies")
        ht = generate_frequency_data(
            mt, freeze, POPS_TO_REMOVE_FOR_POPMAX, args.by_platform
        )

        ht = ht.repartition(args.n_partitions)
        ht = ht.checkpoint(
            var_annotations_ht_path(
                f'ukb_freq_hybrid{"_unrelated" if args.filter_related else ""}',
                data_source,
                freeze,
            ),
            args.overwrite,
        )

        if args.filter_related:
            ht.select_rows("InbreedingCoeff").write(
                var_annotations_ht_path("inbreeding_coeff", data_source, freeze),
                args.overwrite,
            )

    if args.join_gnomad:
        ht = hl.read_table(
            var_annotations_ht_path(
                f'ukb_freq_hybrid{"_unrelated" if args.filter_related else ""}',
                data_source,
                freeze,
            )
        )

        logger.info("Joining UKBB ht to gnomAD exomes and genomes liftover hts")
        ht = join_gnomad(ht, "exomes")
        ht = join_gnomad(ht, "genomes")
        write_temp_gcs(
            ht,
            var_annotations_ht_path(
                f'join_freq_hybrid{"_unrelated" if args.filter_related else ""}',
                data_source,
                freeze,
            ),
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
        "--pops", help="Pops to exclude for popmax", default="asj,fin,oth"
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
        "--by_platform",
        help="Also calculate frequencies by platform",
        action="store_true",
    )
    parser.add_argument(
        "--filter_related",
        help="Filter related individuals before calculating allele frequencies",
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
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
