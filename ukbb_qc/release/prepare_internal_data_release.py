import argparse
import logging
from typing import Union

import hail as hl

from gnomad.resources.grch38.reference_data import dbsnp, lcr_intervals
from gnomad.utils.slack import slack_notifications
from ukbb_qc.resources.basics import (
    capture_ht_path,
    get_ukbb_data,
    logging_path,
    release_ht_path,
    release_mt_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.variant_qc import info_ht_path, var_annotations_ht_path
from ukbb_qc.slack_creds import slack_token
from ukbb_qc.utils.utils import get_age_distributions, get_age_ht


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("data_release")
logger.setLevel(logging.INFO)

AS_FIELDS = [
    "AS_BaseQRankSum",
    "AS_FS",
    "AS_MQ",
    "AS_MQRankSum",
    "AS_pab_max",
    "AS_QD",
    "AS_ReadPosRankSum",
    "AS_SOR",
    "AS_VarDP",
]
"""
Allele-specific variant annotations to be extracted from the info HT.
"""
SITE_FIELDS = [
    "FS",
    "InbreedingCoeff",
    "MQ",
    "MQRankSum",
    "QD",
    "ReadPosRankSum",
    "sibling_singleton",
    "SOR",
    "transmitted_singleton",
    "VarDP",
]
"""
Site level variant annotations to be extracted from the info HT.
"""
RF_FIELDS = [
    "interval_qc_pass",
    "fail_hard_filters",
    "filters",
    "rf_label",
    "rf_train",
    "rf_probability",
    "tp",
]
"""
Annotations specific to the random forest model to be extracted from the random forest HT.
"""
VQSR_FIELDS = ["AS_VQSLOD", "AS_culprit", "NEGATIVE_TRAIN_SITE", "POSITIVE_TRAIN_SITE"]
"""
Annotations specific to VQSR to be extracted from the VQSR HT.
"""


def flag_problematic_regions(
    t: Union[hl.Table, hl.MatrixTable]
) -> hl.expr.StructExpression:
    """
    Creates `region_flag` struct. 

    Struct contains flags for problematic regions and whether variant is in capture interval.
    
    .. note:: 

        No hg38 resources for decoy, segdup, or self chain yet.

    :param Table/MatrixTable t: Input Table/MatrixTable.
    :return: `region_flag` struct row annotation.
    :rtype: hl.expr.StructExpression
    """
    lcr_ht = lcr_intervals.ht()
    capture_ht = hl.read_table(capture_ht_path("broad"))

    return hl.struct(
        lcr=hl.is_defined(lcr_ht[t.locus]),
        fail_interval_qc=~(t.rf.interval_qc_pass),
        in_capture_region=hl.is_defined(capture_ht[t.locus]),
    )


def prepare_annotations(
    mt: hl.MatrixTable,
    info_ht: hl.Table,
    freq_ht: hl.Table,
    rf_ht: hl.Table,
    vep_ht: hl.Table,
    allele_ht: hl.Table,
    vqsr_ht: hl.Table,
) -> hl.MatrixTable:
    """
    Join all Tables with variant annotations.

    Also removes variants on chrM.

    :param MatrixTable mt: MatrixTable to be annotated.
    :param Table freq_ht: Table with frequency annotations.
    :param Table rf_ht: Table with random forest variant annotations.
    :param Table vep_ht: Table with VEP variant annotations.
    :param Table allele_ht: Table containing allele annotations.
    :param Table vqsr_ht: Table containing VQSR annotations.
    :return: Table containing joined annotations.
    :rtype: hl.Table
    """
    logger.info("Filtering out low QUAL variants...")
    mt = mt.filter_rows(
        (hl.agg.any(mt.GT.is_non_ref() | hl.is_defined(mt.END)))
        & (~info_ht[mt.row_key].AS_lowqual)
    )
    logger.info(f"Count after filtering out low QUAL variants {mt.count()}")

    logger.info("Selecting annotations from annotation Tables...")
    info_fields = SITE_FIELDS + AS_FIELDS
    info_ht = info_ht.transmute(info=info_ht.info.select(*info_fields)).select(
        "info", "qual"
    )
    info_ht = info_ht.transmute(
        info=info_ht.info.annotate(
            sibling_singleton=rf_ht[info_ht.key].sibling_singleton,
            transmitted_singleton=rf_ht[info_ht.key].transmitted_singleton,
        )
    )
    rf_ht = rf_ht.select(*RF_FIELDS)
    vep_ht = vep_ht.transmute(vep=vep_ht.vep.drop("colocated_variants"))
    vqsr_ht = vqsr_ht.transmute(info=vqsr_ht.info.select(*VQSR_FIELDS)).select("info")
    # NOTE: will need to nest qual hist fields under qual_hists struct for 300K
    # gq_hist_alt, gq_hist_all, dp_hist_alt, dp_hist_all, ab_hist_alt
    freq_ht = freq_ht.drop(
        "InbreedingCoeff"
    )  # NOTE: Will need to drop InbreedingCoeff annotation here in 500K, but this isn't present in 300K
    dbsnp_ht = dbsnp.ht().select("rsid")

    logger.info(f"Frequency HT before filtering out chrM: {freq_ht.count()}")
    freq_ht = hl.filter_intervals(
        freq_ht, [hl.parse_locus_interval("chrM")], keep=False
    )
    logger.info(f"Frequency HT count after filtering out chrM: {freq_ht.count()}")

    logger.info(
        "Annotating MT with random forest HT and flagging problematic regions..."
    )
    mt = mt.annotate_rows(rf=rf_ht[mt.row_key])
    mt = mt.annotate_rows(region_flag=flag_problematic_regions(mt))
    mt = mt.transmute_rows(rf=mt.rf.drop("interval_qc_pass"))

    logger.info("Annotating MT with frequency information...")
    mt = mt.annotate_rows(**freq_ht[mt.row_key])
    mt = mt.annotate_globals(**freq_ht.index_globals())

    logger.info("Annotating MT info, vep, allele_info, vqsr, rsid, and qual...")
    mt = mt.annotate_rows(
        info=info_ht[mt.row_key].info,
        vep=vep_ht[mt.row_key].vep,
        allele_info=allele_ht[mt.row_key].allele_data,
        vqsr=vqsr_ht[mt.row_key].info,
        rsid=dbsnp_ht[mt.row_key].rsid,
        qual=info_ht[mt.row_key].qual,
    )
    mt = mt.annotate_globals(rf=rf_ht.index_globals())
    logger.info(f"Final MT count: {mt.count()}")
    return mt


def main(args):

    hl.init(log="/release.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze

    try:
        if args.prepare_internal_mt:
            logger.info("Getting raw MT and dropping all unnecessary entries...")
            # NOTE: reading in raw MatrixTable to be able to return all samples/variants
            mt = get_ukbb_data(
                data_source,
                freeze,
                key_by_locus_and_alleles=args.key_by_locus_and_alleles,
                split=False,
                raw=True,
                repartition=args.repartition,
                n_partitions=args.raw_partitions,
            ).select_entries(
                "DP",
                "GQ",
                "LA",
                "LAD",
                "LGT",
                "LPGT",
                "LPL",
                "MIN_DP",
                "PID",
                "RGQ",
                "SB",
            )

            logger.info("Splitting raw MT...")
            mt = hl.experimental.sparse_split_multi(mt)
            mt = mt.select_entries(
                "GT", "GQ", "DP", "AD", "MIN_DP", "PGT", "PID", "PL", "SB"
            )

            logger.info("Reading in all variant annotation tables...")
            info_ht = hl.read_table(info_ht_path(data_source, freeze))
            freq_ht = hl.read_table(
                "join_freq", var_annotations_ht_path(data_source, freeze)
            )
            rf_ht = hl.read_table(var_annotations_ht_path("rf", data_source, freeze))
            vep_ht = hl.read_table(var_annotations_ht_path("vep", data_source, freeze))
            allele_ht = hl.read_table(
                var_annotations_ht_path("allele_data", data_source, freeze)
            )
            vqsr_ht = hl.read_table(
                var_annotations_ht_path("vqsr", data_source, freeze)
            )

            logger.info("Adding annotations...")
            mt = prepare_annotations(
                mt, info_ht, freq_ht, rf_ht, vep_ht, allele_ht, vqsr_ht,
            )

            logger.info("Getting age hist data...")
            age_hist_data = get_age_distributions(get_age_ht(freeze))
            mt = mt.annotate_globals(age_distribution=age_hist_data)
            mt.write(release_mt_path(data_source, freeze), args.overwrite)

        if args.prepare_browser_ht:
            mt = hl.read_matrix_table(release_mt_path(data_source, freeze))
            logger.info(f"Release MT count: {mt.count()}")

            ht = mt.rows()
            logger.info(f"Release HT count: {ht.count()}")

            ht = ht.naive_coalesce(args.n_partitions)
            ht.write(release_ht_path(data_source, freeze), args.overwrite)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(data_source, freeze))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--key_by_locus_and_alleles",
        help="Re-key raw MatrixTable by locus and alleles. REQUIRED only for tranche 3/300K dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--repartition",
        help="Repartition raw MatrixTable on read. REQUIRED only for tranche 3/300K dataset",
        action="store_true",
    )
    parser.add_argument(
        "--raw_partitions",
        help="Number of desired partitions for raw MatrixTable. Only used for tranche 3/300K dataset",
        default=30000,
        type=int,
    )
    parser.add_argument(
        "--n_partitions",
        help="Number of desired partitions for output Tables",
        default=10000,
        type=int,
    )
    parser.add_argument(
        "--prepare_internal_mt",
        help="Prepare internal MatrixTable",
        action="store_true",
    )
    parser.add_argument(
        "--prepare_browser_ht", help="Prepare sites HT for browser", action="store_true"
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
