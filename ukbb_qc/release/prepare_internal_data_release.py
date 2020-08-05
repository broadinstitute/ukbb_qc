import argparse
import logging
from typing import Union

import hail as hl

from gnomad.resources.grch38.reference_data import dbsnp, lcr_intervals
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import (
    AS_FIELDS,
    RF_FIELDS,
    SITE_FIELDS,
    VQSR_FIELDS,
)
from ukbb_qc.resources.basics import (
    get_ukbb_data,
    logging_path,
    release_ht_path,
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


# Add interval QC pass to RF fields (to pull field from RF HT)
ADDITIONAL_RF_FIELDS = ["filters", "in_capture_interval", "interval_qc_pass"]
RF_FIELDS.extend(ADDITIONAL_RF_FIELDS)

# Remove InbreedingCoeff from allele-specific fields (processed separately from other fields)
AS_FIELDS.remove("InbreedingCoeff")


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

    return hl.struct(
        lcr=hl.is_defined(lcr_ht[t.locus]),
        fail_interval_qc=~(t.rf.interval_qc_pass),
        in_capture_region=t.rf.in_capture_interval,
    )


def prepare_annotations(
    ht: hl.Table,
    info_ht: hl.Table,
    freq_ht: hl.Table,
    rf_ht: hl.Table,
    vep_ht: hl.Table,
    allele_ht: hl.Table,
    vqsr_ht: hl.Table,
) -> hl.MatrixTable:
    """
    Join all Tables with variant annotations.

    :param Table ht: Table to be annotated.
    :param Table freq_ht: Table with frequency annotations.
    :param Table rf_ht: Table with random forest variant annotations.
    :param Table vep_ht: Table with VEP variant annotations.
    :param Table allele_ht: Table containing allele annotations.
    :param Table vqsr_ht: Table containing VQSR annotations.
    :return: Table containing joined annotations.
    :rtype: hl.Table
    """
    logger.info("Selecting annotations from annotation Tables...")
    info_fields = SITE_FIELDS + AS_FIELDS
    info_ht = info_ht.transmute(info=info_ht.info.select(*info_fields)).select(
        "info", "qual"
    )
    info_ht = info_ht.transmute(
        info=info_ht.info.annotate(
            sibling_singleton=rf_ht[info_ht.key].sibling_singleton,
            transmitted_singleton=rf_ht[info_ht.key].transmitted_singleton,
            InbreedingCoeff=freq_ht[info_ht.key].InbreedingCoeff,
        )
    )
    rf_ht = rf_ht.select(*RF_FIELDS)
    vep_ht = vep_ht.transmute(vep=vep_ht.vep.drop("colocated_variants"))
    vqsr_ht = vqsr_ht.transmute(info=vqsr_ht.info.select(*VQSR_FIELDS)).select("info")
    dbsnp_ht = dbsnp.ht().select("rsid")

    logger.info(
        "Annotating HT with random forest HT and flagging problematic regions..."
    )
    ht = ht.annotate(rf=rf_ht[ht.key])
    ht = ht.annotate_globals(rf_globals=rf_ht.index_globals())
    ht = ht.annotate(region_flag=flag_problematic_regions(ht))
    ht = ht.transmute(rf=ht.rf.drop("interval_qc_pass", "in_capture_interval"))

    logger.info("Annotating HT with frequency information...")
    ht = ht.annotate(**freq_ht[ht.key])
    ht = ht.annotate_globals(**freq_ht.index_globals())

    logger.info("Annotating HT info, vep, allele_info, vqsr, rsid, and qual...")
    ht = ht.annotate(
        info=info_ht[ht.key].info,
        vep=vep_ht[ht.key].vep,
        allele_info=allele_ht[ht.key].allele_data,
        vqsr=vqsr_ht[ht.key].info,
        rsid=dbsnp_ht[ht.key].rsid,
        qual=info_ht[ht.key].qual,
    )
    return ht


def main(args):

    hl.init(log="/release.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)

    try:
        logger.info("Reading in hardcalls MT and subsetting to rows...")
        # NOTE: reading in raw MatrixTable to be able to return all samples/variants
        ht = get_ukbb_data(*tranche_data).rows()

        logger.info("Filtering out low QUAL variants...")
        info_ht = hl.read_table(info_ht_path(*tranche_data))
        info_ht = info_ht.filter(~info_ht.AS_lowqual)
        ht = ht.filter(hl.is_defined(info_ht[ht.key]))

        logger.info("Removing chrM...")
        ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chrM")], keep=False)

        logger.info("Reading in all variant annotation tables...")
        freq_ht = hl.read_table(var_annotations_ht_path("join_freq", *tranche_data))
        rf_ht = hl.read_table(var_annotations_ht_path("rf", *tranche_data))
        vep_ht = hl.read_table(var_annotations_ht_path("vep", *tranche_data))
        allele_ht = hl.read_table(var_annotations_ht_path("allele_data", *tranche_data))
        vqsr_ht = hl.read_table(var_annotations_ht_path("vqsr", *tranche_data))

        logger.info("Adding annotations...")
        ht = prepare_annotations(
            ht, info_ht, freq_ht, rf_ht, vep_ht, allele_ht, vqsr_ht,
        )

        logger.info("Getting age hist data...")
        age_hist_data = get_age_distributions(get_age_ht(freeze))
        ht = ht.annotate_globals(age_distribution=age_hist_data)
        ht = ht.naive_coalesce(args.n_partitions)
        ht.write(release_ht_path(*tranche_data), args.overwrite)
        logger.info(f"Final HT count: {ht.count()}")

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(*tranche_data))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--n_partitions",
        help="Number of desired partitions for output Tables",
        default=10000,
        type=int,
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
