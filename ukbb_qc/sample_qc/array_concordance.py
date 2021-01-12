import argparse
import logging
from typing import Tuple

import hail as hl

from gnomad.utils.file_utils import file_exists
from gnomad.utils.liftover import (
    get_liftover_genome,
    default_lift_data,
    annotate_snp_mismatch,
)
from gnomad.utils.slack import slack_notifications
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.basics import (
    array_mt_path,
    array_sample_map_ht_path,
    get_array_data_path,
    get_checkpoint_path,
    get_ukbb_data,
    logging_path,
)
from ukbb_qc.resources.sample_qc import (
    array_concordance_results_path,
    array_concordance_sites_path,
)
from ukbb_qc.slack_creds import slack_token
from ukbb_qc.utils.utils import (
    remove_hard_filter_samples,
    get_sites,
    annotate_interval_qc_filter,
)


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("array_concordance")
logger.setLevel(logging.INFO)


def prepare_array_and_exome_mt(
    freeze: int,
    array_mt: hl.MatrixTable,
    exome_mt: hl.MatrixTable,
    call_rate_cutoff: float,
    af_cutoff: float,
) -> Tuple[hl.MatrixTable, hl.MatrixTable]:
    """
    Prepares array and exome MatrixTables for calculating concordance using `get_array_exome_concordance`.

    Maps array to exome IDs and filters both matrix tables to overlapping samples and high callrate, common sites.

    NOTE: Assumes input exome MatrixTable is split (hl.experimental.sparse_split_multi)
    
    :param int freeze: One of the data freezes
    :param MatrixTable array_mt: Array MatrixTable
    :param MatrixTable exome_mt: Exome MatrixTable
    :param float call_rate_cutoff: Minimum call rate cutoff for sites
    :param float af_cutoff: Minimum alternate allele frequency cutoff for sites
    :return: Prepared array and exome MatrixTables
    :rtype: hl.MatrixTable, hl.MatrixTable
    """
    logger.info("Mapping array sample names to exome sample names...")
    sample_map = hl.read_table(array_sample_map_ht_path(freeze))
    exome_mt = exome_mt.annotate_cols(
        ukbb_app_26041_id=sample_map[exome_mt.col_key].ukbb_app_26041_id
    )
    exome_mt = exome_mt.key_cols_by("ukbb_app_26041_id")
    exome_mt = exome_mt.semi_join_cols(array_mt.cols())
    array_mt = array_mt.semi_join_cols(exome_mt.cols())

    exome_array_count = array_mt.count_cols()
    logger.info(
        f"Total number of IDs in the exome data that are also in the array data: {exome_array_count}..."
    )

    logger.info("Filtering sparse exome MT to non ref sites...")
    exome_mt = exome_mt.filter_rows(hl.len(exome_mt.alleles) > 1)

    logger.info(
        "Getting high callrate, common sites in autosomes only from tranche 2/freeze 5..."
    )
    if not file_exists(array_concordance_sites_path()):
        sites_ht = get_sites(af_cutoff, call_rate_cutoff)
    else:
        sites_ht = hl.read_table(array_concordance_sites_path())

    logger.info("Filtering exome and array MT to tranche 2 sites...")
    exome_mt = exome_mt.filter_rows(hl.is_defined(sites_ht[exome_mt.row_key]))
    array_mt = array_mt.filter_rows(hl.is_defined(sites_ht[array_mt.row_key]))
    return array_mt, exome_mt


def get_array_exome_concordance(
    array_mt: hl.MatrixTable, exome_mt: hl.MatrixTable
) -> Tuple[hl.Table, hl.Table]:
    """
    Runs concordance on input array and exome MatrixTables.

    Returns both sample and variant concordance Tables.

    :param array_mt: Array MatrixTable
    :param exome_mt: Exome MatrixTable
    :return: Sample and variant concordance Tables
    """
    summary, samples, variants = hl.concordance(array_mt, exome_mt)
    variants = variants.annotate(
        num_gt_con_non_ref=(variants.concordance[3][3] + variants.concordance[4][4]),
        # NOTE: no homref calls in sparse (only check index 3 and higher in second matrix)
        # Confirmed that concordance calculated using only sparse MT has high correlation with concordance
        # calculated on dense MT using tranche 2/freeze 5/200K
        num_gt_non_ref=(
            hl.sum(variants.concordance[2][3:])
            + hl.sum(variants.concordance[3][3:])
            + hl.sum(variants.concordance[4][3:])
        ),
    )
    variants = variants.annotate(
        prop_gt_con_non_ref=variants.num_gt_con_non_ref / variants.num_gt_non_ref
    )
    samples = samples.annotate(
        num_gt_con_non_ref=(samples.concordance[3][3] + samples.concordance[4][4]),
        num_gt_non_ref=(
            hl.sum(samples.concordance[2][3:])
            + hl.sum(samples.concordance[3][3:])
            + hl.sum(samples.concordance[4][3:])
        ),
    )
    samples = samples.annotate(
        prop_gt_con_non_ref=samples.num_gt_con_non_ref / samples.num_gt_non_ref
    )

    return samples, variants


def main(args):

    hl.init(log="/array_concordance.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze
    call_rate_cutoff = args.call_rate_cutoff
    af_cutoff = args.af_cutoff
    overwrite = args.overwrite

    try:
        if args.import_arrays:
            logger.info(
                "Importing array data (autosomes only) into hail Matrix Table..."
            )

            array_files = [
                (
                    get_array_data_path("bed", f"chr{chrom}"),
                    get_array_data_path("bim", f"chr{chrom}"),
                )
                for chrom in range(1, 23)
            ]

            mts = [
                hl.import_plink(
                    bed=f[0],
                    bim=f[1],
                    fam=get_array_data_path("fam", "chr22"),
                    reference_genome="GRCh37",
                    a2_reference=False,
                    skip_invalid_loci=True,
                )
                for f in array_files
            ]

            array_mt = mts[0].union_rows(*mts[1:])
            array_mt = array_mt.checkpoint(array_mt_path(), overwrite=overwrite)
            array_variants, array_samples = array_mt.count()
            logger.info(
                f"{array_variants} variants and {array_samples} samples found in array data"
            )

        if args.liftover_arrays:
            logger.info("Lifting over array data from GRCh37 to GRCh38...")
            array_mt = hl.read_matrix_table(array_mt_path())

            logger.info("Preparing reference genomes for liftover")
            source, target = get_liftover_genome(array_mt)

            logger.info(f"Lifting data to {target.name}")
            array_mt = default_lift_data(array_mt)
            array_mt = array_mt.checkpoint(
                array_mt_path(liftover=True, checkpoint=True), overwrite=overwrite,
            )

            logger.info("Checking SNPs for reference mismatches")
            array_mt = annotate_snp_mismatch(array_mt, data_type=None, rg=target)
            array_mt = array_mt.checkpoint(
                array_mt_path(liftover=True), overwrite=overwrite
            )
            array_variants, array_samples = array_mt.count()
            logger.info(
                f"{array_variants} variants and {array_samples} samples found in liftover array data"
            )

        if args.array_concordance:
            logger.info("Checking concordance between exome and array data...")
            array_mt = hl.read_matrix_table(array_mt_path(liftover=True))
            exome_mt = get_ukbb_data(data_source, freeze, adj=True, split=True)

            logger.info("Removing hard filtered samples...")
            exome_mt = remove_hard_filter_samples(
                data_source, freeze, exome_mt, gt_field="GT"
            )
            logger.info(
                f"Count after removing hard filtered samples: {exome_mt.count()}"
            )

            if args.interval_qc_filter:
                logger.info(
                    f"Filtering exome MT to intervals with at least {args.pct_samples * 100}% of samples"
                    f"for coverage field {args.cov_filter_field}"
                )
                exome_mt = annotate_interval_qc_filter(
                    data_source,
                    freeze,
                    exome_mt,
                    cov_filter_field=args.cov_filter_field,
                    pct_samples=args.pct_samples,
                    autosomes_only=True,
                )
                exome_mt = exome_mt.filter_rows(exome_mt.interval_qc_pass)

            array_mt, exome_mt = prepare_array_and_exome_mt(
                freeze, array_mt, exome_mt, call_rate_cutoff, af_cutoff
            )

            # NOTE: for freeze 6 (300K), had to remove two samples with duplicate IDs
            # removed UKB_4048554_0301608642 and column index 204526 (higher column index for duplicate sample UKB_1223807_0330880742)
            exome_mt = exome_mt.checkpoint(
                get_checkpoint_path(
                    data_source,
                    freeze,
                    name=f"exome_subset_concordance_callrate_{call_rate_cutoff}_af_{af_cutoff}",
                    mt=True,
                ),
                overwrite=overwrite,
            )
            array_mt = array_mt.checkpoint(
                get_checkpoint_path(
                    data_source,
                    freeze,
                    name=f"array_subset_concordance_callrate_{call_rate_cutoff}_af_{af_cutoff}",
                    mt=True,
                ),
                overwrite=overwrite,
            )

            samples, variants = get_array_exome_concordance(array_mt, exome_mt)

            variants = variants.annotate_globals(
                callrate_cutoff=call_rate_cutoff, af_cutoff=af_cutoff
            )
            samples = samples.annotate_globals(
                callrate_cutoff=call_rate_cutoff, af_cutoff=af_cutoff
            )
            variants = variants.naive_coalesce(args.n_partitions)
            variants.write(
                array_concordance_results_path(data_source, freeze, sample=False),
                overwrite=overwrite,
            )
            samples = samples.repartition(args.n_partitions)
            samples.write(
                array_concordance_results_path(data_source, freeze), overwrite=overwrite
            )

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(logging_path(data_source, freeze))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--n_partitions",
        help="Number of partitions for desired output Tables. Used when writing variant and sample concordance HTs",
        default=5000,
        type=int,
    )
    array_import = parser.add_argument_group(
        "Import array plink files and write Matrix Table"
    )
    array_import.add_argument(
        "-i", "--import_arrays", help="Import array data", action="store_true"
    )
    array_import.add_argument(
        "-l",
        "--liftover_arrays",
        help="Liftover array Matrix Table to GRCh38",
        action="store_true",
    )

    concordance = parser.add_argument_group("Compute array concordance with exomes")
    concordance.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    concordance.add_argument(
        "-c",
        "--array_concordance",
        help="Compute array concordance",
        action="store_true",
    )
    concordance.add_argument(
        "--call_rate_cutoff", help="Call rate cutoff.", type=float, default=0.95
    )
    concordance.add_argument(
        "--af_cutoff",
        help="Allele frequency cutoff (also used for variant QC HT creation).",
        type=float,
        default=0.0001,
    )
    concordance.add_argument(
        "--interval_qc_filter",
        help="Should interval QC be applied",
        action="store_true",
    )
    concordance.add_argument(
        "--cov_filter_field",
        help="Field used to filter intervals",
        default="pct_samples_20x",
    )
    concordance.add_argument(
        "--pct_samples",
        help="Percent samples at specified coverage field to filter intervals",
        default=0.85,
        type=float,
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
