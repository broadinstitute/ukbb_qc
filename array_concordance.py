import hail as hl
import ukbb_qc.resources as urc
from gnomad_hail import try_slack
from gnomad_hail.utils.liftover import get_liftover_genome, lift_data, annotate_snp_mismatch
from ukbb_qc.utils import interval_qc_filter, remove_hard_filter_samples
import argparse
import logging


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("load_array_data")
logger.setLevel(logging.INFO)


def prepare_array_and_exome_mt(data_source, freeze, array_mt, exome_mt, call_rate_cutoff, af_cutoff):
    logger.info("Mapping array sample names to exome sample names...")
    sample_map = hl.read_table(urc.array_sample_map_ht("broad", freeze))
    exome_mt = exome_mt.annotate_cols(**sample_map[exome_mt.col_key])
    exome_mt = exome_mt.filter_cols(hl.is_defined(exome_mt.ukbb_app_26041_id))

    sample_map = sample_map.key_by(exome_s=sample_map.ukbb_app_26041_id)
    array_mt = array_mt.key_cols_by(s=sample_map[array_mt.s].s)
    array_mt = array_mt.filter_cols(hl.is_defined(array_mt.s))

    exome_array_count = array_mt.count_cols()
    logger.info(
        f"Total number of IDs in the sample map that are also in the array data: {exome_array_count}..."
    )

    array_mt = hl.variant_qc(array_mt)
    array_mt = array_mt.checkpoint(
        urc.get_mt_checkpoint_path(data_source, freeze, name=f"array_mt_variant_qc"),
        overwrite=True,
    )
    array_mt = array_mt.filter_rows(
        (array_mt.variant_qc.call_rate > call_rate_cutoff)
        & (array_mt.variant_qc.AF[1] > af_cutoff)
    )

    # NOTE: Filter to autosomes because of adjusted sex ploidies in hardcalls mt (hail throws a ploidy 0 error)
    exome_mt = exome_mt.filter_rows(exome_mt.locus.in_autosome())
    exome_mt = exome_mt.filter_rows(
        hl.is_defined(array_mt.index_rows(exome_mt.row_key))
    )

    exome_mt = hl.variant_qc(exome_mt)
    exome_mt = exome_mt.checkpoint(
        urc.get_mt_checkpoint_path(data_source, freeze, name=f"exome_mt_variant_qc"),
        overwrite=True,
    )
    exome_mt = exome_mt.filter_rows(
        (exome_mt.variant_qc.call_rate > call_rate_cutoff)
        & (exome_mt.variant_qc.AF[1] > af_cutoff)
    )

    array_mt = array_mt.filter_rows(
        hl.is_defined(exome_mt.index_rows(array_mt.row_key))
    )

    return array_mt, exome_mt


def get_array_exome_concordance(array_mt, exome_mt):
    summary, samples, variants = hl.concordance(array_mt, exome_mt)
    variants = variants.annotate(
        num_gt_con=(
            variants.concordance[2][2]
            + variants.concordance[3][3]
            + variants.concordance[4][4]
        ),
        num_gt=(
            hl.sum(variants.concordance[2][2:])
            + hl.sum(variants.concordance[3][2:])
            + hl.sum(variants.concordance[4][2:])
        ),
        num_gt_con_non_ref=(variants.concordance[3][3] + variants.concordance[4][4]),
        num_gt_non_ref=(
            hl.sum(variants.concordance[2][3:])
            + hl.sum(variants.concordance[3][2:])
            + hl.sum(variants.concordance[4][2:])
        ),
    )
    variants = variants.annotate(
        prop_gt_con=variants.num_gt_con / variants.num_gt,
        prop_gt_con_non_ref=variants.num_gt_con_non_ref / variants.num_gt_non_ref,
    )
    samples = samples.annotate(
        num_gt_con=(
            samples.concordance[2][2]
            + samples.concordance[3][3]
            + samples.concordance[4][4]
        ),
        num_gt=(
            hl.sum(samples.concordance[2][2:])
            + hl.sum(samples.concordance[3][2:])
            + hl.sum(samples.concordance[4][2:])
        ),
        num_gt_con_non_ref=(samples.concordance[3][3] + samples.concordance[4][4]),
        num_gt_non_ref=(
            hl.sum(samples.concordance[2][3:])
            + hl.sum(samples.concordance[3][2:])
            + hl.sum(samples.concordance[4][2:])
        ),
    )
    samples = samples.annotate(
        prop_gt_con=samples.num_gt_con / samples.num_gt,
        prop_gt_con_non_ref=samples.num_gt_con_non_ref / samples.num_gt_non_ref,
    )

    return samples, variants


def main(args):
    hl.init(log="/array_concordance.log")
    data_source = args.data_source
    freeze = args.freeze
    call_rate_cutoff = args.call_rate_cutoff
    af_cutoff = args.af_cutoff
    overwrite = args.overwrite

    if args.import_arrays:
        logger.info("Importing array data (autosomes only) into hail Matrix Table...")

        array_files = [
            (
                urc.get_array_data_path("bed", f"chr{chrom}"),
                urc.get_array_data_path("bim", f"chr{chrom}"),
            )
            for chrom in range(1, 23)
        ]

        mts = [
            hl.import_plink(
                bed=f[0],
                bim=f[1],
                fam=urc.get_array_data_path("fam", "chr22"),
                reference_genome="GRCh37",
                a2_reference=False,
                skip_invalid_loci=True,
            )
            for f in array_files
        ]

        array_mt = mts[0].union_rows(*mts[1:])
        array_mt = array_mt.checkpoint(urc.array_mt_path(), overwrite=overwrite)
        array_variants, array_samples = array_mt.count()
        logger.info(
            f"{array_variants} variants and {array_samples} samples found in array data"
        )

    if args.liftover_arrays:
        logger.info("Lifting over array data from GRCh37 to GRCh38...")
        array_mt = hl.read_matrix_table(urc.array_mt_path())

        logger.info("Preparing reference genomes for liftover")
        source, target = get_liftover_genome(array_mt)

        logger.info(f"Lifting data to {target.name}")
        array_mt = lift_data(
            array_mt,
            gnomad=False,
            data_type=None,
            path=urc.array_mt_path(liftover=True, checkpoint=True),
            rg=target,
            overwrite=overwrite
        )

        logger.info("Checking SNPs for reference mismatches")
        array_mt = annotate_snp_mismatch(array_mt, data_type=None, rg=target)
        array_mt = array_mt.checkpoint(
            urc.array_mt_path(liftover=True), overwrite=overwrite
        )
        array_variants, array_samples = array_mt.count()
        logger.info(
            f"{array_variants} variants and {array_samples} samples found in liftover array data"
        )

    if args.array_concordance:
        logger.info("Checking concordance between exome and array data...")
        array_mt = hl.read_matrix_table(urc.array_mt_path(liftover=True))
        exome_mt = urc.get_ukbb_data(data_source, freeze, adj=True, split=True)

        exome_mt = remove_hard_filter_samples(data_source, freeze, exome_mt)
        logger.info(f"Count after removing hard filtered samples: {exome_mt.count()}")

        if args.interval_qc_filter:
            exome_mt = interval_qc_filter(
                data_source, freeze, exome_mt, args.pct_samples_20x
            )

        array_mt, exome_mt = prepare_array_and_exome_mt(
            data_source, freeze, array_mt, exome_mt, call_rate_cutoff, af_cutoff
        )

        exome_mt = exome_mt.checkpoint(
            urc.get_mt_checkpoint_path(
                data_source,
                freeze,
                name=f"exome_subset_concordance_callrate_{call_rate_cutoff}_af_{af_cutoff}"
            ),
            overwrite=overwrite
        )
        array_mt = array_mt.checkpoint(
            urc.get_mt_checkpoint_path(
                data_source,
                freeze,
                name=f"array_subset_concordance_callrate_{call_rate_cutoff}_af_{af_cutoff}"
            ),
            overwrite=overwrite
        )

        samples, variants = get_array_exome_concordance(array_mt, exome_mt)

        variants = variants.annotate_globals(
            callrate_cutoff=call_rate_cutoff, af_cutoff=af_cutoff
        )
        samples = samples.annotate_globals(
            callrate_cutoff=call_rate_cutoff, af_cutoff=af_cutoff
        )

        variants.write(
            urc.array_variant_concordance_path(data_source, freeze), overwrite=overwrite
        )
        samples.write(
            urc.array_sample_concordance_path(data_source, freeze), overwrite=overwrite
        )

    if args.make_variant_qc_ht:
        variants = hl.read_table(urc.array_variant_concordance_path(data_source, freeze))
        callrate_cutoff = variants.callrate_cutoff.take(1)[0]
        af_cutoff = variants.af_cutoff.take(1)[0]

        exome_mt = hl.read_matrix_table(
            urc.get_mt_checkpoint_path(
                data_source,
                freeze,
                name=f"exome_subset_concordance_callrate_{callrate_cutoff}_af_{af_cutoff}"
            )
        )

        variants = variants.annotate(AF=exome_mt.rows()[variants.key].variant_qc.AF[1])
        variants = variants.filter(
            (variants.prop_gt_con_non_ref > args.concordance_cutoff)
            & (variants.AF > args.variant_qc_af_cutoff)
        )
        variants = variants.repartition(1000)
        variants.write(
            f"{urc.variant_qc_prefix(data_source, freeze)}/array_variant_concordance_callrate_{callrate_cutoff}_non_ref_con_{args.concordance_cutoff}_AF_{args.variant_qc_af_cutoff}.ht",
            overwrite=overwrite
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true"
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
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
        action="store_true"
    )

    concordance = parser.add_argument_group("Compute array concordance with exomes")
    concordance.add_argument(
        "-s",
        "--data_source",
        help="Data source",
        choices=["regeneron", "broad"],
        default="broad"
    )
    concordance.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=urc.CURRENT_FREEZE, type=int
    )
    concordance.add_argument(
        "-c",
        "--array_concordance",
        help="Compute array concordance",
        action="store_true"
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
        action="store_true"
    )
    concordance.add_argument(
        "--pct_samples_20x",
        help="Percent samples at 20X to filter intervals",
        default=0.85,
        type=float
    )
    concordance.add_argument(
        "--make_variant_qc_ht",
        help="Make HT for concordance HT for variant QC",
        action="store_true"
    )
    concordance.add_argument(
        "--concordance_cutoff",
        help="Array exome concordance cutoff for variant QC HT.",
        type=float,
        default=0.9
    )
    concordance.add_argument(
        "--variant_qc_af_cutoff",
        help="Allele frequency cutoff used for variant QC HT creation, must be equal to or greater than af_cutoff.",
        type=float,
        default=0.001
    )

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)