import argparse
import logging

import hail as hl

from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import SPARSE_ENTRIES

from ukbb_qc.resources.basics import get_checkpoint_path, get_ukbb_data, release_ht_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.variant_qc import info_ht_path
from ukbb_qc.slack_creds import slack_token


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("get_impacted_variants")
logger.setLevel(logging.INFO)


def get_het_non_ref_impacted_var(
    mt: hl.MatrixTable, info_ht: hl.Table, freq_ht: hl.Table
) -> None:
    """
    Filter to variants where homalt hotfix incorrectly adjusts het nonref genotype calls.

    Code for this function was copied from this notebook: gs://broad-ukbb/broad.freeze_7/notebooks/ukbb_multi_allelic_check.ipynb

    :param hl.MatrixTable mt: Raw, split MatrixTable annotated with original genotype het nonref status.
    :param hl.Table info_ht: Info Hail Table containing AS_lowqual information.
    :param hl.Table freq_ht: Hail Table containing 455K frequency information.
    :return: None 
    """
    # NOTE: This step is redundant: we remove low qual variants from the release sites HT,
    # and this code is pulling frequencies from the release sites HT
    # Kept in the code just to show how it was run
    logger.info("Removing low qual variants and star alleles...")
    mt = mt.filter_rows(
        (~info_ht[mt.row_key].AS_lowqual)
        & ((hl.len(mt.alleles) > 1) & (mt.alleles[1] != "*"))
    )

    logger.info("Filtering to common (AF > 0.01) variants...")
    mt = mt.annotate_rows(AF=freq_ht[mt.row_key].freq[0].AF)
    mt = mt.filter_rows(mt.AF > 0.01)

    logger.info(
        "Filtering to variants with at least one het nonref call and checkpointing..."
    )
    mt = mt.filter_rows(hl.agg.sum(mt.het_non_ref) >= 1)
    # NOTE: This MT was written on 5/2/21 in ukbb_multi_allelic_check.ipynb
    mt = mt.checkpoint("gs://broad-ukbb/broad.freeze_7/temp/temp_het_nonref_sites.mt")

    logger.info(
        "Adding new genotype entries with original homalt hotfix and het nonref-aware hotfix..."
    )
    mt = mt.annotate_entries(
        GT_hotfix=hl.if_else(
            mt.GT.is_het() & (mt.AF > 0.01) & (mt.AD[1] / mt.DP > 0.9),
            hl.call(1, 1),
            mt.GT,
        ),
        GT_hetnonref_fix=hl.if_else(
            mt.GT.is_het()
            & ~mt.het_non_ref
            & (mt.AF > 0.01)
            & (mt.AD[1] / mt.DP > 0.9),
            hl.call(1, 1),
            mt.GT,
        ),
    )
    logger.info("Checking homozygote counts using the different genotypes...")
    # Taking index 1 of gt_stats because it contains an entry for the reference allele at index 0
    # homozygote_count (tarray of tint32) -
    # Homozygote genotype counts for each allele, including the reference. Only diploid genotype calls are counted.
    # (from hail docs, https://hail.is/docs/0.2/aggregators.html#hail.expr.aggregators.call_stats)
    mt = mt.annotate_rows(gt_stats=hl.agg.call_stats(mt.GT, mt.alleles))
    mt = mt.transmute_rows(hom_count=mt.gt_stats.homozygote_count[1])

    mt = mt.annotate_rows(gt_stats=hl.agg.call_stats(mt.GT_hotfix, mt.alleles))
    mt = mt.transmute_rows(hom_count_hotfix=mt.gt_stats.homozygote_count[1])

    mt = mt.annotate_rows(gt_stats=hl.agg.call_stats(mt.GT_hetnonref_fix, mt.alleles))
    mt = mt.transmute_rows(hom_count_hetnonref_fix=mt.gt_stats.homozygote_count[1])

    logger.info(
        "Filtering to rows where homalt hotfix erroneously adjusts het nonref calls..."
    )
    mt = mt.filter_rows(mt.hom_count_hotfix != mt.hom_count_hetnonref_fix)
    # This MT was also written on 5/20/21 in ukbb_multi_allelic_check.ipynb
    # This MT has 56,444 rows
    mt.write("gs://broad-ukbb/broad.freeze_7/temp/het_nonref_sites.mt", overwrite=True)


def get_freq_impacted_var(mt: hl.MatrixTable, ht: hl.Table, ht_300K: hl.Table) -> None:
    """
    Filter to variants that have an AF > 0.01 only in the 300K tranche or only in the 455K tranche.

    :param hl.MatrixTable mt: Raw split MatrixTable.
    :param hl.Table ht: 455K release HT.
    :param hl.Table ht_300K: 300K release HT.
    :return: None
    """
    logger.info("Filtering both frequency HTs to AF > 0.01...")
    ht_300K = ht_300K.filter(ht_300K.freq[0].AF > 0.01)
    ht = ht.filter(ht.freq[0].AF > 0.01)

    logger.info("Getting variants that are only common in the 455K...")
    not_in_300k = ht.anti_join(ht_300K)

    logger.info("Getting variants that are only common in the 300K...")
    not_in_455k = ht_300K.anti_join(ht)

    mt = mt.annotate_rows(
        freq_455k=not_in_300k[mt.row_key].freq, freq_300k=not_in_455k[mt.row_key].freq
    )
    mt = mt.filter_rows(hl.is_defined(mt.freq_455k) | hl.is_defined(mt.freq_300k))
    # NOTE: Wrote this MT on 6/8/21 in gs://broad-ukbb/broad.freeze_7/notebooks/homalt_hotfix_counts.ipynb
    # MT written on 6/8/21 had 16529 rows and 454774 columns
    # I also need to overwrite this MT because I did not filter any samples (so this contains withdrawn samples and known dups)
    mt.write(
        "gs://broad-ukbb/broad.freeze_7/temp/homalt_hotfix_variants.mt", overwrite=True
    )


def main(args):
    """
    Script used to get variants impacted by homalt hotfix.

    Most of the work was done in notebooks, and any results created in notebooks are noted with comments.
    """
    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)

    hl.init(log="/get_impacted_variants.log", default_reference="GRCh38")

    logger.info("Reading in raw MT...")
    mt = get_ukbb_data(*tranche_data, raw=True).select_entries(*SPARSE_ENTRIES)

    logger.info(
        "Adding annotation for whether original genotype (pre-splitting multiallelics) is het nonref..."
    )
    mt = mt.annotate_entries(het_non_ref=mt.LGT.is_het_non_ref())

    logger.info("Splitting multiallelics...")
    mt = hl.experimental.sparse_split_multi(mt)

    logger.info("Reading in 455K release HT (to get frequency information)...")
    # TODO: Move this table back to release HT path
    # Overwrote HT at release HT path with frequencies recalculated from first patch release, but
    # we should keep the original 455K frequencies
    freq_ht = (
        hl.read_table(
            "gs://broad-ukbb/broad.freeze_7/temp/broad.freeze_7.release.sites.ht"
        )
        .select_globals()
        .select("freq")
    )

    if args.get_het_non_ref_impacted_var:
        logger.info("Reading in info HT...")
        info_ht = hl.read_table(info_ht_path(*tranche_data))

        logger.info(
            "Checking for variants with het nonref calls that were incorrectly adjusted with homalt hotfix..."
        )
        get_het_non_ref_impacted_var(mt, info_ht, freq_ht)

    if args.get_freq_impacted_var:
        # NOTE: We used the 300K frequencies to apply the homalt hotfix to the 455K tranche prior to calculating
        # frequencies on the 455K tranche to avoid an extra densify step
        logger.info("Reading in 300K release HT (to get 300K frequency information)...")
        freq_ht_300k = hl.read_table(release_ht_path(data_source, 6))

        logger.info(
            "Checking for variants that are common (AF > 0.01) in only the 300K or the 455K but not both datasets..."
        )
        get_freq_impacted_var(mt, freq_ht, freq_ht_300k)

    if args.join_var:
        # Adding _n_partitions to these reads so that the join operates on two HTs with the same number of partitions
        logger.info("Reading in het non ref impacted variants MT into HT...")
        het_nonref_var_ht = (
            hl.read_matrix_table(
                "gs://broad-ukbb/broad.freeze_7/temp/het_nonref_sites.mt",
                _n_partitions=args.n_partitions,
            )
            .rows()
            .select()
        )

        logger.info(
            "Reading in variants that are common in only one of the 300K and 455K data tranches..."
        )
        freq_var_ht = (
            hl.read_matrix_table(
                "gs://broad-ukbb/broad.freeze_7/temp/homalt_hotfix_variants.mt",
                _n_partitions=args.n_partitions,
            )
            .rows()
            .select()
        )

        ht = het_nonref_var_ht.join(freq_var_ht, how="outer")
        ht.write(
            get_checkpoint_path(*tranche_data, name="all_patch_release_variants"),
            overwrite=args.overwrite,
        )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--overwrite", help="Overwrite existing data", action="store_true"
    )
    parser.add_argument(
        "--get-het-non-ref-impacted-var",
        help="Filter raw MT to variants with het nonref genotypes impacted by homalt hotfix",
        action="store_true",
    )
    parser.add_argument(
        "--get-freq-impacted-var",
        help="Filter raw MT to variants with AF > 0.01 only in the 300K or only in the 455K",
        action="store_true",
    )
    parser.add_argument(
        "--join-var",
        help="Join het nonref and freq variants into a single hail Table",
        action="store_true",
    )
    parser.add_argument(
        "--n-partitions",
        help="Number of partitions to use for joint variants HT",
        default=1000,
        type=int,
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
