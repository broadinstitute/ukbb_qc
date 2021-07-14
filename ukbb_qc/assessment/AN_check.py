import argparse
import logging

import hail as hl

from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import SPARSE_ENTRIES

from ukbb_qc.resources.basics import get_ukbb_data, logging_path, release_ht_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.variant_qc import info_ht_path
from ukbb_qc.slack_creds import slack_token


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("AN_check")
logger.setLevel(logging.INFO)


def main(args):
    """
    Script to check how many variants have AN impacted by bug in hardcalls code.

    Bug in hardcalls code filtered out star alleles, which removed ref block information for some samples.
    This means that AN is lower/AF is higher than should be for a subset of variants.
    """

    hl.init(log="/AN_check.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze

    try:
        logger.info("Reading in raw MT...")
        mt = get_ukbb_data(
            data_source,
            freeze,
            key_by_locus_and_alleles=False,
            split=False,
            raw=True,
            repartition=False,
            meta_root="meta",
        ).select_entries(*SPARSE_ENTRIES)

        logger.info("Filtering raw MT to chr21...")
        mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chr21")])

        logger.info("Densifying...")
        mt = hl.experimental.densify(mt)
        mt = mt.filter_rows(hl.len(mt.alleles) > 1)

        logger.info("Adding het_non_ref annotation...")
        mt = mt.annotate_entries(het_non_ref=mt.LGT.is_het_non_ref())

        logger.info("Splitting densified MT...")
        mt = hl.experimental.sparse_split_multi(mt)

        logger.info("Applying homalt hotfix (het nonref aware)...")
        freq_ht = hl.read_table(release_ht_path(data_source="broad", freeze=6)).select(
            "freq"
        )
        freq_ht = freq_ht.select(AF=freq_ht.freq[0].AF)
        mt = mt.annotate_entries(
            GT=hl.if_else(
                mt.GT.is_het()
                # Skip adjusting genotypes if sample originally had a het nonref genotype
                & ~mt.het_non_ref
                & (freq_ht[mt.row_key].AF > 0.01)
                & (mt.AD[1] / mt.DP > 0.9),
                hl.call(1, 1),
                mt.GT,
            )
        )

        logger.info("Removing low QUAL variants and star alleles...")
        info_ht = hl.read_table(info_ht_path(data_source, freeze))
        mt = mt.filter_rows(
            (~info_ht[mt.row_key].AS_lowqual)
            & ((hl.len(mt.alleles) > 1) & (mt.alleles[1] != "*"))
        )

        logger.info("Filtering to unrelated high quality samples with defined batch...")
        mt = mt.filter_cols(mt.meta.sample_filters.high_quality)
        mt = mt.filter_cols(~mt.meta.sample_filters.related)
        mt = mt.filter_cols(hl.is_defined(mt.meta.ukbb_meta.batch))
        mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

        logger.info("Adding call stats and subsetting to rows...")
        mt = mt.annotate_rows(new_freq=hl.agg.call_stats(mt.GT, mt.alleles))
        ht = mt.rows()

        logger.info("Adding previously calculated 455K frequencies frequencies...")
        release_ht = (
            hl.read_table(
                "gs://broad-ukbb/broad.freeze_7/temp/broad.freeze_7.release.sites.ht"
            )
            .select_globals()
            .select("freq")
        )
        ht = ht.select(new_AN=ht.new_freq.AN, prev_AN=release_ht[ht.key].freq[1].AN)
        ht.aggregate(
            hl.struct(
                new_AN_neq_prev_AN=hl.agg.count_where(ht.new_AN != ht.prev_AN),
                new_AN_def=hl.agg.count_where(hl.is_defined(ht.new_AN)),
                prev_AN_def=hl.agg.count_where(hl.is_defined(ht.prev_AN)),
                prev_AN_def_and_neq_new_AN=hl.agg.count_where(
                    (hl.is_defined(ht.prev_AN)) & (ht.new_AN != ht.prev_AN)
                ),
            )
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
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
