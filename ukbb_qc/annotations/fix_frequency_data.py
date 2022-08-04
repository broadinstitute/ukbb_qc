"""
This script generates allele frequency data for all high quality European samples (including relateds).

The previous frequency results that included relatednes ("cohort" frequency)
did not filter duplicate samples. The previous results also had some slightly deflated ANs
due to an accidental filtering of star (*) alleles pre-densification.

This script is intended to fix two things:
    - Remove double counting (from duplicate samples) 
    - Correct ANs (by using the 455k VDS, which contains all variants)
"""
import argparse
import logging

import hail as hl

from gnomad.utils.annotations import annotate_adj, annotate_freq
from gnomad.utils.slack import slack_notifications

from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds

from ukbb_qc.resources.basics import logging_path, release_ht_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import meta_ht_path
from ukbb_qc.resources.variant_qc import var_annotations_ht_path
from ukbb_qc.slack_creds import slack_token


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("fix_frequency_data")
logger.setLevel(logging.INFO)


POP = "EUR"
"""
Population to include in frequency calculations.
"""

N_SAMPLES = 394841
"""
Total number of high quality EUR samples (minus duplicates) in 455k VDS.
"""

TRANCHE_DATA = ("broad", CURRENT_FREEZE)
"""
UKB tranche data (data source and data freeze number).
"""


def main(args):
    """
    Recalculate frequency on high quality, European samples in the 455k VDS.
    """
    hl.init(log="/fix_455_frequency.log", default_reference="GRCh38")

    try:
        logger.info("Reading in 455k sample meta...")
        meta_ht = hl.read_table(meta_ht_path(*TRANCHE_DATA))

        # Filter to high quality, EUR UKBB samples only, also remove control samples and duplicates
        meta_ht = meta_ht.filter(
            meta_ht.sample_filters.high_quality
            & ~meta_ht.sample_filters.duplicate
            & (meta_ht.pan_ancestry_meta.pop == POP)
        )

        logger.info(
            "Reading in v4 VDS and filtering to high quality EUR UKB samples..."
        )
        vds = get_gnomad_v4_vds(split=False, remove_hard_filtered_samples=False)
        vds = hl.vds.filter_samples(vds, meta_ht, remove_dead_alleles=True)
        var_mt = vds.variant_data
        ref_mt = vds.reference_data
        assert (
            var_mt.count_cols() == N_SAMPLES
        ), f"Number of samples is {var_mt.count_cols()} but expected {N_SAMPLES}!"

        logger.info("Adding het non-ref, sex, and population annotations...")
        # Adding a Boolean for whether a sample had a heterozygous non-reference genotype
        # Need to add this prior to splitting multiallelics to make sure these genotypes
        # are not adjusted by the homalt hotfix downstream
        assert "LGT" in var_mt.entry, "VDS must be unsplit!"
        var_mt = var_mt.annotate_entries(het_non_ref=var_mt.LGT.is_het_non_ref())
        vds = hl.vds.VariantDataset(ref_mt, var_mt)
        vds = hl.vds.split_multi(vds)

        logger.info("Densifying...")
        mt = hl.vds.to_dense_mt(vds)
        mt = mt.filter_rows(hl.len(mt.alleles) > 1)
        # Add adj annotation (to get raw and adj frequencies)
        mt = annotate_adj(mt)

        # Temporary hotfix for depletion of homozygous alternate genotypes
        logger.info(
            "Setting het genotypes at sites with >1% AF (using previous frequencies) and > 0.9 AB to homalt..."
        )

        # Load previously calculated allele frequencies to avoid an extra frequency calculation
        freq_ht = (
            hl.read_table(release_ht_path(*TRANCHE_DATA))
            .select_globals()
            .select("freq")
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

        logger.info("Calculating frequencies")
        logger.info("Generating frequency data...")
        mt = annotate_freq(mt)
        ht = mt.select_rows("freq").select_globals("freq_meta")
        ht = ht.naive_coalesce(args.n_partitions)
        ht.write(var_annotations_ht_path("ukb_freq_fix", *TRANCHE_DATA), args.overwrite)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(*TRANCHE_DATA))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--n_partitions",
        help="Desired number of partitions for output HTs",
        default=5000,
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
