import argparse
import logging
from pprint import pformat

import hail as hl

from gnomad.resources.grch38.reference_data import clinvar
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications
from gnomad.variant_qc.evaluation import compute_grouped_binned_ht
from gnomad.variant_qc.pipeline import create_binned_ht, score_bin_agg
from ukbb_qc.load_data.utils import load_clinvar_path
from ukbb_qc.resources.basics import get_checkpoint_path, logging_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.variant_qc import (
    clinvar_pathogenic_ht_path,
    info_ht_path,
    rf_annotated_path,
    rf_path,
    score_bin_path,
    var_annotations_ht_path,
)
from ukbb_qc.slack_creds import slack_token
from ukbb_qc.utils.utils import vqsr_run_check


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("bin_rf")
logger.setLevel(logging.INFO)


def create_rank_bin_ht(
    data_source: str, freeze: int, metric: str, n_bins: int, overwrite: bool = False
) -> None:
    """
    Creates a table with bin annotations added for a RF run and writes it to its correct location in annotations.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param str metric: Which data/run hash is being created
    :param n_bins: Number of bins to bin the data into
    :param bool overwrite: Should output files be overwritten if present
    :return: Nothing
    :rtype: None
    """
    logger.info(f"Annotating variant QC results for {metric} with bins")
    info_ht = hl.read_table(info_ht_path(data_source, freeze))
    if metric.endswith("vqsr"):
        rf_ht = hl.read_table(rf_annotated_path(data_source, freeze))
        ht = hl.read_table(var_annotations_ht_path(metric, data_source, freeze))

        ht = ht.filter(
            ~info_ht[ht.key].AS_lowqual
            & ~hl.is_infinite(ht.info.AS_VQSLOD)
            & hl.is_defined(rf_ht[ht.key])
        )
        ht = ht.annotate(
            **rf_ht[ht.key],
            score=ht.info.AS_VQSLOD,
            negative_train_site=ht.info.NEGATIVE_TRAIN_SITE,
            positive_train_site=ht.info.POSITIVE_TRAIN_SITE,
            AS_culprit=ht.info.AS_culprit,
        )

    else:
        ht = hl.read_table(rf_path(data_source, freeze, "rf_result", run_hash=metric))
        ht = ht.annotate(
            info=info_ht[ht.key].info,
            positive_train_site=ht.tp,
            negative_train_site=ht.fp,
            score=ht.rf_probability["TP"],
        )

    bin_ht = create_binned_ht(
        ht, n_bins, add_substrat={"interval": ht.interval_qc_pass}
    )

    bin_ht.write(
        score_bin_path(metric, data_source, freeze), overwrite=overwrite,
    )


def create_grouped_bin_ht(
    data_source: str, freeze: int, metric: str, overwrite: bool
) -> None:
    """
    Creates binned data from a bin annotated Table grouped by bin_id (rank, bi-allelic, etc.), contig, snv,
    bi_allelic and singleton containing the information needed for evaluation plots.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param str metric: Which data/run hash is being created
    :param bool overwrite: Should output files be overwritten if present
    :return: None
    :rtype: None
    """

    ht = hl.read_table(score_bin_path(metric, data_source, freeze))

    # Count variants for ranking
    count_expr = {
        x: hl.agg.filter(
            hl.is_defined(ht[x]),
            hl.agg.counter(
                hl.cond(hl.is_snp(ht.alleles[0], ht.alleles[1]), "snv", "indel")
            ),
        )
        for x in ht.row
        if x.endswith("bin")
    }
    bin_variant_counts = ht.aggregate(hl.struct(**count_expr))
    logger.info(f"Found the following variant counts:\n {pformat(bin_variant_counts)}")
    ht = ht.annotate_globals(bin_variant_counts=bin_variant_counts)

    # Load ClinVar pathogenic data
    clinvar_version = clinvar.default_version
    if not file_exists(clinvar_pathogenic_ht_path(clinvar_version)):
        clinvar_path_ht = load_clinvar_path()
        clinvar_path_ht.write(clinvar_pathogenic_ht_path(clinvar_version))

    clinvar_path_ht = hl.read_table(clinvar_pathogenic_ht_path(clinvar_version))
    ht = ht.annotate(clinvar_path=hl.is_defined(clinvar_path_ht[ht.key]))
    trio_stats_ht = hl.read_table(
        var_annotations_ht_path("trio_stats", data_source, freeze)
    )

    logger.info(f"Creating grouped bin table...")
    grouped_binned_ht = compute_grouped_binned_ht(
        ht,
        checkpoint_path=get_checkpoint_path(
            data_source, freeze, f"grouped_bin_{metric}"
        ),
    )

    logger.info(f"Aggregating grouped bin table...")
    parent_ht = grouped_binned_ht._parent
    agg_ht = grouped_binned_ht.aggregate(
        n_clinvar_path=hl.agg.count_where(parent_ht.clinvar_path),
        n_ukbb_array_con_common=hl.agg.count_where(parent_ht.ukbb_array_con_common),
        n_sib_singletons=hl.agg.count_where(parent_ht.sibling_singleton),
        **score_bin_agg(grouped_binned_ht, fam_stats_ht=trio_stats_ht),
    )

    agg_ht.write(
        score_bin_path(metric, data_source, freeze, aggregated=True),
        overwrite=overwrite,
    )


def main(args):
    hl.init(log="/create_rank_bins.log")

    data_source = "broad"
    freeze = args.freeze
    if args.vqsr:
        vqsr_run_check(data_source, freeze, args.vqsr_type)
        metric = "vqsr" if args.vqsr_type == "AS" else f"{args.vqsr_type}_vqsr"
    else:
        metric = args.run_hash

    try:

        if args.create_rank_bin_ht:
            create_rank_bin_ht(data_source, freeze, metric, args.n_bins, args.overwrite)
        if args.run_sanity_checks:
            ht = hl.read_table(score_bin_path(metric, data_source, freeze))
            logger.info("Running sanity checks...")
            print(
                ht.aggregate(
                    hl.struct(
                        was_split=hl.agg.counter(ht.was_split),
                        has_biallelic_rank=hl.agg.counter(
                            hl.is_defined(ht.biallelic_bin)
                        ),
                        was_singleton=hl.agg.counter(ht.singleton),
                        has_singleton_rank=hl.agg.counter(
                            hl.is_defined(ht.singleton_bin)
                        ),
                        was_split_singleton=hl.agg.counter(
                            ht.singleton & ~ht.was_split
                        ),
                        has_biallelic_singleton_rank=hl.agg.counter(
                            hl.is_defined(ht.biallelic_singleton_bin)
                        ),
                    )
                )
            )
        # Note: Use only workers, it typically crashes with premptibles
        if args.create_aggregated_bin_ht:
            create_grouped_bin_ht(data_source, freeze, metric, args.overwrite)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(data_source, freeze))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument("--run_hash", help="Run hash for RF results to be ranked.")
    parser.add_argument(
        "--vqsr", help="When set, creates the VQSR rank file.", action="store_true"
    )
    parser.add_argument(
        "--vqsr_type", help="What type of VQSR was run", type=str, default="AS",
    )
    parser.add_argument(
        "--create_rank_bin_ht",
        help="When set, creates a variant HT annotated with ranked bins of approximately equal size based on vqsr/RF run hash score.",
        action="store_true",
    )
    parser.add_argument(
        "--run_sanity_checks",
        help="When set, runs ranking sanity checks.",
        action="store_true",
    )
    parser.add_argument(
        "--create_aggregated_bin_ht",
        help="When set, creates a file with aggregate counts of variants based on ranked bins.",
        action="store_true",
    )
    parser.add_argument(
        "--n_bins",
        help="Number of bins for the binned file (default: 100)",
        default=100,
        type=int,
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
