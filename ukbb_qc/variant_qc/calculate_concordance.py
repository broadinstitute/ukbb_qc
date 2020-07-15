import argparse
import logging
import hail as hl
from gnomad.utils.filtering import filter_low_conf_regions
from gnomad.utils.slack import slack_notifications
from gnomad.variant_qc.evaluation import compute_quantile_bin, create_truth_sample_ht
from ukbb_qc.resources.basics import capture_ht_path, get_checkpoint_path, get_ukbb_data
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.variant_qc import (
    binned_concordance_path,
    get_truth_sample_data,
    info_ht_path,
    score_quantile_bin_path,
    TRUTH_SAMPLES,
    TRUTH_SAMPLE_NAMES,
    truth_sample_mt_path,
    var_annotations_ht_path,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("calculate_concordance")
logger.setLevel(logging.INFO)


def compute_binned_truth_sample_concordance(
    ht: hl.Table, binned_score_ht: hl.Table, n_bins: int = 100
) -> hl.Table:
    """
    Determines the concordance (TP, FP, FN) between a truth sample within the callset and the samples truth data
    grouped by bins computed using `compute_quantile_bin`.

    .. note::
        The input 'ht` should contain three row fields:
            - score: value to use for quantile binning
            - GT: a CallExpression containing the genotype of the evaluation data for the sample
            - truth_GT: a CallExpression containing the genotype of the truth sample
        The input `binned_score_ht` should contain:
             - score: value used to bin the full callset
             - bin: the full callset quantile bin

    The table is grouped by global/truth sample bin and variant type and contains TP, FP and FN.

    :param ht: Input HT
    :param binned_score_ht: Table with the an annotation for quantile bin for each variant
    :param n_bins: Number of bins to bin the data into
    :return: Binned truth sample concordance HT
    """
    # Annotate score and global bin
    indexed_binned_score_ht = binned_score_ht[ht.key]
    ht = ht.annotate(
        score=indexed_binned_score_ht.score,
        snv=hl.is_snp(ht.alleles[0], ht.alleles[1]),
        global_bin=indexed_binned_score_ht.bin,
        global_interval_bin=indexed_binned_score_ht.interval_bin,
    )

    # Annotate the truth sample quantile bin
    bin_ht = compute_quantile_bin(
        ht,
        score_expr=ht.score,
        bin_expr={
            "truth_sample_bin": hl.expr.bool(True),
            "truth_sample_interval_bin": ht.interval_qc_pass,
        },
        n_bins=n_bins,
    )
    ht = ht.join(bin_ht, how="left")

    # Explode the global and truth sample bins
    ht = ht.annotate(
        bin=[
            hl.tuple(["global_bin", ht.global_bin]),
            hl.tuple(["truth_sample_bin", ht.truth_sample_bin]),
            hl.tuple(["global_interval_bin", ht.global_interval_bin]),
            hl.tuple(["truth_sample_interval_bin", ht.truth_sample_interval_bin]),
        ]
    )

    ht = ht.explode(ht.bin)
    ht = ht.annotate(bin_id=ht.bin[0], bin=hl.int(ht.bin[1]))

    # Compute TP, FP and FN by bin_id, variant type and bin
    return (
        ht.group_by("bin_id", "snv", "bin")
        .aggregate(
            # TP => allele is found in both data sets
            tp=hl.agg.count_where(ht.GT.is_non_ref() & ht.truth_GT.is_non_ref()),
            # FP => allele is found only in test data set
            fp=hl.agg.count_where(
                ht.GT.is_non_ref() & hl.or_else(ht.truth_GT.is_hom_ref(), True)
            ),
            # FN => allele is found in truth data only
            fn=hl.agg.count_where(
                ht.GT.is_hom_ref() & hl.or_else(ht.truth_GT.is_non_ref(), True)
            ),
            min_score=hl.agg.min(ht.score),
            max_score=hl.agg.max(ht.score),
            n_alleles=hl.agg.count(),
        )
        .repartition(5)
    )


def main(args):

    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)

    info_ht = hl.read_table(info_ht_path(*tranche_data))

    if args.extract_truth_samples:
        logger.info(f"Extracting truth samples from adj filtered hardcall MT...")
        mt = get_ukbb_data(*tranche_data, adj=True)

        mt = mt.filter_cols(hl.literal(TRUTH_SAMPLES).contains(mt.s))

        # Checkpoint to prevent needing to go through the large table a second time
        mt = mt.checkpoint(
            get_checkpoint_path(*tranche_data, "truth_samples", mt=True),
            overwrite=args.overwrite,
        )

        for truth_sample in TRUTH_SAMPLE_NAMES:
            truth_sample_mt = mt.filter_cols(mt.s == get_truth_sample_data(*tranche_data, truth_sample, "s"))
            # Filter to variants in truth data
            truth_sample_mt = truth_sample_mt.filter_rows(hl.agg.any(truth_sample_mt.GT.is_non_ref()))
            truth_sample_mt.naive_coalesce(args.n_partitions).write(
                truth_sample_mt_path(truth_sample, *tranche_data), overwrite=args.overwrite
            )

    if args.merge_with_truth_data:
        for truth_sample in TRUTH_SAMPLE_NAMES:
            logger.info(f"Creating a merged table with callset truth sample and truth data for {truth_sample}...")

            # Load truth data
            mt = get_truth_sample_data(*tranche_data, truth_sample, "callset_truth_mt")
            truth_hc_intervals = get_truth_sample_data(*tranche_data, truth_sample, "hc_intervals")
            truth_mt = get_truth_sample_data(*tranche_data, truth_sample, "truth_mt")
            truth_mt = truth_mt.key_cols_by(s=hl.str(get_truth_sample_data(*tranche_data, truth_sample, "s")))

            # remove low quality sites
            mt = mt.filter_rows(~info_ht[mt.row_key].lowqual)

            ht = create_truth_sample_ht(mt, truth_mt, truth_hc_intervals)
            ht.write(
                var_annotations_ht_path(truth_sample, *tranche_data),
                overwrite=args.overwrite,
            )

    if args.bin_truth_sample_concordance:
        metrics = (
            []
            if not args.run_hash
            else [args.run_hash]
            if isinstance(args.run_hash, str)
            else args.run_hash
        )

        if args.vqsr:
            metrics.append("vqsr" if args.vqsr_type == "AS" else "AS_TS_vqsr")

        for truth_sample in TRUTH_SAMPLE_NAMES:
            for metric in metrics:
                logger.info(
                    f"Creating binned concordance table for {truth_sample} for metric {metric}"
                )
                ht = hl.read_table(
                    var_annotations_ht_path(truth_sample, *tranche_data),
                )

                logger.info(
                    "Filtering out low confidence regions and sites outside the calling intervals..."
                )
                ht = filter_low_conf_regions(
                    ht,
                    filter_lcr=True,
                    # TODO: Uncomment when we have segdup path
                    filter_decoy=False,  # True,
                    filter_segdup=False,  # True,
                )
                calling_intervals = hl.read_table(capture_ht_path(data_source))
                ht = ht.filter(hl.is_defined(calling_intervals[ht.locus]))

                logger.info("Loading HT containing RF or VQSR scores annotated with a bin based on metric quantiles...")
                metric_ht = hl.read_table(
                    score_quantile_bin_path(metric, *tranche_data)
                )
                ht = ht.filter(hl.is_defined(metric_ht[ht.key]))

                ht = ht.annotate(
                    score=metric_ht[ht.key].score,
                    interval_qc_pass=metric_ht[ht.key].interval_qc_pass,
                )

                # TODO: add interval_global_rank and interval_truth_sample_rank using ht.interval_qc_pass
                ht = compute_binned_truth_sample_concordance(
                    ht, metric_ht, args.n_bins,
                )
                ht.write(
                    binned_concordance_path(truth_sample, metric, *tranche_data),
                    overwrite=args.overwrite,
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--extract_truth_samples",
        help="Extract truth samples from matrix table",
        action="store_true",
    )
    parser.add_argument(
        "--n_partitions",
        help="Desired number of partitions for output Table/MatrixTable",
        default=5000,
        type=int,
    )
    parser.add_argument(
        "--merge_with_truth_data",
        help="Computes a table for each truth sample comparing the truth sample in the callset vs the truth.",
        action="store_true",
    )
    parser.add_argument(
        "--bin_truth_sample_concordance",
        help="Merges individual concordance results with specified metric binned files.",
        action="store_true",
    )
    parser.add_argument("--run_hash", help="RF hash(es) for annotation.", nargs="+")
    parser.add_argument(
        "--vqsr", help="When set, annotates with VQSR rank file.", action="store_true"
    )
    parser.add_argument(
        "--vqsr_type",
        help="What type of VQSR was run: allele-specific, or allele-specific with transmitted singletons",
        type=str,
        choices=["AS", "AS_TS"],
        default="AS",
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
