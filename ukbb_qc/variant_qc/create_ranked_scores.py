import argparse
import logging
import sys
from pprint import pformat
import hail as hl
import gnomad_hail.resources.grch38 as grch38_resources
from gnomad_hail.utils.slack import try_slack
from gnomad_hail.utils.variant_qc import default_create_binned_ht, compute_grouped_binned_ht, default_score_bin_agg
from ukbb_qc.utils import annotate_interval_qc_filter
from ukbb_qc.resources.variant_qc import score_ranking_path, rf_annotated_path, rf_path, var_annotations_ht_path, get_checkpoint_path, CURRENT_FREEZE


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("quantile_bin_rf")
logger.setLevel(logging.INFO)


def run_sanity_checks(ht: hl.Table) -> None:
    """
    Runs and prints sanity checks on rank table.

    :param Table ht: input ranks Table
    :return: Nothing
    :rtype: None
    """
    print(
        ht.aggregate(
            hl.struct(
                was_split=hl.agg.counter(ht.was_split),
                has_biallelic_rank=hl.agg.counter(hl.is_defined(ht.biallelic_rank)),
                was_singleton=hl.agg.counter(ht.singleton),
                has_singleton_rank=hl.agg.counter(hl.is_defined(ht.singleton_rank)),
                was_split_singleton=hl.agg.counter(ht.singleton & ~ht.was_split),
                has_biallelic_singleton_rank=hl.agg.counter(
                    hl.is_defined(ht.biallelic_singleton_rank)
                ),
            )
        )
    )


def create_quantile_bin_ht(data_source: str, freeze: int, metric: str, n_bins: int, overwrite: bool) -> None:
    """
    Creates a table with quantile bin annotations added for a RF run and writes it to its correct location in annotations.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param str metric: Which data/run hash is being created
    :param n_bins: Number of bins to bin the data into
    :param bool overwrite: Should output files be overwriten if present
    :return: Nothing
    :rtype: None
    """
    # TODO: Should I add VQSR negative_train_site, default expects it ...
    logger.info(
        f"Annotating RF file for {data_source}.freeze_{freeze} RF run {metric} with quantile bins"
    )
    if metric.endswith("vqsr"):
        rf_ht = hl.read_table(rf_annotated_path(data_source, freeze))
        ht = hl.read_table(
            var_annotations_ht_path(
                data_source,
                freeze,
                metric
            )
        ).repartition(1000)

        # TODO: Should filter sites filtered by Laurent's lowqual, is this info in the info HT ?
        ht = ht.filter(~ht.filters.contains("LowQual"))
        ht = ht.annotate(
            AS_VQSLOD=hl.float(ht.info.AS_VQSLOD[ht.a_index - 1]),
            AS_culprit=ht.info.AS_culprit[ht.a_index - 1],
        )

        # TODO: NEGATIVE_TRAIN_SITE & POSITIVE_TRAIN_SITE is there an AS?
        ht = ht.annotate(
            **rf_ht[ht.key],
            score=ht.AS_VQSLOD,
            negative_train_site=ht.info.NEGATIVE_TRAIN_SITE,
            positive_train_site=ht.info.POSITIVE_TRAIN_SITE,
            culprit=ht.AS_culprit,
        )

        # TODO: Do we want to use the tranche 2 interval QC
        ht = annotate_interval_qc_filter(data_source, freeze, ht)
        ht = ht.filter(ht.n_nonref > 0)
    else:
        ht = hl.read_table(rf_path(data_source, freeze, 'rf_result', run_hash=metric))
        # TODO: Might need to rethink the score annotation 1-rf_probability, but default expects score annotation
        ht = ht.annotate(
            score=1-ht.rf_probability['TP']
        )

    # TODO: Should I add VQSR negative_train_site, default expects it ...
    bin_ht = default_create_binned_ht(
        ht,
        n_bins,
        add_substrat={"interval": ht.interval_qc_pass}
    )

    bin_ht.write(
        score_ranking_path(
            data_source,
            freeze,
            metric
        ),
        overwrite=overwrite,
    )


def filter_clinvar_path():
    clinvar_ht = grch38_resources.reference_data.clinvar.ht()
    logger.info(f"Found {clinvar_ht.count()} variants in Clinvar before filtering")
    no_star_assertions = hl.literal(
        {
            "no_assertion_provided",
            "no_assertion_criteria_provided",
            "no_interpretation_for_the_single_variant",
        }
    )
    clinvar_ht = clinvar_ht.filter(
        hl.set(clinvar_ht.info.CLNREVSTAT).intersection(no_star_assertions).length()
        > 0,
        keep=False,
    )
    logger.info(
        f"Found {clinvar_ht.count()} variants after removing variants without assertions"
    )
    clinvar_ht = clinvar_ht.filter(
        clinvar_ht.info.CLNSIG.map(lambda x: x.lower())
        .map(lambda x: x.contains("pathogenic"))
        .any(lambda x: x),
        keep=True,
    )
    logger.info(
        f"Found {clinvar_ht.count()} variants after filtering to (likely) pathogenic variants"
    )
    clinvar_ht = clinvar_ht.filter(
        hl.is_defined(clinvar_ht.info.CLNSIGCONF), keep=False
    )
    logger.info(
        f"Found {clinvar_ht.count()} variants after filtering to variants without CLNSIGCONF (conflicting clinical interpretations)"
    )

    return clinvar_ht


def create_grouped_bin_ht(
    data_source: str, freeze: int, metric: str, overwrite: bool
) -> None:
    """
    Creates binned data from a quantile bin annotated Table grouped by bin_id (rank, biallelic, etc.), contig, snv,
    bi_allelic and singleton containing the information needed for evaluation plots.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param str metric: Which data/run hash is being created
    :param bool overwrite: Should output files be overwriten if present
    :return: Binned Table
    :rtype: None
    """

    ht = hl.read_table(
        score_ranking_path(
            data_source,
            freeze,
            metric
        ),
        overwrite=overwrite,
    )

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
    bin_variant_counts = ht.aggregate(hl.Struct(**count_expr))
    logger.info(f"Found the following variant counts:\n {pformat(bin_variant_counts)}")
    ht = ht.annotate_globals(bin_variant_counts=bin_variant_counts)

    # Load ClinVar path data
    if not hl.utils.hadoop_exists(
        "gs://broad-ukbb/resources/clinvar_20190609.pathogenic.ht"
    ):
        clinvar_path_ht = filter_clinvar_path()
        clinvar_path_ht.write(
            "gs://broad-ukbb/resources/clinvar_20190609.pathogenic.ht"
        )

    clinvar_path_ht = filter_clinvar_path()
    ht = ht.annotate(clinvar_path=hl.is_defined(clinvar_path_ht[ht.key]))
    trio_stats_ht = hl.read_table(
        var_annotations_ht_path(data_source, freeze, "trio_stats")
    )

    logger.info(f"Creating binned rank table")
    grouped_binned_ht = compute_grouped_binned_ht(
        ht,
        checkpoint_path=get_checkpoint_path(data_source, freeze, f'grouped_bin_{metric}'),
    )

    parent_ht = ht._parent
    agg_ht = grouped_binned_ht.aggregate(
        n_clinvar_path=hl.agg.count_where(parent_ht.clinvar_path),
        n_ukbb_array_con_common=hl.agg.count_where(parent_ht.ukbb_array_con_common),
        n_sib_singletons=hl.agg.count_where(parent_ht.sib_singletons),
        **default_score_bin_agg(grouped_binned_ht, fam_stats_ht=trio_stats_ht),
    )

    agg_ht.write(score_ranking_path(data_source, freeze, metric, binned=True), overwrite=overwrite)


def main(args):
    hl.init(log="/create_quantile_rank_bins.log")

    if args.vqsr:
        metric = "vqsr" if args.vqsr_type == "AS" else "AS_TS_vqsr"
    else:
        metric = args.run_hash

    if args.create_quantile_bin_ht:
        create_quantile_bin_ht(args.data_source, args.freeze, metric, args.overwrite)
    if args.run_sanity_checks:
        run_sanity_checks(hl.read_table(score_ranking_path(args.data_source, args.freeze, metric)))
    if args.create_aggregated_bin_ht:
        create_grouped_bin_ht(args.data_source, args.freeze, metric, args.overwrite)


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
        "-s",
        "--data_source",
        help="Data source",
        choices=["regeneron", "broad"],
        default="broad",
    )
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )

    parser.add_argument("--run_hash", help="Run hash for RF results to be ranked.")
    parser.add_argument(
        "--vqsr", help="When set, creates the VQSR rank file.", action="store_true"
    )
    parser.add_argument(
        "--vqsr_type",
        help="What type of VQSR was run: allele-specific, or allele-specific with transmitted singletons",
        type=str,
        choices=["AS", "AS_TS"],
        default="AS",
    )
    parser.add_argument(
        "--create_quantile_bin_ht",
        help="When set, creates file annotated with quantile bin based on vqsr/ RF run hash score.",
        action="store_true",
    )
    parser.add_argument(
        "--run_sanity_checks",
        help="When set, runs ranking sanity checks.",
        action="store_true",
    )
    parser.add_argument(
        "--create_aggregated_bin_ht",
        help="When set, creates a file with aggregate counts of variants based on quantile bins.",
        action="store_true",
    )
    parser.add_argument(
        "--n_bins",
        help="Number of bins for the binned file (default: 100)",
        default=100,
        type=int,
    )
    args = parser.parse_args()

    if (
        not args.create_quantile_bin_ht
        and not args.run_sanity_checks
        and not args.create_aggregated_bin_ht
    ):
        sys.exit(
            "Error: At least one of --create_quantile_bin_ht, --run_sanity_checks or --create_aggregated_bin_ht must be specified."
        )

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
