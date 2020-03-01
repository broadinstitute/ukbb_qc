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
logger = logging.getLogger("rank_rf")
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


def create_rf_rank(data_source: str, freeze: int, run_hash: str) -> None:
    """
    Creates a ranked table for a RF run and writes it to its correct location in annotations.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param str run_hash: RF run hash
    :return: Nothing
    :rtype: None
    """
    # TODO: Should I add VQSR negative_train_site, default expects it ...
    logger.info(
        f"Annotating RF file for {data_source}.freeze_{freeze} RF run {run_hash} with quantile bins"
    )
    ht = ht.annotate(score=1 - ht.rf_probability["TP"])
    bin_ht = default_create_binned_ht(
        ht, n_bins, add_substrat={"interval": ht.interval_qc_pass}
    )

    ht.write(
        rf_path(data_source, freeze, "rf_result", run_hash=run_hash), overwrite=True
    )

    grouped_binned_ht = compute_grouped_binned_ht(
        binned_ht, checkpoint_path=f"gs://broad-ukbb-jgoodric/broad_5_vqsr_binned.ht",
    )
    # binned_ht.write(score_ranking_path(data_source, freeze, 'vqsr', binned=True), overwrite=args.overwrite)

    agg_ht = grouped_binned_ht.aggregate(
        **default_score_bin_agg(
            grouped_binned_ht,
            fam_stats_ht=hl.read_table(
                var_annotations_ht_path(data_source, freeze, "trio_stats")
            ),
        )
    )

    agg_ht.write(
        f"gs://broad-ukbb-jgoodric/broad_5_vqsr_binned.ht", overwrite=args.overwrite
    )


def create_vqsr_annotated_ht(data_source: str, freeze: int, vqsr_type: str):
    """
    Creates a rank table for VQSR and writes it to its correct location in annotations.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Nothing
    :rtype: None
    """
    rf_ht = hl.read_table(rf_annotated_path(data_source, freeze))
    ht = hl.read_table(
        var_annotations_ht_path(
            data_source, freeze, "vqsr" if vqsr_type == "AS" else "AS_TS_vqsr"
        )
    ).repartition(1000)
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
    ht = annotate_interval_qc_filter(data_source, freeze, ht)
    ht = ht.filter(ht.n_nonref > 0)
    ht.write(
        f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/{"vqsr" if vqsr_type == "AS" else "AS_TS_vqsr"}.ht',
        overwrite=True,
    )

    logger.info(f"Creating rank file for {data_source}.freeze_{freeze} VQSR")
    ht = hl.read_table(
        f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/{"vqsr" if vqsr_type == "AS" else "AS_TS_vqsr"}.ht'
    )

    ht.write(
        score_ranking_path(
            data_source, freeze, "vqsr" if vqsr_type == "AS" else "AS_TS_vqsr"
        ),
        overwrite=True,
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


def create_binned_data(
    ht: hl.Table, data: str, data_source: str, freeze: int, n_bins: int
) -> hl.Table:
    """
    Creates binned data from a rank Table grouped by rank_id (rank, biallelic, etc.), contig, snv, bi_allelic and singleton
    containing the information needed for evaluation plots.

    :param Table ht: Input rank table
    :param str data: Which data/run hash is being created
    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param int n_bins: Number of bins.
    :return: Binned Table
    :rtype: Table
    """

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
    ht = ht.annotate_globals(rank_variant_counts=bin_variant_counts)

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
    grouped_binned_ht = compute_grouped_binned_ht(ht, checkpoint_path=f"",)

    parent_ht = ht._parent
    agg_ht = grouped_binned_ht.aggregate(
        n_clinvar_path=hl.agg.count_where(parent_ht.clinvar_path),
        n_ukbb_array_con_common=hl.agg.count_where(parent_ht.ukbb_array_con_common),
        n_sib_singletons=hl.agg.count_where(parent_ht.sib_singletons),
        **default_score_bin_agg(grouped_binned_ht, fam_stats_ht=trio_stats_ht),
    )

    agg_ht.write(
        f"gs://broad-ukbb-jgoodric/broad_5_vqsr_binned.ht", overwrite=args.overwrite
    )


def main(args):
    def run_data(
        rank_func: Callable[..., None],
        rank_func_args: List[Any],
        score_type: str,
        data_source: str,
        freeze: int,
    ) -> None:
        """
        Wrapper for running script actions on given data.

        :param callable rank_func: Function for creating ranking file
        :param list of Any rank_func_args: Arguments to pass to the ranking function
        :param str score_type: Score being processed
        :param str data_source: 'regeneron' or 'broad'
        :param int freeze: One of the data freezes
        :return: Nothing -- this runs the script actions
        :rtype:  None
        """
        print(score_type)
        if score_type in ["vqsr", "AS_TS_vqsr"]:
            rank_file_path = score_ranking_path(data_source, freeze, score_type)
        else:
            rank_file_path = rf_path(
                data_source, freeze, "rf_result", run_hash=args.run_hash
            )

        if args.create_rank_file:
            rank_func(*rank_func_args)
        if args.run_sanity_checks:
            run_sanity_checks(hl.read_table(rank_file_path))
        if args.create_binned_file:
            ht = hl.read_table(rank_file_path)
            binned_ht = create_binned_data(
                ht, score_type, data_source, freeze, args.n_bins
            )
            binned_ht.write(
                score_ranking_path(data_source, freeze, score_type, binned=True),
                overwrite=True,
            )

    hl.init(log="/create_rank.log")
    data_source = args.data_source
    freeze = args.freeze

    if args.run_hash:
        run_data(
            create_rf_rank,
            [data_source, freeze, args.run_hash],
            args.run_hash,
            data_source,
            freeze,
        )

    if args.vqsr:
        run_data(
            create_vqsr_rank_ht,
            [data_source, freeze, args.vqsr_type],
            "vqsr" if args.vqsr_type == "AS" else "AS_TS_vqsr",
            data_source,
            freeze,
        )


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
        "--create_rank_file",
        help="When set, creates ranking file.",
        action="store_true",
    )
    parser.add_argument(
        "--run_sanity_checks",
        help="When set, runs ranking sanity checks.",
        action="store_true",
    )
    parser.add_argument(
        "--create_binned_file",
        help="When set, creates binned ranked file.",
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
        not args.create_rank_file
        and not args.run_sanity_checks
        and not args.create_binned_file
    ):
        sys.exit(
            "Error: At least one of --create_rank_file, --run_sanity_checks or --create_binned_file must be specified."
        )

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
