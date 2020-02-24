from gnomad_hail import *
from ukbb_qc.resources import *


def compute_concordance(
    mt: hl.MatrixTable, other_mt: hl.MatrixTable, name: str
) -> Tuple[hl.Table, hl.Table]:
    # Filter to sites present in mt samples
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    other_mt = other_mt.filter_rows(hl.agg.any(other_mt.GT.is_non_ref()))

    summary, sample_concordance_ht, sites_concordance_ht = hl.concordance(mt, other_mt)
    logger.info(f"{name} concordance summary: %s", pformat(summary))

    return sample_concordance_ht, sites_concordance_ht


def write_truth_concordance(
    data_source: str, freeze: int, truth_sample: str, overwrite: bool
) -> None:
    truth_samples = get_truth_sample_info(data_source, freeze)
    mt = hl.read_matrix_table(truth_samples[truth_sample]["mt"])
    mt = mt.filter_rows(mt.alleles[1] != "*")
    vqsr_ht = hl.read_table(var_annotations_ht_path(data_source, freeze, "vqsr"))
    mt = mt.filter_rows(~vqsr_ht[mt.row_key].filters.contains("LowQual"))

    mt = mt.select_rows(a_index=mt.a_index, was_split=mt.was_split)

    mt = unphase_mt(mt)

    truth_mt = hl.read_matrix_table(truth_samples[truth_sample]["truth_mt"])

    truth_mt = truth_mt.key_cols_by(s=hl.str(truth_samples[truth_sample]["s"]))

    calling_intervals = hl.read_table(capture_ht_path(data_source, freeze))
    truth_mt = truth_mt.filter_rows(hl.is_defined(calling_intervals[truth_mt.locus]))
    mt = mt.filter_rows(hl.is_defined(calling_intervals[mt.locus]))
    truth_mt = unphase_mt(hl.split_multi_hts(truth_mt, left_aligned=False))

    sample_concordance_ht, sites_concordance_ht = compute_concordance(
        mt, truth_mt, name=truth_sample
    )
    sites_concordance_ht.write(
        var_annotations_ht_path(data_source, freeze, f"{truth_sample}_concordance"),
        overwrite=overwrite,
    )
    sample_concordance_ht.write(
        sample_annotations_table_path(
            data_source, freeze, f"{truth_sample}_concordance"
        ),
        overwrite=overwrite,
    )


def create_binned_concordance(
    data_source: str,
    freeze: int,
    truth_sample: str,
    metric: str,
    nbins: int,
    overwrite: bool,
) -> None:
    """
    Creates and writes a concordance table binned by rank (both absolute and relative) for a given data type, truth sample and metric.

    :param str data_type: One 'exomes' or 'genomes'
    :param str truth_sample: Which truth sample concordance to load
    :param str metric: One of the evaluation metrics (or a RF hash)
    :param int nbins: Number of bins for the rank
    :param bool overwrite: Whether to overwrite existing table
    :return: Nothing -- just writes the table
    :rtype: None
    """
    truth_samples = get_truth_sample_info(data_source, freeze)

    if (
        hl.hadoop_exists(
            binned_concordance_path(data_source, freeze, truth_sample, metric)
            + "/_SUCCESS"
        )
        and not overwrite
    ):
        logger.warn(
            f"Skipping binned concordance creation as {binned_concordance_path(data_source, freeze, truth_sample, metric)} exists and overwrite=False"
        )
    else:
        ht = hl.read_table(
            var_annotations_ht_path(data_source, freeze, f"{truth_sample}_concordance")
        )
        # Remove 1bp indels for syndip as cannot be trusted
        if truth_sample == "syndip":
            ht = ht.filter(
                hl.is_indel(ht.alleles[0], ht.alleles[1])
                & (hl.abs(hl.len(ht.alleles[0]) - hl.len(ht.alleles[1])) == 1),
                keep=False,
            )

        high_conf_intervals = hl.import_bed(
            truth_samples[truth_sample]["bed"],
            reference_genome="GRCh38",
            skip_invalid_intervals=True,
        )
        lcr = get_lcr_intervals()
        # TODO: change this back when we have segdup path
        # segdup = hl.import_locus_intervals(segdup_intervals_path)
        ht = ht.filter(
            hl.is_defined(high_conf_intervals[ht.locus])
            & hl.is_missing(lcr[ht.locus])  # &
            # hl.is_missing(segdup[ht.locus])
        )

        if metric in ["vqsr", "AS_TS_vqsr"]:
            metric_ht = hl.read_table(score_ranking_path(data_source, freeze, metric))
        else:
            metric_ht = hl.read_table(
                rf_path(data_source, freeze, "rf_result", run_hash=metric)
            )

        metric_ht.describe()

        ht = ht.filter(ht.alleles[1] != "*")
        metric_ht = metric_ht.filter(metric_ht.alleles[1] != "*")
        vqsr_ht = hl.read_table(var_annotations_ht_path(data_source, freeze, "vqsr"))
        ht = ht.filter(~vqsr_ht[ht.key].filters.contains("LowQual"))
        metric_ht = metric_ht.filter(
            ~vqsr_ht[metric_ht.key].filters.contains("LowQual")
        )

        (
            metric_snvs,
            metrics_indels,
            interval_metric_snvs,
            interval_metrics_indels,
        ) = metric_ht.aggregate(
            [
                hl.agg.count_where(
                    hl.is_snp(metric_ht.alleles[0], metric_ht.alleles[1])
                ),
                hl.agg.count_where(
                    ~hl.is_snp(metric_ht.alleles[0], metric_ht.alleles[1])
                ),
                hl.agg.count_where(
                    hl.is_snp(metric_ht.alleles[0], metric_ht.alleles[1])
                    & metric_ht.interval_qc_pass
                ),
                hl.agg.count_where(
                    ~hl.is_snp(metric_ht.alleles[0], metric_ht.alleles[1])
                    & metric_ht.interval_qc_pass
                ),
            ]
        )

        snvs, indels, interval_snvs, interval_indels = ht.aggregate(
            [
                hl.agg.count_where(hl.is_snp(ht.alleles[0], ht.alleles[1])),
                hl.agg.count_where(~hl.is_snp(ht.alleles[0], ht.alleles[1])),
                hl.agg.count_where(
                    hl.is_snp(ht.alleles[0], ht.alleles[1])
                    & metric_ht[ht.key].interval_qc_pass
                ),
                hl.agg.count_where(
                    ~hl.is_snp(ht.alleles[0], ht.alleles[1])
                    & metric_ht[ht.key].interval_qc_pass
                ),
            ]
        )

        ht = ht.annotate_globals(
            global_counts=hl.struct(
                snvs=metric_snvs,
                indels=metrics_indels,
                interval_snvs=interval_metric_snvs,
                interval_indels=interval_metrics_indels,
            ),
            counts=hl.struct(
                snvs=snvs,
                indels=indels,
                interval_snvs=interval_snvs,
                interval_indels=interval_indels,
            ),
        )

        ht = ht.annotate(
            snv=hl.is_snp(ht.alleles[0], ht.alleles[1]),
            score=metric_ht[ht.key].score,
            global_rank=metric_ht[ht.key].rank,
            interval_global_rank=metric_ht[ht.key].interval_rank,
            interval_qc_pass=metric_ht[ht.key].interval_qc_pass,
            # TP => allele is found in both data sets
            n_tp=ht.concordance[3][3]
            + ht.concordance[3][4]
            + ht.concordance[4][3]
            + ht.concordance[4][4],
            # FP => allele is found only in test data set
            n_fp=hl.sum(ht.concordance[3][:2]) + hl.sum(ht.concordance[4][:2]),
            # FN => allele is found only in truth data set
            n_fn=hl.sum(ht.concordance[:2].map(lambda x: x[3] + x[4])),
        )

        ht = add_rank(
            ht, -1.0 * ht.score, subrank_expr={"interval_rank": ht.interval_qc_pass}
        )

        ht = ht.annotate(
            rank=[
                hl.tuple(
                    [
                        "global_rank",
                        (ht.global_rank + 1)
                        / hl.cond(
                            ht.snv,
                            ht.globals.global_counts.snvs,
                            ht.globals.global_counts.indels,
                        ),
                    ]
                ),
                hl.tuple(
                    [
                        "truth_sample_rank",
                        (ht.rank + 1)
                        / hl.cond(
                            ht.snv, ht.globals.counts.snvs, ht.globals.counts.indels
                        ),
                    ]
                ),
                hl.tuple(
                    [
                        "interval_global_rank",
                        (ht.interval_global_rank + 1)
                        / hl.cond(
                            ht.snv,
                            ht.globals.global_counts.interval_snvs,
                            ht.globals.global_counts.interval_indels,
                        ),
                    ]
                ),
                hl.tuple(
                    [
                        "interval_truth_sample_rank",
                        (ht.interval_rank + 1)
                        / hl.cond(
                            ht.snv,
                            ht.globals.counts.interval_snvs,
                            ht.globals.counts.interval_indels,
                        ),
                    ]
                ),
            ]
        )

        ht = ht.explode(ht.rank)
        ht = ht.annotate(rank_name=ht.rank[0], bin=hl.int(ht.rank[1] * nbins))

        ht = (
            ht.group_by("rank_name", "snv", "bin")
            .aggregate(
                # Look at site-level metrics -> tp > fp > fn -- only important for multi-sample comparisons
                tp=hl.agg.count_where(ht.n_tp > 0),
                fp=hl.agg.count_where((ht.n_tp == 0) & (ht.n_fp > 0)),
                fn=hl.agg.count_where((ht.n_tp == 0) & (ht.n_fp == 0) & (ht.n_fn > 0)),
                min_score=hl.agg.min(ht.score),
                max_score=hl.agg.max(ht.score),
                n_alleles=hl.agg.count(),
            )
            .repartition(5)
        )

        ht.write(
            binned_concordance_path(data_source, freeze, truth_sample, metric),
            overwrite=overwrite,
        )


def main(args):
    data_source = args.data_source
    freeze = args.freeze

    if args.compute_concordance:
        if args.na12878:
            logger.info(f"Computing  NA12878 concordance")
            write_truth_concordance(
                data_source, freeze, "na12878", overwrite=args.overwrite
            )

        if args.syndip:
            logger.info(f"Computing syndip concordance")
            write_truth_concordance(data_source, freeze, "syndip", args.overwrite)

    if args.bin_concordance:
        truth_samples = []
        if args.na12878:
            truth_samples.append("na12878")
        if args.syndip:
            truth_samples.append("syndip")

        metrics = (
            []
            if not args.run_hash
            else [args.run_hash]
            if isinstance(args.run_hash, str)
            else args.run_hash
        )
        if args.vqsr:
            metrics.append("vqsr" if args.vqsr_type == "AS" else "AS_TS_vqsr")

        for truth_sample in truth_samples:
            for metric in metrics:
                logger.info(
                    f"Creating binned concordance table for {truth_sample} for metric {metric}"
                )
                create_binned_concordance(
                    data_source,
                    freeze,
                    truth_sample,
                    metric,
                    args.n_bins,
                    args.overwrite,
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", help="Prints debug statements", action="store_true")
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
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

    action_grp = parser.add_argument_group("Actions")
    action_grp.add_argument(
        "--compute_concordance",
        help="Computes concordance for each of the dataset specified separately.",
        action="store_true",
    )
    action_grp.add_argument(
        "--bin_concordance",
        help="Merges individual concordance results with specified metrics ranks.",
        action="store_true",
    )

    data_grp = parser.add_argument_group("Datasets available")
    data_grp.add_argument("--na12878", help="NA12878 ", action="store_true")
    data_grp.add_argument("--syndip", help="Syndip", action="store_true")

    bin_grp = parser.add_argument_group("Bin concordance options")
    bin_grp.add_argument("--run_hash", help="RF hash(es) for annotation.", nargs="+")
    bin_grp.add_argument(
        "--vqsr", help="When set, annotates with VQSR rank file.", action="store_true"
    )
    parser.add_argument(
        "--vqsr_type",
        help="What type of VQSR was run: allele-specific, or allele-specific with transmitted singletons",
        type=str,
        choices=["AS", "AS_TS"],
        default="AS",
    )
    bin_grp.add_argument(
        "--n_bins",
        help="Number of bins for the binned file (default: 100)",
        default=100,
        type=int,
    )

    args = parser.parse_args()
    if args.compute_concordance and not (
        args.na12878 or args.syndip or args.omes or args.omes_by_platform
    ):
        sys.exit(
            "Error At least one dataset should be specified when running --compute_concordance."
        )

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
