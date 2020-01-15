import hail as hl
import ukbb_qc.resources as urc
from gnomad_hail import try_slack, filter_to_autosomes
import argparse
import logging

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("interval_qc")
logger.setLevel(logging.INFO)


def main(args):
    hl.init(log="/interval_qc.log")

    data_source = args.data_source
    freeze = args.freeze

    mt = urc.get_ukbb_data(data_source, freeze, raw=True, adj=False, split=False)
    mt = mt.key_rows_by("locus", "alleles")
    mt = mt.select_entries(mt.DP)

    if args.chromosome is not None:
        mt = hl.filter_intervals(
            mt,
            [
                hl.parse_locus_interval(
                    f"chr{args.chromosome}", reference_genome="GRCh38"
                )
            ],
        )
    else:
        mt = filter_to_autosomes(mt)

    ht = hl.read_table(urc.capture_ht_path(data_source, freeze))
    ht = ht.annotate(interval_str=hl.str(ht.interval))
    mt = mt.annotate_rows(**ht[mt.locus])
    ht = ht.key_by("interval_str")
    mt = mt.checkpoint(
        urc.get_mt_checkpoint_path(
            data_source, freeze, "DP_only_checkpoint_interval_qc"
        ),
        overwrite=args.overwrite,
    )

    target_mt = (
        mt.group_rows_by("interval_str")
        .partition_hint(5000)
        .aggregate_entries(
            mean_dp=hl.agg.filter(
                (hl.is_defined(mt.DP)) & (~hl.is_nan(mt.DP)), hl.agg.mean(mt.DP)
            ),
            pct_gt_20x=hl.agg.filter(
                hl.is_defined(mt.DP), hl.agg.fraction(mt.DP >= 20)
            ),
            pct_dp_defined=hl.agg.count_where(hl.is_defined(mt.DP)) / hl.agg.count(),
        )
        .result()
    )
    target_mt = target_mt.naive_coalesce(5000)
    target_mt = target_mt.checkpoint(
        urc.get_mt_checkpoint_path(data_source, freeze, "coverage_by_target_samples"),
        overwrite=args.overwrite,
    )
    target_mt = target_mt.annotate_rows(
        target_num_mean_dp_defined=hl.agg.count_where(~hl.is_nan(target_mt.mean_dp)),
        target_mean_dp=hl.agg.filter(
            ~hl.is_nan(target_mt.mean_dp), hl.agg.mean(target_mt.mean_dp)
        ),
        target_pct_gt_20x=hl.agg.mean(target_mt.pct_gt_20x),
        target_mean_pct_dp_defined=hl.agg.mean(target_mt.pct_dp_defined),
        pct_samples_15x=hl.agg.fraction(target_mt.mean_dp >= 15),
        pct_samples_20x=hl.agg.fraction(target_mt.mean_dp >= 20),
        pct_samples_25x=hl.agg.fraction(target_mt.mean_dp >= 25),
        pct_samples_30x=hl.agg.fraction(target_mt.mean_dp >= 30),
    )

    target_ht = target_mt.rows()
    target_ht = target_ht.annotate(**ht[target_ht.interval_str])
    target_ht = target_ht.key_by("interval")
    target_ht = target_ht.naive_coalesce(100)

    target_ht.write(
        urc.interval_qc_path(data_source, freeze, chrom=args.chromosome),
        overwrite=args.overwrite,
    )
    target_ht = hl.read_table(urc.interval_qc_path(data_source, freeze))

    target_ht.describe()


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
        "-s",
        "--data_source",
        help="Data source",
        choices=["regeneron", "broad"],
        default="broad",
    )
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=urc.CURRENT_FREEZE, type=int
    )

    parser.add_argument(
        "-c",
        "--chromosome",
        type=int,
        choices=range(1, 23),
        help="If set it will only run the given chromosome",
    )

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
