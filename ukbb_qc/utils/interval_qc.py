import hail as hl
import ukbb_qc.resources as urc
from gnomad_hail import try_slack, filter_to_autosomes
import argparse
import logging
from typing import *

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("interval_qc")
logger.setLevel(logging.INFO)


def interval_qc(
    mt: hl.MatrixTable,
    interval_ht: hl.Table,
    split_by_sex: bool = False,
    checkpoint_path: Optional[str] = None,
):
    """
    Determines the percent of samples reaching 10x, 15x, 20x, 25x, and 30x mean coverage in each interval

    :param MatrixTable mt: Matrix table to use for interval QC
    :param Table interval_ht: Table containing the intervals that interval QC is performed over
    :param bool split_by_sex: Should the interval QC be stratified by sex, if True, mt must be annotated with sex_karyotype
    :param str checkpoint_path: Optional path to a file to checkpoint the mean_dp per sample across each interval MT
    :return: Table with percent samples with coverage at 10x, 15x, 20x, 25x, and 30x coverage
    :rtype: Table
    """

    interval_ht = interval_ht.annotate(
        interval_label=interval_ht.interval
    )

    mt = mt.annotate_rows(interval=interval_ht[mt.locus].interval_label)

    target_mt = (
        mt.group_rows_by(mt.interval)
        .aggregate_entries(
            mean_dp=hl.agg.mean(hl.cond(hl.is_defined(mt.DP), mt.DP, 0)),
            pct_gt_20x=hl.agg.fraction(hl.cond(hl.is_defined(mt.DP), mt.DP, 0) >= 20),
            pct_dp_defined=hl.agg.count_where(hl.is_defined(mt.DP)) / hl.agg.count(),
        )
        .result()
    )
    target_mt = target_mt.naive_coalesce(5000)
    target_mt = target_mt.checkpoint(checkpoint_path, overwrite=True)
    if split_by_sex:
        target_mt = target_mt.annotate_rows(
            target_mean_dp=hl.agg.group_by(
                target_mt.sex_karyotype,
                hl.agg.filter(
                    ~hl.is_nan(target_mt.mean_dp), hl.agg.mean(target_mt.mean_dp)
                ),
            ),
            target_pct_gt_20x=hl.agg.group_by(
                target_mt.sex_karyotype, hl.agg.mean(target_mt.pct_gt_20x)
            ),
            pct_samples_10x=hl.agg.group_by(
                target_mt.sex_karyotype, hl.agg.fraction(target_mt.mean_dp >= 10)
            ),
            pct_samples_15x=hl.agg.group_by(
                target_mt.sex_karyotype, hl.agg.fraction(target_mt.mean_dp >= 15)
            ),
            pct_samples_20x=hl.agg.group_by(
                target_mt.sex_karyotype, hl.agg.fraction(target_mt.mean_dp >= 20)
            ),
            pct_samples_25x=hl.agg.group_by(
                target_mt.sex_karyotype, hl.agg.fraction(target_mt.mean_dp >= 25)
            ),
            pct_samples_30x=hl.agg.group_by(
                target_mt.sex_karyotype, hl.agg.fraction(target_mt.mean_dp >= 30)
            ),
        )
    else:
        target_mt = target_mt.annotate_rows(
            target_mean_dp=hl.agg.filter(
                ~hl.is_nan(target_mt.mean_dp), hl.agg.mean(target_mt.mean_dp)
            ),
            target_pct_gt_20x=hl.agg.mean(target_mt.pct_gt_20x),
            pct_samples_10x=hl.agg.fraction(target_mt.mean_dp >= 10),
            pct_samples_15x=hl.agg.fraction(target_mt.mean_dp >= 15),
            pct_samples_20x=hl.agg.fraction(target_mt.mean_dp >= 20),
            pct_samples_25x=hl.agg.fraction(target_mt.mean_dp >= 25),
            pct_samples_30x=hl.agg.fraction(target_mt.mean_dp >= 30),
        )

    target_ht = target_mt.rows()
    target_ht = target_ht.naive_coalesce(100)
    target_ht.describe()

    return target_ht


def main(args):
    hl.init(log="/interval_qc.log")

    if not args.autosomes and not args.sex_chr:
        logger.warning("Must choose one of autosomes or sex_chr options")

    data_source = args.data_source
    freeze = args.freeze
    ht = hl.read_table(urc.capture_ht_path(data_source, freeze))

    checkpoint_path = urc.get_mt_checkpoint_path(
        data_source, freeze, "DP_only_checkpoint_interval_qc_broad_tranche1" #"DP_only_checkpoint_interval_qc"
    )
    if not hl.utils.hadoop_exists(f"{checkpoint_path}/_SUCCESS"):
        mt = urc.get_ukbb_data(data_source, freeze, raw=True, adj=False, split=False)
        mt = mt.key_rows_by("locus", "alleles")
        mt = mt.select_entries(mt.DP)
        mt.write(checkpoint_path, overwrite=True)

    if args.autosomes:
        mt = hl.read_matrix_table(checkpoint_path)
        mt = filter_to_autosomes(mt)
        target_ht = interval_qc(
            mt,
            ht,
            checkpoint_path=urc.get_mt_checkpoint_path(
                data_source, freeze, "coverage_by_target_samples"
            ),
        )
        target_ht.write(
            urc.interval_qc_path(data_source, freeze, "autosomes"),
            overwrite=args.overwrite,
        )
    if args.sex_chr:
        mt = hl.read_matrix_table(checkpoint_path)
        mt = hl.filter_intervals(
            mt,
            [
                hl.parse_locus_interval(f"chrX", reference_genome="GRCh38"),
                hl.parse_locus_interval(f"chrY", reference_genome="GRCh38")
            ]
        )
        sex_ht = hl.read_table(urc.sex_ht_path(data_source, freeze))
        mt = mt.annotate_cols(**sex_ht[mt.col_key])
        mt.describe()
        mt = mt.filter_cols((mt.sex_karyotype == "XX") | (mt.sex_karyotype == "XY"))
        target_ht = interval_qc(
            mt,
            ht,
            split_by_sex=True,
            checkpoint_path=urc.get_mt_checkpoint_path(
                data_source, freeze, "coverage_by_target_samples"
            ),
        )
        target_ht.write(
            urc.interval_qc_path(data_source, freeze, "sex_chr"),
            overwrite=args.overwrite,
        )


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
        "-f",
        "--freeze",
        help="Data freeze to use",
        default=urc.CURRENT_FREEZE,
        type=int,
    )

    parser.add_argument(
        "--autosomes",
        action="store_true",
        help="If set it will only run the autosomes",
    )

    parser.add_argument(
        "--sex_chr",
        action="store_true",
        help="If set it will only run the sex chromosomes",
    )

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)