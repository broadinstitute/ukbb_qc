import argparse
import hail as hl
import logging
from gnomad_hail.utils.generic import filter_to_autosomes
from gnomad_hail.utils.slack import try_slack
from ukbb_qc.resources.basics import (
    CURRENT_FREEZE,
    capture_ht_path,
    get_checkpoint_path,
    get_ukbb_data,
)
from ukbb_qc.resources.sample_qc import interval_qc_path, sex_ht_path
from typing import Optional


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("interval_qc")
logger.setLevel(logging.INFO)


def interval_qc(
    mt: hl.MatrixTable,
    interval_ht: hl.Table,
    split_by_sex: bool = False,
    mt_n_partitions: int = 5000,
    ht_n_partitions: int = 100,
    checkpoint_path: Optional[str] = None,
) -> hl.Table:
    """
    Determines the percent of samples reaching 10x, 15x, 20x, 25x, and 30x mean coverage in each interval

    :param MatrixTable mt: Matrix table to use for interval QC
    :param Table interval_ht: Table containing the intervals that interval QC is performed over
    :param bool split_by_sex: Whether the interval QC should be stratified by sex. if True, mt must be annotated with sex_karyotype.
    :param int mt_n_partitions: Number of desired partitions for intermediate MatrixTable
    :param int ht_n_partitions: Number of desired partitions for output Table
    :param str checkpoint_path: Optional path to a file to checkpoint the mean_dp per sample across each interval MT
    :return: Table with percent samples with coverage at 10x, 15x, 20x, 25x, and 30x coverage
    :rtype: Table
    """

    interval_ht = interval_ht.annotate(
        interval_label=interval_ht.interval,
        interval_size=(interval_ht.interval.end - interval_ht.interval.start),
    )

    mt = mt.annotate_rows(
        interval=interval_ht[mt.locus].interval_label,
        interval_size=interval_ht[mt.locus].interval_size,
    )

    target_mt = (
        mt.group_rows_by(mt.interval)
        .aggregate_entries(
            mean_dp=hl.agg.mean(
                hl.if_else(
                    mt.LGT.is_hom_ref(),
                    mt.DP * (mt.END - mt.locus.position),
                    hl.if_else(hl.is_defined(mt.DP), mt.DP, 0),
                )
            ),
            pct_gt_20x=hl.agg.fraction(hl.cond(hl.is_defined(mt.DP), mt.DP, 0) >= 20),
            pct_dp_defined=hl.agg.count_where(hl.is_defined(mt.DP)) / hl.agg.count(),
        )
        .result()
    )
    target_mt = target_mt.naive_coalesce(mt_n_partitions)
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
    target_ht = target_ht.naive_coalesce(ht_n_partitions)
    target_ht.describe()
    return target_ht


def main(args):
    hl.init(log="/interval_qc.log")

    if not args.autosomes and not args.sex_chr:
        logger.warning("Must choose one of autosomes or sex_chr options")

    data_source = args.data_source
    freeze = args.freeze
    ht = hl.read_table(capture_ht_path(data_source))

    checkpoint_path = get_checkpoint_path(
        data_source,
        freeze,
        name=f"DP_only_checkpoint_interval_qc_broad_freeze{freeze}",
        mt=True,
    )
    if not hl.utils.hadoop_exists(f"{checkpoint_path}/_SUCCESS"):
        mt = get_ukbb_data(data_source, freeze, raw=True, adj=False, split=False)
        mt = mt.key_rows_by("locus", "alleles")
        mt = mt.select_entries("DP", "END", "LGT")
        mt.write(checkpoint_path, overwrite=True)

    if args.autosomes:
        mt = hl.read_matrix_table(checkpoint_path)
        mt = filter_to_autosomes(mt)
        target_ht = interval_qc(
            mt,
            ht,
            checkpoint_path=get_checkpoint_path(
                data_source, freeze, "coverage_by_target_samples", mt=True
            ),
        )
        target_ht.write(
            interval_qc_path(data_source, freeze, "autosomes"),
            overwrite=args.overwrite,
        )

    if args.sex_chr:
        mt = hl.read_matrix_table(checkpoint_path)
        mt = hl.filter_intervals(
            mt,
            [
                hl.parse_locus_interval("chrX", reference_genome="GRCh38"),
                hl.parse_locus_interval("chrY", reference_genome="GRCh38"),
            ],
        )
        sex_ht = hl.read_table(sex_ht_path(data_source, freeze)).select("sex_karyotype")
        mt = mt.annotate_cols(**sex_ht[mt.col_key])
        mt.describe()
        mt = mt.filter_cols((mt.sex_karyotype == "XX") | (mt.sex_karyotype == "XY"))
        target_ht = interval_qc(
            mt,
            ht,
            split_by_sex=True,
            checkpoint_path=get_checkpoint_path(
                data_source, freeze, "coverage_by_target_samples", mt=True
            ),
        )
        target_ht.write(
            interval_qc_path(data_source, freeze, "sex_chr"), overwrite=args.overwrite,
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
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int,
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
