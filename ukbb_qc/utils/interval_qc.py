import argparse
import logging
from typing import List

import hail as hl

from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.slack import try_slack
from ukbb_qc.resources.basics import (
    capture_ht_path,
    get_ukbb_data,
    logging_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import callrate_mt_path, interval_qc_path, sex_ht_path
from ukbb_qc.utils.sparse_utils import compute_interval_callrate_dp_mt


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("interval_qc")
logger.setLevel(logging.INFO)


def interval_qc(
    target_mt: hl.MatrixTable,
    target_pct_gt_cov: List = [10, 20],
    pct_sample_cov: List = [10, 20],
    split_by_sex: bool = False,
    n_partitions: int = 100,
) -> hl.Table:
    """
    Determines the percent of samples reaching mean coverage at specified levels in each interval.
    Assumes that the field for target_pct_gt_cov exists in the input MatrixTable.

    :param MatrixTable target_mt: Input MatrixTable for interval QC, output by compute_callrate_dp_mt. Annotated with interval. 
    :param int target_pct_gt_cov: Coverage levels to check for each target. Default is 10 and 20. This field must be in target_mt!
    :param List pct_sample_cov: Desired coverage levels at which to check sample coverage. Default 10x, 15x, 20x, 25x, 30x.
    :param bool split_by_sex: Whether the interval QC should be stratified by sex. if True, mt must be annotated with sex_karyotype.
    :param int n_partitions: Number of desired partitions for output Table.
    :return: Table with mean DP and percent samples with coverage at pct_sample_cov across target.
    :rtype: Table
    """
    if split_by_sex:
        target_mt = target_mt.annotate_rows(
            target_mean_dp=hl.agg.group_by(
                target_mt.sex_karyotype,
                hl.agg.filter(
                    ~hl.is_nan(target_mt.mean_dp), hl.agg.mean(target_mt.mean_dp)
                ),
            ),
            **{
                f"target_pct_gt_{cov}x": hl.agg.group_by(
                    target_mt.sex_karyotype, hl.agg.mean(target_mt[f"pct_gt_{cov}x"]),
                )
                for cov in target_pct_gt_cov
            },
            **{
                f"pct_samples_{cov}x": hl.agg.group_by(
                    target_mt.sex_karyotype, hl.agg.fraction(target_mt.mean_dp >= cov)
                )
                for cov in pct_sample_cov
            },
        )
    else:
        target_mt = target_mt.annotate_rows(
            target_mean_dp=hl.agg.filter(
                ~hl.is_nan(target_mt.mean_dp), hl.agg.mean(target_mt.mean_dp)
            ),
            **{
                f"target_pct_gt_{cov}x": hl.agg.mean(target_mt[f"pct_gt_{cov}x"])
                for cov in target_pct_gt_cov
            },
            **{
                f"pct_samples_{cov}x": hl.agg.fraction(target_mt.mean_dp >= cov)
                for cov in pct_sample_cov
            },
        )

    target_ht = target_mt.rows()
    target_ht = target_ht.naive_coalesce(n_partitions)
    return target_ht


def main(args):

    data_source = "broad"
    freeze = args.freeze
    log = (
        "interval_qc_callrate_mt.log"
        if args.compute_interval_callrate_mt
        else "interval_qc.log"
    )
    hl.init(log=f"/{log}", default_reference="GRCh38")

    try:
        if not args.autosomes and not args.sex_chr:
            logger.warning("Must choose one of autosomes or sex_chr options")

        target_pct_gt_cov = list(map(int, args.target_cov.split(",")))
        pct_sample_cov = list(map(int, args.sample_cov.split(",")))
        n_partitions = args.n_partitions

        if (args.compute_interval_callrate_mt and args.autosomes) or (
            args.compute_interval_callrate_mt and args.sex_chr
        ):
            logger.warning(
                "n_partitions arg will be used for both callrate MT and interval QC"
                "if args.compute_interval_callrate_mt and any interval QC args are set!"
            )

        if args.compute_interval_callrate_mt:
            logger.warning("Computing the call rate MT requires a densify!\n")
            logger.info("Reading in raw MT...")
            mt = get_ukbb_data(
                data_source,
                freeze,
                split=False,
                key_by_locus_and_alleles=True,
                raw=True,
                repartition=args.repartition,
                n_partitions=n_partitions,
            )
            capture_ht = hl.read_table(capture_ht_path(data_source))
            compute_interval_callrate_dp_mt(
                data_source,
                freeze,
                mt,
                capture_ht,
                autosomes_only=False,
                target_pct_gt_cov=target_pct_gt_cov,
            )

        logger.info("Reading in call rate MT...")
        mt = hl.read_matrix_table(
            callrate_mt_path(data_source, freeze, interval_filtered=False)
        )
        mt = mt.annotate_rows(locus=mt.interval.start)
        mt = mt.key_rows_by("locus")

        if args.autosomes:
            logger.info("Filtering to autosomes...")
            mt = filter_to_autosomes(mt)
            mt = mt.key_rows_by("interval").drop("locus")

            logger.info("Starting interval QC...")
            target_ht = interval_qc(
                mt,
                target_pct_gt_cov=target_pct_gt_cov,
                pct_sample_cov=pct_sample_cov,
                split_by_sex=False,
                n_partitions=n_partitions,
            )
            target_ht.write(
                interval_qc_path(data_source, freeze, "autosomes"),
                overwrite=args.overwrite,
            )

        if args.sex_chr:
            logger.info("Filtering to sex chromosomes...")
            mt = hl.filter_intervals(
                mt,
                [
                    hl.parse_locus_interval("chrX", reference_genome="GRCh38"),
                    hl.parse_locus_interval("chrY", reference_genome="GRCh38"),
                ],
            )
            mt = mt.key_rows_by("interval").drop("locus")

            logger.info("Filtering to XX and XY samples...")
            sex_ht = hl.read_table(sex_ht_path(data_source, freeze)).select(
                "sex_karyotype"
            )
            mt = mt.annotate_cols(**sex_ht[mt.col_key])
            mt = mt.filter_cols((mt.sex_karyotype == "XX") | (mt.sex_karyotype == "XY"))

            logger.info("Starting interval QC...")
            target_ht = interval_qc(
                mt,
                target_pct_gt_cov=target_pct_gt_cov,
                pct_sample_cov=pct_sample_cov,
                split_by_sex=True,
                n_partitions=n_partitions,
            )
            target_ht.write(
                interval_qc_path(data_source, freeze, "sex_chr"),
                overwrite=args.overwrite,
            )

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(logging_path(data_source, freeze))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int,
    )
    parser.add_argument(
        "--n_partitions", help="Desired number of partitions for output", type=int,
    )
    parser.add_argument(
        "--repartition",
        help="Repartition raw MT on read. Needs to be true for tranche 3/freeze 6/300K.",
        action="store_true",
    )
    qc_type = parser.add_mutually_exclusive_group()
    qc_type.add_argument(
        "--autosomes",
        action="store_true",
        help="If set it will only run the autosomes",
    )
    qc_type.add_argument(
        "--sex_chr",
        action="store_true",
        help="If set it will only run the sex chromosomes",
    )
    parser.add_argument(
        "--compute_interval_callrate_mt",
        help="Computes an interval by sample mt of callrate and depth",
        action="store_true",
    )
    parser.add_argument(
        "--target_cov",
        help="Coverage levels to check per target. Default is 10x and 20x (10x to account for sex chromosomes)",
        default="10,20",
    )
    parser.add_argument(
        "--sample_cov",
        help="Comma separated list of desired coverage levels at which to check sample coverage. Default is 10x and 20x",
        default="10,20",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
