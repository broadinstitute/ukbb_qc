import argparse
import hail as hl
import logging
from gnomad.utils.generic import filter_to_autosomes
from gnomad.utils.slack import try_slack
from ukbb_qc.resources.basics import (
    capture_ht_path,
    get_checkpoint_path,
    get_ukbb_data,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import callrate_mt_path, interval_qc_path, sex_ht_path
from ukbb_qc.utils.sparse_utils import compute_callrate_dp_mt


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("interval_qc")
logger.setLevel(logging.INFO)


def interval_qc(
    target_mt: hl.MatrixTable,
    target_coverage: List = [10, 20],
    coverage_levels: List = [10, 15, 20, 25, 30],
    split_by_sex: bool = False,
    n_partitions: int = 100,
) -> hl.Table:
    """
    Determines the percent of samples reaching mean coverage at specified levels in each interval.
    Assumes that the field for target_coverage exists in the input MatrixTable,

    :param MatrixTable target_mt: Input MatrixTable for interval QC, output by compute_callrate_dp_mt. Annotated with interval. 
    :param int target_coverage: Coverage levels to check for each target. Default is 10 and 20. This field must be in target_mt!
    :param List coverage_levels: Desired coverage levels at which to check sample coverage. Default 10x, 15x, 20x, 25x, 30x.
    :param bool split_by_sex: Whether the interval QC should be stratified by sex. if True, mt must be annotated with sex_karyotype.
    :param int n_partitions: Number of desired partitions for output Table.
    :return: Table with mean DP and percent samples with coverage at coverage_levels across target.
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
                for cov in target_coverage
            },
            **{
                f"pct_samples_{cov}x": hl.agg.group_by(
                    target_mt.sex_karyotype, hl.agg.fraction(target_mt.mean_dp >= cov)
                )
                for cov in coverage_levels
            },
        )
    else:
        target_mt = target_mt.annotate_rows(
            target_mean_dp=hl.agg.filter(
                ~hl.is_nan(target_mt.mean_dp), hl.agg.mean(target_mt.mean_dp)
            ),
            **{
                f"target_pct_gt_{cov}x": hl.agg.mean(target_mt[f"pct_gt_{cov}x"])
                for cov in target_coverage
            },
            **{
                f"pct_samples_{cov}x": hl.agg.fraction(target_mt.mean_dp >= cov)
                for cov in coverage_levels
            },
        )

    target_ht = target_mt.rows()
    target_ht = target_ht.naive_coalesce(n_partitions)
    return target_ht


def main(args):
    hl.init(log="/interval_qc.log")

    if not args.autosomes and not args.sex_chr:
        logger.warning("Must choose one of autosomes or sex_chr options")

    data_source = "broad"
    freeze = args.freeze
    target_coverage = list(map(int, args.target_cov.split(",")))
    coverage_levels = list(map(int, args.cov_levels.split(",")))
    n_partitions = args.n_partitions

    if args.compute_callrate_mt:
        logger.warning(
            "Computing the call rate MT requires a densify!\n"
            "Make sure you are using an autoscaling policy."
        )
        logger.info("Reading in raw MT...")
        mt = get_ukbb_data(
            data_source, freeze, raw=True, split=False, key_by_locus_and_alleles=True,
        )
        logger.info(
            f"Total number of variants in raw unsplit matrix table: {mt.count_rows()}"
        )
        capture_ht = hl.read_table(capture_ht_path(data_source))
        compute_callrate_dp_mt(
            data_source,
            freeze,
            mt,
            capture_ht,
            autosomes_only=False,
            target_coverage=target_coverage,
        )

    logger.info("Reading in call rate MT...")
    mt = hl.read_matrix_table(
        callrate_mt_path(data_source, freeze, interval_filtered=False)
    )
    if args.autosomes:
        logger.info("Filtering to autosomes...")
        mt = filter_to_autosomes(mt)

        logger.info("Starting interval QC...")
        target_ht = interval_qc(
            mt,
            target_coverage=target_coverage,
            coverage_levels=coverage_levels,
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

        logger.info("Filtering to XX and XY samples...")
        sex_ht = hl.read_table(sex_ht_path(data_source, freeze)).select("sex_karyotype")
        mt = mt.annotate_cols(**sex_ht[mt.col_key])
        mt = mt.filter_cols((mt.sex_karyotype == "XX") | (mt.sex_karyotype == "XY"))

        logger.info("Starting interval QC...")
        target_ht = interval_qc(
            mt,
            target_coverage=target_coverage,
            coverage_levels=coverage_levels,
            split_by_sex=True,
            n_partitions=n_partitions,
        )
        target_ht.write(
            interval_qc_path(data_source, freeze, "sex_chr"), overwrite=args.overwrite,
        )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int,
    )
    parser.add_argument(
        "--n_partitions", help="Desired number of partitions for output", type=int,
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
    parser.add_argument(
        "--compute_callrate_mt",
        help="Computes an interval by sample mt of callrate and depth",
        action="store_true",
    )
    parser.add_argument(
        "--target_cov",
        help="Coverage levels to check per target. Default is 10x and 20x (10x to account for sex chromosomes)",
        default="10,20",
    )
    parser.add_argument(
        "--cov_levels",
        help="Comma separated list of desired coverage levels at which to check sample coverage. Default is 10, 15, 20, 25, 30",
        default="10,15,20,25,30",
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
