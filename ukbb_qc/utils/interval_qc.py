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
from ukbb_qc.utils.sparse_utils import compute_callrate_dp_mt


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("interval_qc")
logger.setLevel(logging.INFO)


def interval_qc(
    target_mt: hl.MatrixTable,
    split_by_sex: bool = False,
    ht_n_partitions: int = 100
) -> hl.Table:
    """
    Determines the percent of samples reaching 10x, 15x, 20x, 25x, and 30x mean coverage in each interval

    :param MatrixTable target_mt: Input MatrixTable for interval QC, output by compute_callrate_dp_mt. Annotated with interval. 
    :param bool split_by_sex: Whether the interval QC should be stratified by sex. if True, mt must be annotated with sex_karyotype.
    :param int ht_n_partitions: Number of desired partitions for output Table
    :return: Table with mean DP and percent samples with coverage at 10x, 15x, 20x, 25x, 30x across target.
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

    # NOTE: this function will densify
    if args.compute_callrate_mt:
        logger.info("Reading in raw MT...")
        mt = get_ukbb_data(
            data_source, freeze, raw=True, adj=False, key_by_locus_and_alleles=True
        )
        logger.info(
            f"Total number of variants in raw unsplit matrix table: {mt.count_rows()}"
        )
        capture_ht = hl.read_table(capture_ht_path(data_source))
        compute_callrate_dp_mt(data_source, freeze, mt, capture_ht, autosomes_only=False)

    logger.info("Reading in call rate MT...")
    mt = hl.read_matrix_table(callrate_mt_path(data_source, freeze, interval_filtered=False))
    if args.autosomes:
        logger.info("Filtering to autosomes...")
        mt = filter_to_autosomes(mt)

        logger.info("Starting interval QC...")
        target_ht = interval_qc(mt)
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
            split_by_sex=True
        )
        target_ht.write(
            interval_qc_path(data_source, freeze, "sex_chr"), overwrite=args.overwrite,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
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
        "--n_partitions", help="Desired number of partitions for output", type=int
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
