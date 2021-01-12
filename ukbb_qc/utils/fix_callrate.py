import argparse
import hail as hl

from gnomad.utils.annotations import bi_allelic_expr
from gnomad.utils.slack import slack_notifications
from ukbb_qc.resources.basics import (
    capture_ht_path,
    get_checkpoint_path,
    get_ukbb_data,
    logging_path,
)
from ukbb_qc.resources.sample_qc import callrate_mt_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.slack_creds import slack_token


def main():
    """
    We noticed after running the first round of hard filters on freeze 7/the 450k MT that every sample had a
    callrate of 100%. Upon closer inspection of the callrate MT, we found that the `total` field for every sample
    matched their `n_defined` field, meaning that the `total` field had incorrect values.

    Both fields are calculated using:
    ```
    mt = mt.group_rows_by(mt.interval).aggregate(
        n_defined=hl.agg.count_where(hl.is_defined(mt.GT)),
        total=hl.agg.count(),
        dp_sum=hl.agg.sum(mt.DP),
        mean_dp=hl.agg.mean(mt.DP),
        **{
            f"pct_gt_{cov}x": hl.agg.fraction(mt.DP >= cov) for cov in target_pct_gt_cov
        },
        pct_dp_defined=hl.agg.count_where(mt.DP > 0) / hl.agg.count(),
    )
    ```
    The reason why `total` is not calculated correctly is because the entry aggregation + `count()` 
    does not count filtered entries. Not every genotype is populated during a densify, which means the 
    unpopulated entries would appear to be "filtered", so this `count()` is undercounting genotypes.
    For example, if a sample has a ref block from chr1:1-5 and a variant at chr1:7, and the MT has rows 
    chr:1 and then chr:7, then after the densify, genotypes for that sample from chr1:2-5 don’t actually get populated 
    (even though we know they are defined).

    As an important note, the reason the callrate MT code worked on freeze 6/the 300k MT is because of a hail bug.
    The bug was that previous entry aggregations counted filtered entries. We ran the 300k callrate MT code in April 2020,
    and the bug was fixed three months later in June.

    The fix for `total` for the 450k MT and moving forwards is to group the callrate MT by interval and 
    aggregate the *rows* to get the number of variants per interval.
    """

    hl.init(log="/fix_callrate.log", default_reference="GRCh38")

    data_source = "broad"
    freeze = CURRENT_FREEZE

    try:
        # Read in callrate MT
        # Rows in 450K callrate MT (from .describe()): 192977
        callrate_mt = hl.read_matrix_table(
            callrate_mt_path(data_source, freeze, interval_filtered=False)
        )

        # Read in raw MT and capture HT
        mt = get_ukbb_data(
            data_source, freeze, split=False, key_by_locus_and_alleles=False, raw=True,
        )
        capture_ht = hl.read_table(capture_ht_path(data_source))

        # Lines 70-84 are copied from the original callrate MT code
        # We are running them in this fix to make sure we are producing the same grouped MT
        # Remove ref block rows
        mt = mt.filter_rows(hl.len(mt.alleles) > 1)

        # Annotate raw MT with capture intervals and filter to rows that are within capture regions
        capture_ht = capture_ht.annotate(interval_label=capture_ht.interval)
        mt = mt.annotate_rows(interval=capture_ht.index(mt.locus).interval_label)
        mt = mt.filter_rows(hl.is_defined(mt.interval))

        # Filter to bi-allelic variants
        mt = mt.filter_rows(bi_allelic_expr(mt))

        # Key rows by interval, preventing shuffle
        mt = mt.key_rows_by(mt.interval)
        mt = mt._filter_partitions(
            [], keep=False
        )  # hack to prevent optimization back to the original execution

        # Group MT by interval and get count of variants per interval -- necessary for callrate
        mt = mt.group_rows_by(mt.interval).aggregate_rows(n_var=hl.agg.count()).result()

        # Create HT with intervals and number of variants per interval
        # This number is specifically for the 450K MT
        ht = mt.rows()
        print(f"Check if HT count matches 192977: {ht.count == 192977}")

        # Annotate current callrate MT with n_var annotation
        callrate_mt = callrate_mt.annotate_rows(n_var=ht[callrate_mt.interval].n_var)

        # Write current callrate MT to temp and overwrite MT at callrate mt path
        callrate_mt = callrate_mt.checkpoint(
            get_checkpoint_path(data_source, freeze, name="temp_callrate", mt=True)
        )
        callrate_mt = callrate_mt.write(
            callrate_mt_path(data_source, freeze), overwrite=True
        )

    finally:
        print("Copying log...")
        hl.copy_log(logging_path(data_source, freeze))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int,
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to.",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
