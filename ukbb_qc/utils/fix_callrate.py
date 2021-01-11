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

    hl.init(log="/fix_callrate.log", default_reference="GRCh38")

    data_source = "broad"
    freeze = CURRENT_FREEZE

    try:
        # Read in callrate MT
        # Rows in callrate MT (from .describe()): 192977
        callrate_mt = hl.read_matrix_table(
            callrate_mt_path(data_source, freeze, interval_filtered=False)
        )

        # Read in raw MT and capture HT
        mt = get_ukbb_data(
            data_source, freeze, split=False, key_by_locus_and_alleles=False, raw=True,
        )
        capture_ht = hl.read_table(capture_ht_path(data_source))

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
        mt.describe()

        # Create HT with intervals and number of variants per interval
        ht = mt.rows()
        print(f"Check if HT count matches 192977: {ht.count == 192977}")

        # Annotate current callrate MT with n_var annotation
        callrate_mt = callrate_mt.annotate_rows(n_var=ht[callrate_mt.interval].n_var)
        callrate_mt.describe()

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


with slack_notifications(slack_token, "@kc"):
    main()
