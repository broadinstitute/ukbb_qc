import argparse
import logging
from typing import List

import hail as hl

from gnomad.utils.slack import slack_notifications

from ukbb_qc.resources.basics import (
    dup_mt_path,
    logging_path,
    raw_mt_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.slack_creds import slack_token


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("prep_dup_map")
logger.setLevel(logging.INFO)


def main(args):
    hl.init(log="/prep_dup_map.log", default_reference="GRCh38")

    data_source = "broad"
    freeze = args.freeze
    temp_path = f"gs://broad-ukbb/broad.freeze_{freeze}/temp"

    try:
        logger.info("Reading in raw MT and filtering to chr20...")
        path = raw_mt_path(data_source, freeze)
        mt = hl.read_matrix_table(path)
        mt = hl.filter_intervals(
            mt, [hl.parse_locus_interval("chr20", reference_genome="GRCh38")]
        )
        intervals = mt._calculate_new_partitions(args.n_partitions)
        mt = hl.read_matrix_table(path, _intervals=intervals)
        ht = mt.cols()
        logger.info(f"Raw MT has {ht.count()} samples...")

        logger.info("Adding column index annotation and splitting multiallelics...")
        mt = mt.add_col_index()
        mt = hl.experimental.sparse_split_multi(mt)
        mt = mt.select_entries("GT", "END", "DP")
        mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref() | hl.is_defined(mt.END)))

        logger.info("Extracting duplicate sample IDs...")
        samples_dict = ht.aggregate(hl.agg.counter(ht.s))
        dups = []
        for s, n in samples_dict.items():
            if n > 1:
                logger.info(f"{s}: found {n} times")
                dups.append(s)

        logger.info(
            "Filtering to duplicate samples + their variants and checkpointing..."
        )
        mt = mt.filter_cols(hl.literal(dups).contains(mt.s))
        mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
        mt = mt.checkpoint(f"{temp_path}/dup_samples.mt", overwrite=True)

        logger.info("Running sample QC...")
        mt = hl.sample_qc(mt)
        mt = mt.annotate_cols(
            call_rate=mt.sample_qc.call_rate,
            n_called=mt.sample_qc.n_called,
            mean_dp=mt.sample_qc.dp_stats.mean,
        )
        ht = mt.cols().drop("sample_qc")

        logger.info("Exporting sample QC results to TSV...")
        ht.export(f"{temp_path}/dup_samples_call_rate.tsv")

        logger.info("Parsing sample QC results...")
        dups_dict = {}
        for s in set(dups):
            dups_dict[s] = {}
        with hl.hadoop_open(f"{temp_path}/dup_samples_call_rate.tsv", "r") as i:
            _ = i.readline()
            for line in i:
                sample, idx, callrate, ncalled, meandp = line.strip().split("\t")
                dups_dict[sample][idx] = (callrate, ncalled, meandp)

        indices = []
        # Check both versions of each dup sample
        for s in dups_dict:

            temp_idx = []
            temp_dp = []

            # Check each column index for each sample
            # Store depth associated with each index
            for i in dups_dict[s]:
                temp_idx.append(i)
                temp_dp.append(dups_dict[s][i][2])

            # Store first index if its depth is greater
            if temp_dp[0] > temp_dp[1]:
                indices.append(int(temp_idx[0]))
            # Otherwise, store second index if that depth is greater
            if temp_dp[1] > temp_dp[0]:
                indices.append(int(temp_idx[1]))

        logger.info("Filtering MT to version of samples with higher depths...")

        def _get_cond_expr(
            col_idx_expr: hl.expr.Int64Expression, indices: List[int]
        ) -> hl.expr.BooleanExpression:
            """
            Returns criteria to filter MT columsn based on input column index expression and list of column indices

            :param hl.expr.Int64Expression col_idx_expr: Expression containing column index.
            :param List[str] indices: List of column indices to keep.
            :return: Expression used to filter MT to desired columns.
            :rtype: hl.expr.BooleanExpression
            """
            cond_expr = col_idx_expr == int(indices[0])
            for i in range(1, len(indices)):
                cond_expr |= col_idx_expr == int(indices[i])
            return cond_expr

        higher_depth_mt = mt.filter_cols(_get_cond_expr(mt.col_idx, indices))
        higher_depth_mt = higher_depth_mt.filter_rows(
            hl.agg.any(higher_depth_mt.GT.is_non_ref())
        )

        logger.info("Reading in MT containing newest version of duplicate samples...")
        new_path = dup_mt_path(freeze)
        new_mt = hl.read_matrix_table(new_path)
        new_mt = hl.filter_intervals(
            new_mt, [hl.parse_locus_interval("chr20", reference_genome="GRCh38")]
        )
        new_intervals = new_mt._calculate_new_partitions(args.n_partitions)
        new_mt = hl.read_matrix_table(new_path, _intervals=new_intervals)
        new_mt = hl.experimental.sparse_split_multi(new_mt)
        new_mt = new_mt.filter_rows(hl.agg.any(new_mt.GT.is_non_ref()))

        logger.info(
            "Testing how many rows are in the MT with the most recent version vs filtered MT..."
        )
        test = new_mt.anti_join_rows(higher_depth_mt.rows())
        logger.info(f"Found {test.count()} rows unique to the new MT...")

        # NOTE: When I ran an anti_join using the new_mt against the checkpointed mt, I found 10 unique rows
        # I'm not sure why there were 10 unique rows, but they were all rows with star alleles

        # NOTE: I ran the above sample filtering + anti_join until there weren't any rows unique to the new MT
        # If there were unique rows, I checked for which sample's variants were unique. I then created another
        # version of the filtered MT including the version of the sample with the other column index
        # before re-running the filter_rows and anti_join
        # I wrote the final column indices and the corresponding sample IDs to the `dup_map_path` in basics.py
        # NB: gs://broad-ukbb/broad.freeze_7/notebooks/resolve_dups.ipynb

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(data_source, freeze))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int,
    )
    parser.add_argument(
        "--n_partitions",
        help="Number of desired partitions for input MT. Default is 500",
        type=int,
        default=500,
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
