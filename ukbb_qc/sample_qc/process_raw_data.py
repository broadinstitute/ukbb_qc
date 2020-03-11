import argparse
import hail as hl
import logging
from gnomad.utils.generic import filter_to_autosomes
from gnomad.utils.gnomad_functions import filter_to_adj
from gnomad.utils.sample_qc import get_qc_mt
from gnomad.utils.sparse_mt import densify_sites
from ukbb_qc.resources.basics import (
    get_checkpoint_path,
    get_ukbb_data,
    last_END_positions_ht_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import qc_mt_path
from ukbb_qc.resources.variant_qc import get_truth_sample_info
from ukbb_qc.utils.utils import get_qc_mt_sites


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("process_raw_data")
logger.setLevel(logging.INFO)


def main(args):
    hl.init(log="/process_raw_data.log", default_reference="GRCh38")

    data_source = "broad"
    freeze = args.freeze

    if args.compute_qc_mt:
        logger.info("Reading in raw MT...")
        mt = get_ukbb_data(
            data_source, freeze, raw=True, adj=False, key_by_locus_and_alleles=True,
        )
        mt = mt.select_entries(
            "DP", GT=mt.LGT, adj=get_adj_expr(mt.LGT, mt.GQ, mt.DP, mt.LAD)
        )
        logger.info(
            f"Total number of variants in raw unsplit matrix table: {mt.count_rows()}"
        )

        logger.info("Reading in QC MT sites from tranche 2/freeze 5...")
        if not hl.utils.hadoop_exists(f"{qc_sites_path()}/_SUCCESS"):
            get_qc_mt_sites()
        qc_sites_ht = hl.read_table(qc_sites_path())

        logger.info("Densifying sites...")
        last_END_ht = hl.read_table(last_END_positions_ht_path(data_source, freeze))
        mt = densify_sites(mt, qc_sites_ht, last_END_ht)

        logger.info("Checkpointing densified MT")
        mt = mt.checkpoint(
            get_checkpoint_path(data_source, freeze, name="dense", mt=True,),
            overwrite=True,
        )

        logger.info("Repartitioning densified MT")
        mt = mt.naive_coalesce(args.n_partitions)
        mt = mt.checkpoint(
            get_checkpoint_path(
                data_source, freeze, name="dense.repartitioned", mt=True,
            ),
            overwrite=True,
        )
        logger.info("Adding info and low QUAL annotations and filtering to adj...")
        # NOTE: Need MQ, QD, FS for hard filters
        info_expr = get_site_info_expr(mt)
        info_expr = info_expr.annotate(**get_as_info_expr(mt))
        mt = mt.annotate_rows(info=info_expr)
        mt = mt.annotate_rows(
            lowqual=get_lowqual_expr(
                mt.alleles, mt.info.QUALapprox, indel_phred_het_prior=40,
            )
        )
        mt = filter_to_adj(mt)

        logger.info("Checkpointing MT...")
        mt = mt.checkpoint(
            get_checkpoint_path(
                data_source,
                freeze,
                name=f"{data_source}.freeze_{freeze}.qc_sites.mt",
                mt=True,
            ),
            overwrite=True,
        )

        # NOTE: adding decoy and segdup hg38 resources/filters are still pending
        qc_mt = get_qc_mt(
            mt,
            min_af=args.min_af,
            min_callrate=args.min_callrate,
            apply_hard_filters=True,
            ld_r2=None,
            filter_lcr=False,
            filter_decoy=False,
            filter_segdup=False,
        )
        qc_mt = qc_mt.checkpoint(
            qc_mt_path(data_source, freeze, ld_pruned=True), overwrite=args.overwrite,
        )
        logger.info(
            f"Total number of bi-allelic, high-callrate, common SNPs for sample QC: {qc_mt.count_rows()}"
        )

    if args.compute_sample_qc_ht:
        qc_mt = hl.read_matrix_table(qc_mt_path(data_source, freeze))
        qc_mt = filter_to_autosomes(qc_mt)
        qc_ht = hl.sample_qc(qc_mt).cols().select("sample_qc")
        qc_ht = qc_ht.transmute(
            sample_qc=qc_ht.sample_qc.annotate(
                dp_mean=qc_ht.sample_qc.dp_stats.mean,
                dp_stdev=qc_ht.sample_qc.dp_stats.stdev,
            ).select("call_rate", "dp_mean", "dp_stdev")
        )
        qc_ht.write(qc_ht_path(data_source, freeze), overwrite=args.overwrite)

    if args.extract_truth_samples:
        mt = get_ukbb_data(data_source, freeze, raw=True, adj=False, split=False,)
        truth_samples = get_truth_sample_info(data_source, freeze)
        samples = [s["s"] for s in truth_samples.values()]
        mt = mt.filter_cols(hl.literal(samples).contains(mt.s))
        logger.info(f"Extracting {mt.count_cols()} truth samples.")
        mt.checkpoint(
            get_mt_checkpoint_path(data_source, freeze, "truth_samples"),
            overwrite=args.overwrite,
        )

        for truth_sample in get_truth_sample_info(data_source, freeze):
            mt = hl.read_matrix_table(
                get_mt_checkpoint_path(data_source, freeze, "truth_samples")
            )
            mt = mt.filter_cols(mt.s == truth_samples[truth_sample]["s"])
            mt = mt.key_rows_by("locus", "alleles")
            mt = hl.experimental.sparse_split_multi(mt)
            mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
            mt.naive_coalesce(args.n_partitions).write(
                truth_samples[truth_sample]["mt"], overwrite=args.overwrite,
            )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int,
    )
    parser.add_argument(
        "--compute_qc_mt",
        help="Compute matrix to be used in sample qc",
        action="store_true",
    )
    parser.add_argument(
        "--min_af",
        help="Minimum variant allele frequency to retain variant in qc matrix table.",
        default=0.001,
        type=float,
    )
    parser.add_argument(
        "--min_callrate",
        help="Minimum variant callrate to retain variant in qc matrix table.",
        default=0.99,
        type=float,
    )
    parser.add_argument(
        "--compute_sample_qc_ht",
        help="Compute sample qc on qc matrix table",
        action="store_true",
    )
    parser.add_argument(
        "--extract_truth_samples",
        help="Extract truth samples from matrix table",
        action="store_true",
    )
    parser.add_argument(
        "--n_partitions",
        help="Desired number of partitions for output Table/MatrixTable",
        default=5000,
        type=int,
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to.",
    )

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
