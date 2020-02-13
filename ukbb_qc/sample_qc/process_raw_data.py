import argparse
import hail as hl
import logging
from gnomad_hail.utils.generic import filter_to_autosomes
from gnomad_hail.utils.gnomad_functions import filter_to_adj
from gnomad_hail.utils.sample_qc import get_qc_mt
from gnomad_hail.utils.sparse_mt import compute_last_ref_block_end, densify_sites
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE, get_checkpoint_path, get_ukbb_data
from ukbb_qc.resources.sample_qc import array_sample_map_ht_path, qc_mt_path 
from ukbb_qc.utils.utils import get_qc_mt_sites


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("process_raw_data")
logger.setLevel(logging.INFO)


def main(args):
    hl.init(log="/process_raw_data.log", default_reference="GRCh38")

    data_source = args.data_source
    freeze = args.freeze

    if args.create_exome_array_id_map_ht:
        logger.info("Checking for sample discrepancies between MatrixTable and linking file...")
        sample_map_ht = import_array_exome_id_map_ht(data_source, freeze)
        exome_ht = get_ukbb_data(data_source, freeze, raw=True, split=False, key_by_locus_and_alleles=True).cols()
        exome_ht = exome_ht.annotate(
            ukbb_app_26041_id=sample_map_ht[exome_ht.s.split("_")[1]].eid_26041,
            **sample_map_ht[exome_ht.s.split("_")[1]],
        )
        logger.info(f"Total number of samples in the exome data: {exome_ht.count()}...")

        # Check for samples that are in exome file but not in array file or vice versa
        s_exome_not_in_map = exome_ht.filter(
            hl.is_missing(exome_ht.ukbb_app_26041_id)
        ).select()
        s_map_not_in_exomes = sample_map_ht.anti_join(
            exome_ht.key_by(i=exome_ht.s.split("_")[1])
        )

        logger.info(
            f"Total number of IDs in the sample map that are not in the exome data: {s_map_not_in_exomes.count()}..."
        )
        s_map_not_in_exomes.show(s_map_not_in_exomes.count())
        logger.info(
            f"Total number of IDs in the exome data that are not in the sample map: {s_exome_not_in_map.count()}..."
        )
        s_exome_not_in_map.show(s_exome_not_in_map.count())

        exome_ht.write(
            array_sample_map_ht_path(data_source, freeze), overwrite=args.overwrite
        )

    if args.compute_last_END_positions:
        mt = get_ukbb_data(data_source, freeze, raw=True, adj=False, split=False)
        last_END_positions_ht = compute_last_ref_block_end(mt)
        last_END_positions_ht.write(last_END_positions_ht_path(data_source, freeze))

    if args.compute_qc_mt:
        mt = get_ukbb_data(data_source, freeze, raw=True, adj=False, split=False, key_by_locus_and_alleles=True)
        mt = mt.select_entries("DP", "END", "LGT")
        logger.info(
            f"Total number of variants in raw unsplit matrix table: {mt.count_rows()}"
        )

        logger.info("Filtering to QC MT sites from tranche 2/freeze 5 and removing LCR intervals")
        if not hl.utils.hadoop_exists(f"{qc_sites_path()}/_SUCCESS"):
            get_qc_mt_sites()
        qc_sites_ht = hl.read_table(qc_sites_path())
        lcr = get_lcr_intervals()
        qc_sites_ht = qc_sites_ht.filter(
            hl.is_missing(lcr[qc_sites_ht.key])
        )
        mt = mt.annotate_rows(**sites_ht[mt.row_key])
        mt = mt.filter_rows(hl.is_defined(mt.info)
        
        logger.info("Splitting multiallelics and filtering to adj")
        mt = hl.experimental.sparse_split_multi(mt)
        mt = filter_to_adj(mt)

        mt = mt.checkpoint(
            get_checkpoint_path(
                data_source,
                freeze,
                name=f"{data_source}.freeze_{freeze}.rekey.biallelic.adj.lcr.mt",
                mt=True
            ),
            overwrite=True,
        )
        logger.info(
            f"Total number of variants after splitting and LCR+site filtering: {mt.count_rows()}"
        )

        if data_source == "regeneron":
            apply_hard_filters = False
        else:
            apply_hard_filters = True

        # NOTE: adding decoy and segdup hg38 resources/filters are still pending
        qc_mt = get_qc_mt(
            mt,
            min_af=args.min_af,
            min_callrate=args.min_callrate,
            apply_hard_filters=apply_hard_filters,
            ld_r2=None,
            filter_lcr=False,
            filter_decoy=False,
            filter_segdup=False,
        )
        qc_mt = qc_mt.naive_coalesce(args.n_partitions)
        qc_mt = qc_mt.checkpoint(
            qc_mt_path(data_source, freeze), overwrite=args.overwrite
        )
        logger.info(
            f"Total number of bi-allelic, high-callrate, common SNPs for sample QC: {qc_mt.count_rows()}"
        )

    if args.compute_sample_qc_ht:
        qc_mt = hl.read_matrix_table(qc_mt_path(data_source, freeze))
        qc_mt = filter_to_autosomes(qc_mt)
        qc_ht = hl.sample_qc(qc_mt).cols().select("sample_qc")
        qc_ht = qc_ht.transmute(
            sample_qc=qc_ht.sample_qc.select("call_rate", "gq_stats", "dp_stats")
        )
        qc_ht.write(qc_ht_path(data_source, freeze), overwrite=args.overwrite)

    if args.extract_truth_samples:
        mt = get_ukbb_data(data_source, freeze, raw=True, adj=False, split=False)
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
            mt.naive_coalesce(100).write(
                truth_samples[truth_sample]["mt"], overwrite=args.overwrite
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
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )

    parser.add_argument(
        "--create_exome_array_id_map_ht",
        help="Load exome to array id mapping file into hail Table",
        action="store_true",
    )
    parser.add_argument(
        "--compute_last_END_positions",
        help="Compute last END position for each line in sparse matrix table",
        action="store_true"
    )
    parser.add_argument(
        "--compute_qc_mt",
        help="Compute matrix to be used in sample qc",
        action="store_true",
    )
    parser.add_argument(
        "--min_callrate",
        help="Minimum variant callrate to retain variant in qc matrix table.",
        default=0.99,
        type=float,
    )
    parser.add_argument(
        "--min_af",
        help="Minimum variant allele frequency to retain variant in qc matrix table.",
        default=0.001,
        type=float,
    )
    parser.add_argument(
        "--pct_samples_20x",
        help="Percent samples at 20X to filter intervals",
        default=0.85,
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
        "--n_partitions", help="Desired number of partitions for output Table/MatrixTable",
        default=5000, type=int
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
