import argparse
import hail as hl
import logging
from gnomad_hail.utils.generic import filter_to_autosomes
from gnomad_hail.utils.gnomad_functions import filter_to_adj
from gnomad_hail.utils.sample_qc import get_qc_mt
from gnomad_hail.utils.sparse_mt import compute_last_ref_block_end, densify_sites
from ukbb_qc.resources.basics import get_checkpoint_path, get_ukbb_data
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import array_sample_map_ht_path, qc_mt_path
from ukbb_qc.resources.variant_qc import get_truth_sample_info
from ukbb_qc.sample_qc.apply_hard_filters import hard_filter_samples
from ukbb_qc.utils.utils import annotate_interval_qc_filter, get_qc_mt_sites


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("process_raw_data")
logger.setLevel(logging.INFO)


def main(args):
    hl.init(log="/process_raw_data.log", default_reference="GRCh38")

    data_source = args.data_source
    freeze = args.freeze

    if args.create_exome_array_id_map_ht:
        logger.info(
            "Checking for sample discrepancies between MatrixTable and linking file..."
        )
        sample_map_ht = import_array_exome_id_map_ht(data_source, freeze)
        exome_ht = get_ukbb_data(
            data_source, freeze, raw=True, split=False, key_by_locus_and_alleles=True
        ).cols()
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

    if args.apply_hard_filters:
        logger.info("Reading in raw MT...")
        mt = get_ukbb_data(
            data_source, freeze, raw=True, adj=False, key_by_locus_and_alleles=True
        )
        mt = mt.annotate_entries(GT=hl.experimental.lgt_to_gt(mt.LGT, mt.LA))
        mt = mt.select_entries("DP", "END", "GT")
        logger.info(
            f"Total number of variants in raw unsplit matrix table: {mt.count_rows()}"
        )

        last_END_ht = hl.read_table(last_END_positions_ht_path(data_source, freeze))
        sex_ht = hl.read_table(sex_ht_path(data_source, freeze))
        interval_qc_ht = hl.read_table(
            interval_qc_path(data_source, freeze, "autosomes")
        )
        intervals_qc_ht = interval_qc_ht.filter(
            interval_qc_ht[args.cov_filter_field] > args.pct_samples
        )

        # Rename interval so that there are no name collisions when calling densify_sites
        intervals_qc_ht = interval_qc_ht.key_by()
        intervals_qc_ht = interval_qc_ht.transmute(
            qc_interval=interval_qc.interval
        ).key_by("qc_interval")

        hard_filters_ht = hard_filter_samples(
            data_source,
            freeze,
            mt,
            interval_qc_ht,
            last_END_ht,
            sex_ht,
            args.min_callrate,
            args.min_dp
        )
        ht = ht.naive_coalesce(args.n_partitions)
        ht = ht.checkpoint(
            hard_filters_ht_path(data_source, freeze), overwrite=args.overwrite
        )
        logger.info("Checking number of samples flagged with hard filters...")
        ht = ht.explode(ht.hard_filters)
        filters = ht.aggregate(hl.agg.counter(ht.hard_filters))
        for filt in filters:
            logger.info(f"Samples flagged due to {filt}: {filters[filt]}")

    if args.compute_qc_mt:
        logger.info("Reading in raw MT...")
        mt = get_ukbb_data(
            data_source, freeze, raw=True, adj=False, key_by_locus_and_alleles=True
        )
        mt = mt.annotate_entries(GT=hl.experimental.lgt_to_gt(mt.LGT, mt.LA))
        mt = mt.select_entries("DP", "END", "GT")
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
            get_checkpoint_path(data_source, freeze, name="dense", mt=True),
            overwrite=True,
        )

        logger.info("Repartitioning densified MT")
        mt = mt.naive_coalesce(n_partitions)
        mt = mt.checkpoint(
            get_checkpoint_path(
                data_source, freeze, name="dense.repartitioned", mt=True
            ),
            overwrite=True,
        )
        logger.info(
            "Adding info annotations, filtering low QUAL variants, and filtering to adj..."
        )
        # NOTE: Need MQ, QD, FS for hard filters
        info_expr = get_site_info_expr(mt)
        info_expr = info_expr.annotate(**get_as_info_expr(mt))
        mt = mt.annotate_rows(info=info_expr)
        mt = mt.annotate_rows(
            lowqual=get_lowqual_expr(
                mt.alleles, mt.info.QUALapprox, indel_phred_het_prior=40
            )
        )
        mt = mt.filter_rows(~mt.lowqual)
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
        qc_mt = qc_mt.checkpoint(
            qc_mt_path(data_source, freeze, ld_pruned=True), overwrite=args.overwrite
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
            mt.naive_coalesce(args.n_partitions).write(
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
        action="store_true",
    )
    parser.add_argument(
        "--apply_hard_filters".help="Apply hard filters to samples", aciton="store_true"
    )
    parser.add_argument("--min-dp", help="Minimum depth", default=20.0, type=float)
    parser.add_argument(
        "--compute_qc_mt",
        help="Compute matrix to be used in sample qc",
        action="store_true",
    )
    parser.add_argument(
        "--apply_interval_qc_filter",  # NOTE not super sure this is necessary
        help="Filter to good intervals from interval QC",
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
        "--densified_partitions",
        help="Desired number of partitions for densified MatrixTable",
        default=25000,
        type=int,
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
        "--slack_channel", help="Slack channel to post results and notifications to."
    )

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
