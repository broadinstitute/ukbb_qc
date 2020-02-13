import argparse
import hail as hl
import logging
from gnomad_hail.utils.sample_qc import *
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE, get_checkpoint_path, get_ukbb_data
from ukbb_qc.resources.sample_qc import qc_mt_path 


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("process_raw")
logger.setLevel(logging.INFO)


def ld_prune_qc_mt(qc_mt: hl.MatrixTable, ld_r2: float = 0.1) -> hl.MatrixTable:
    """
    Returns MatrixTable for sample QC purposes
    :param MatrixTable mt: Raw MatrixTable to be filtered
    :param float ld_r2: r2 to LD prune the QC matrix to
    :return: LD pruned MatrixTable for sample QC purposes
    :rtype: MatrixTable
    """
    pruned_ht = hl.ld_prune(qc_mt.GT, r2=ld_r2)
    qc_mt = qc_mt.filter_rows(hl.is_defined(pruned_ht[qc_mt.row_key]))
    qc_mt = qc_mt.annotate_globals(
        qc_mt_params=qc_mt.qc_mt_params.annotate(ld_r2=ld_r2)
    )
    qc_mt.annotate_cols(sample_callrate=hl.agg.fraction(hl.is_defined(qc_mt.GT)))
    return qc_mt


def main(args):
    hl.init(log="/process_raw.log")

    data_source = args.data_source
    freeze = args.freeze

    if args.create_exome_array_id_map_ht:
        sample_map = hl.import_table(array_sample_map(freeze), delimiter=",", quote='"')
        logger.info(
            f"Total number of IDs in the array to exome sample map: {sample_map.count()}..."
        )

        sample_map = sample_map.key_by(s=sample_map.eid_sample)
        exome_ht = get_ukbb_data(data_source, freeze, raw=True, split=False).cols()
        exome_ht = exome_ht.annotate(
            ukbb_app_26041_id=sample_map[exome_ht.s.split("_")[1]].eid_26041,
            **sample_map[exome_ht.s.split("_")[1]],
        )
        logger.info(f"Total number of samples in the exome data: {exome_ht.count()}...")

        s_exome_not_in_map = exome_ht.filter(
            hl.is_missing(exome_ht.ukbb_app_26041_id)
        ).select()
        s_map_not_in_exomes = sample_map.anti_join(
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
            array_sample_map_ht(data_source, freeze), overwrite=args.overwrite
        )

    if args.compute_qc_mt:
        logger.info(
            "Filtering to bi-allelic, high-callrate, common SNPs for sample QC..."
        )
        mt = get_ukbb_data(data_source, freeze, raw=True, adj=False, split=False)
        mt = mt.filter_rows(
            (hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1])
        )
        mt = mt.key_rows_by("locus", "alleles")
        mt = mt.checkpoint(
            get_mt_checkpoint_path(
                data_source,
                freeze,
                f"{data_source}.freeze_{freeze}.rekey.biallelic.GT.mt",
            )
        )
        logger.info(
            f"Total number of variants in raw unsplit matrix table: {mt.count_rows()}"
        )
        mt = hl.experimental.sparse_split_multi(mt)
        mt = filter_to_adj(mt)

        lcr = get_lcr_intervals()
        mt = mt.filter_rows(hl.is_missing(lcr[mt.locus]))

        if args.interval_qc_filter:
            mt = annotate_interval_qc_filter(data_source, freeze, mt, autosomes_only=True)
            mt = mt.filter_rows(
                mt.interval_qc_pass
            )
        if args.filter_sites:
            qc_mt = hl.read_matrix_table(qc_mt_path(data_source, freeze=5, ld_pruned=True))
            qc_sites = qc_mt.rows().select()
            mt = mt.filter_rows(hl.is_defined(qc_sites[mt.row_key]))

        mt = mt.checkpoint(
            get_mt_checkpoint_path(
                data_source,
                freeze,
                f"{data_source}.freeze_{freeze}.rekey.biallelic.adj.lcr.mt",
            ),
            overwrite=True,
        )
        logger.info(
            f"Total number of variants after LCR and interval/site filtering: {mt.count_rows()}"
        )

        # TODO change this for later tranches
        # sites_ht = hl.import_vcf('gs://broad-ukbb/broad.freeze_5/temp/broad.freeze_5.sites.vcf.bgz',reference_genome='GRCh38')
        sites_ht = hl.read_matrix_table(
            "gs://broad-ukbb/broad.freeze_5/temp/broad.freeze_5.sites.ht"
        ).rows()

        mt = mt.drop(mt.gvcf_info)
        mt = mt.annotate_rows(info=sites_ht[mt.row_key].info)

        if data_source == "regeneron":
            apply_hard_filters = False
        else:
            apply_hard_filters = True

        # Todo: changed to use Laurent's code for the qc_mt, change is that it uses inbreeding_coeff and hardy_weinberg_threshold, is this OK and are defaults OK?
        # Todo: add segdup and decoy filter
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
        qc_mt = qc_mt.naive_coalesce(5000)
        qc_mt = qc_mt.checkpoint(
            qc_mt_path(data_source, freeze), overwrite=args.overwrite
        )
        logger.info(
            f"Total number of variants in bi-allelic, high-callrate, common SNPs for sample QC: {qc_mt.count_rows()}"
        )

    if args.compute_sample_qc_ht:
        qc_mt = hl.read_matrix_table(qc_mt_path(data_source, freeze))
        qc_mt = filter_to_autosomes(qc_mt)
        qc_ht = hl.sample_qc(qc_mt).cols().select("sample_qc")
        qc_ht = qc_ht.transmute(
            sample_qc=qc_ht.sample_qc.select("call_rate", "gq_stats", "dp_stats")
        )
        qc_ht.write(qc_ht_path(data_source, freeze), overwrite=args.overwrite)

    if args.ld_prune_qc_mt:
        logger.info("LD prune qc MatrixTable sample QC...")
        qc_mt = hl.read_matrix_table(qc_mt_path(data_source, freeze))
        qc_mt = qc_mt.unfilter_entries()
        pruned_qc_mt = ld_prune_qc_mt(qc_mt, ld_r2=args.ld_r2)
        pruned_qc_mt = pruned_qc_mt.checkpoint(
            qc_mt_path(data_source, freeze, ld_pruned=True), overwrite=args.overwrite
        )
        logger.info(f"{pruned_qc_mt.count_rows()} variants found in LD-pruned MT")

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
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )

    parser.add_argument(
        "--create_exome_array_id_map_ht",
        help="Load exome to array id mapping file into hail Table",
        action="store_true",
    )
    parser.add_argument(
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
        "--interval_qc_filter",
        help="Should interval QC be applied",
        action="store_true",
    )
    parser.add_argument(
        "--filter_sites",
        help="Whether to filter to QC sites (from tranche 2 QC mt)",
        action="store_true",
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
        "--ld_prune_qc_mt", help="LD prunes the qc matrix", action="store_true"
    )
    parser.add_argument(
        "--ld_r2", help="r2 to LD prune the QC matrix to.", default=0.1, type=float
    )

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
