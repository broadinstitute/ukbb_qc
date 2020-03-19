import argparse
import hail as hl
import logging
from gnomad.utils.slack import try_slack
from gnomad.utils.sparse_mt import compute_last_ref_block_end
from ukbb_qc.assessment.sanity_checks import sample_check, summarize_mt
from ukbb_qc.resources.basics import (
    array_sample_map_ht_path,
    capture_ht_path,
    get_ukbb_data,
    last_END_positions_ht_path,
    raw_mt_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.load_data.utils import (
    import_array_exome_id_map_ht,
    import_capture_intervals,
    import_phenotype_ht,
    import_vqsr,
)


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("load_data")
logger.setLevel(logging.INFO)


def main(args):
    hl.init(log="/load_data.log", default_reference="GRCh38")

    data_source = "broad"
    freeze = args.freeze

    if args.load_exome_array_id_map:
        logger.info("Loading array-exome sample ID map...")
        sample_map_ht = import_array_exome_id_map_ht(freeze)

        logger.info("Loading raw MatrixTable (to get exome IDs)...")
        exome_ht = hl.read_matrix_table(raw_mt_path(data_source, freeze)).cols()

        logger.info("Checking for sample discrepancies between MT and linking file...")
        exome_ht = exome_ht.key_by(eid_sample=exome_ht.s.split("_")[1])
        sample_check(exome_ht, sample_map_ht)

        exome_ht = exome_ht.annotate(**sample_map_ht[exome_ht.eid_sample])
        exome_ht = exome_ht.key_by("s")
        exome_ht.write(array_sample_map_ht_path(freeze), overwrite=args.overwrite)

    if args.load_phenotypes:
        logger.info("Importing phenotype data...")
        import_phenotype_ht()

    if args.load_capture_intervals:
        logger.info("Importing capture intervals...")
        import_capture_intervals(
            args.intervals, args.header, args.overwrite,
        )

    if args.load_vqsr:
        logger.info("Loading VQSR VCF...")
        import_vqsr(
            data_source,
            freeze,
            args.vqsr_vcf_path,
            args.vqsr_type,
            args.n_partitions,
            args.overwrite,
            args.header_path,
        )

    if args.sanity_check_raw_mt:
        logger.info("Reading in raw MT and summarizing variants...")
        mt = get_ukbb_data(
            data_source, freeze, raw=True, split=False, key_by_locus_and_alleles=True
        )
        var_summary = summarize_mt(mt)
        logger.info(f"Variant summary struct: {var_summary}")

        logger.info("Checking for duplicate samples...")
        s_list = mt.aggregate_cols(hl.agg.collect(mt.s))
        s_set = mt.aggregate_cols(hl.agg.collect_as_set(mt.s))

        if len(s_set) != len(s_list):
            logger.warning("There are duplicate sample IDs in the raw MT!")
            dups = []
            for s in s_set:
                if s_list.count(s) > 1:
                    dups.append(s)
            logger.warning(f"Duplicate list: {dups}")

    if args.compute_last_END_positions:
        logger.info("Computing last END position HT...")
        mt = get_ukbb_data(data_source, freeze, raw=True, adj=False, split=False)
        last_END_positions_ht = compute_last_ref_block_end(mt)
        last_END_positions_ht.write(last_END_positions_ht_path(freeze))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int,
    )

    parser.add_argument(
        "--load_exome_array_id_map",
        help="Load exome to array id mapping file into Table",
        action="store_true",
    )
    parser.add_argument(
        "--load_phenotypes", help="Load phenotype file into Table", action="store_true",
    )

    capture_intervals = parser.add_argument_group(
        "capture_intervals",
        description="Arguments relevant to loading capture intervals",
    )
    capture_intervals.add_argument(
        "--load_capture_intervals",
        help="Load capture intervals into Table",
        action="store_true",
    )
    capture_intervals.add_argument(
        "--intervals", help="Path to intervals file",
    )
    capture_intervals.add_argument(
        "--header",
        help="Whether intervals file contains a header line",
        action="store_true",
    )

    vqsr = parser.add_argument_group(
        "vqsr", description="Arguments relevant to loading VQSR",
    )
    vqsr.add_argument(
        "--load_vqsr", help="Load VQSR file into Table", action="store_true",
    )
    vqsr.add_argument("--vqsr_vcf_path", help="Path to VQSR VCF")
    vqsr.add_argument(
        "--vqsr_type",
        help="Path to VQSR VCF. Can be specified as Hadoop glob patterns",
        choices=["AS", "AS_TS"],
        default="AS",
    )
    vqsr.add_argument(
        "--n_partitions",
        help="Number of desired partitions for output Table",
        default=5000,
        type=int,
    )
    vqsr.add_argument(
        "--header_path",
        help="Optional path to a header file to use for importing VQSR VCF",
    )

    parser.add_argument(
        "--sanity_check_raw_mt",
        help="Sanity check raw MatrixTable",
        action="store_true",
    )
    parser.add_argument(
        "--compute_last_END_positions",
        help="Compute last END position for each line in sparse matrix table",
        action="store_true",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to",
    )
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
