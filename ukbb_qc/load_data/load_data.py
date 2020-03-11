from gnomad.utils.slack import try_slack
from gnomad.utils.sparse_mt import compute_last_ref_block_end
from ukbb_qc.assessment.sanity_checks import sample_check, summarize_mt
from ukbb_qc.resources.basics import (
    array_sample_map_ht_path,
    capture_ht_path,
    last_END_positions_ht_path,
)
from ukbb_qc.load_data.utils import (
    import_array_exome_id_map_ht,
    import_capture_intervals,
    import_phenotype_ht,
    import_vqsr,
)
import argparse
import hail as hl
import logging


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("load_data")
logger.setLevel(logging.INFO)


def main(args):
    hl.init(log="/load_data.log", default_reference="GRCh38")

    data_source = args.data_source
    freeze = args.freeze

    if args.load_exome_array_id_map:
        logger.info("Loading array-exome sample ID map...")
        sample_map_ht = import_array_exome_id_map_ht(freeze)
        sample_map_ht.write(
            array_sample_map_ht_path(data_source, freeze), overwrite=args.overwrite
        )

    if args.load_phenotypes:
        logger.info("Importing phenotype data...")
        import_phenotype_ht()

    if args.load_capture_intervals:
        logger.info("Importing capture intervals...")
        import_capture_intervals(
            args.intervals, capture_ht_path(data_source), args.header, args.overwrite
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
        summarize_mt(mt)

        logger.info(
            "Checking for sample discrepancies between MatrixTable and linking file..."
        )
        sample_map_ht = hl.read_table(array_sample_map_ht_path(data_source, freeze))
        samples_ht = mt.cols()
        samples_ht = samples_ht.key_by(array_ID=samples_ht.s.split("_")[1])
        sample_check(samples_ht, sample_map_ht)

    if args.compute_last_END_positions:
        mt = get_ukbb_data(data_source, freeze, raw=True, adj=False, split=False)
        last_END_positions_ht = compute_last_ref_block_end(mt)
        last_END_positions_ht.write(last_END_positions_ht_path(data_source, freeze))


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
        "--load_exome_array_id_map",
        help="Load exome to array id mapping file into Table",
        action="store_true",
    )
    parser.add_argument(
        "--load_phenotypes", help="Load phenotype file into Table", action="store_true",
    )

    parser.add_argument_group(
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

    parser.add_argument_group(
        "vqsr", description="Arguments relevant to loading VQSR",
    )
    vqsr.add_argument(
        "--load_vqsr", help="Load VQSR file into Table", action="store_true",
    )
    vqsr.add_argument("--vqsr_vcf_path", help="Path to VQSR VCF")
    vqsr.add_argument(
        "--vqsr_type", help="Path to VQSR VCF", choices=["AS", "AS_TS"], default="AS",
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
