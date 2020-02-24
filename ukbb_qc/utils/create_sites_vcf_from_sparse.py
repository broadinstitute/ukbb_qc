import hail as hl
from gnomad_hail import *
from gnomad_hail.utils.generic import *
from gnomad_hail.utils.sparse_mt import *
from ukbb_qc.resources import *

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_sites_vcf")
logger.setLevel(logging.INFO)


def main(args):

    hl.init(default_reference="GRCh38", log="/create_sites_vcf.log")
    freeze = args.freeze
    temp_prefix = f"gs://broad-ukbb/broad.freeze_{freeze}/temp"

    logger.info("Reading in input mt (raw sparse mt)")
    mt = hl.read_matrix_table(args.input)
    logger.info(f"Sparse mt count: {mt.count()}")
    logger.info("Filtering out ref blocks")
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)
    logger.info(f"Count after filtration: {mt.count()}")

    if args.test:
        logger.info("Trimming to chr20")
        mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chr20")])
        logger.info(f"Count after filtration: {mt.count()}")

    logger.info("Getting site info")
    mt = mt.annotate_rows(info=get_site_info_expr(mt))

    logger.info("Getting allele specific info")
    mt = mt.annotate_rows(as_info=get_as_info_expr(mt))

    logger.info("Merging info fields")
    mt = mt.transmute_rows(info=mt.info.annotate(**mt.as_info))
    mt.describe()

    # NOTE: Annotations do NOT need densify
    vcf_mt = ht_to_vcf_mt(mt.rows())
    hl.export_vcf(vcf_mt, f"{temp_prefix}/broad.freeze_{freeze}.sites.vcf.bgz")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Path to input mt", required=True)
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument("--test", help="Filter to chr20 to test", action="store_true")
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
