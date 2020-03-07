import argparse
import logging
import hail as hl
from gnomad.utils.annotations import annotate_freq
from gnomad.utils.slack import try_slack
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.basics import get_ukbb_data, array_sample_map_ht_path
from ukbb_qc.resources.variant_qc import var_annotations_ht_path


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("generate_frequency_data")
logger.setLevel(logging.INFO)


def generate_frequency_data(
    mt: hl.MatrixTable, calculate_by_platform: bool = False,
) -> hl.Table:
    """
    Generates frequency struct annotation containing AC, AF, AN, and homozygote count for dataset stratified by population. Optional to stratify by tranche

    :param MatrixTable mt: Input MatrixTable
    :param bool calculate_by_platform: Calculate frequencies per tranche
    :return: Table with frequency annotations in struct named `freq` and metadata in globals named `freq_meta`
    :rtype: Table
    """
    if calculate_by_platform:
        additional_strata_expr = {"tranche": mt.tranche}
    else:
        additional_strata_expr = None

    logger.info("Generating frequency data...")
    mt = annotate_freq(
        mt,
        sex_expr=mt.meta.sex,
        pop_expr=mt.meta.hybrid_pop,
        additional_strata_expr=additional_strata_expr,
    )
    return mt.rows().select("freq")


def main(args):
    hl.init(log="/frequency_generation.log")

    data_source = args.data_source
    freeze = args.freeze
    pops = hl.literal(set(args.pops.split(",")))
    logger.info(f"Filtering to {pops} for freq calculations")

    mt = get_ukbb_data(data_source, freeze, meta_root="meta")
    logger.info(
        f"mt count before filtering out low quality and non-releasable samples: {mt.count()}"
    )
    mt = mt.filter_cols(mt.meta.releasable & mt.meta.high_quality)
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    logger.info(
        f"mt count after filtering out low quality samples and their variants: {mt.count()}"
    )

    mt = mt.filter_cols(pops.contains(mt.meta.hybrid_pop))
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    logger.info(f"mt count after filtering to input pops: {mt.count()}")

    if args.by_platform:
        logger.info(
            "Annotating sample with tranche information (from array sample map ht)"
        )
        sample_map_ht = hl.read_table(array_sample_map_ht_path(data_source, freeze)).select(
            "batch.c"
        )
        mt = mt.annotate_cols(**sample_map_ht[mt.s])
        mt = mt.transmute_cols(tranche=mt["batch.c"])

    if args.calculate_frequencies:
        logger.info("Calculating frequencies")
        ht = generate_frequency_data(mt, args.by_platform)
        ht.write(
            var_annotations_ht_path(data_source, freeze, "cohort_freq"), args.overwrite
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
    parser.add_argument("--pops", help="Pops to include", default="nfe")
    parser.add_argument(
        "--calculate_frequencies",
        help="Calculate most frequency data",
        action="store_true",
    )
    parser.add_argument(
        "--by_platform",
        help="Also calculate frequencies by platform",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
