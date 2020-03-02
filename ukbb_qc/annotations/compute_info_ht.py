import argparse
import hail as hl
import logging
from gnomad_hail.utils.annotations import get_lowqual_expr
from gnomad_hail.utils.generic import bi_allelic_site_inbreeding_expr
from gnomad_hail.utils.slack import try_slack
from gnomad_hail.utils.sparse_mt import (
    default_compute_info,
    split_info_annotation,
    split_lowqual_annotation,
)
from ukbb_qc.resources.sample_qc import related_drop_path
from ukbb_qc.resources.variant_qc import info_ht_path


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("compute_info_ht")
logger.setLevel(logging.INFO)


def main(args):

    # NOTE: This script should be run after relatedness inference is complete
    hl.init(log="/compute_info_ht.log", default_reference="GRCh38")
    data_source = args.data_source
    freeze = args.freeze

    logger.info("Reading in raw MT and removing hard filtered samples...")
    mt = get_ukbb_data(
        data_source, freeze, split=False, raw=True, key_by_locus_and_alleles=True
    )
    mt = mt.filter_rows((hl.len(mt.alleles) > 1))
    mt = remove_hard_filter_samples(data_source, freeze, mt, mt.LGT, non_refs_only=True)

    logger.info("Removing related samples...")
    related_ht = related_drop_path(data_source, freeze)
    related_ht = related_ht.filter(
        (related_ht.relationship == "Unrelated")
    )
    mt = mt.filter_cols(hl.is_defined(related_ht[mt.col_key]))
    mt = mt.filter_rows(hl.agg.any(mt.LGT.is_non_ref()))

    logger.info("Adding inbreeding coefficient...")
    mt = mt.annotate_rows(InbreedingCoeff=bi_allelic_site_inbreeding_expr(mt.LGT))
    ic_ht = mt.rows().select("InbreedingCoeff")

    logger.info("Computing info HT...")
    info_ht = default_compute_info(mt)
    info_ht = info_ht.annotate(**ic_ht[info_ht.key])

    logger.info("Dropping default lowqual expr and reannotating...")
    info_ht = info_ht.drop("lowqual")
    info_ht = info_ht.annotate(
        lowqual=get_lowqual_expr(
            info_ht.alleles, info_ht.info.QUALapprox, indel_phred_het_prior=40
        )
    )
    info_ht = info_ht.checkpoint(info_ht_path(data_source, freeze, split=False))

    logger.info("Splitting info ht...")
    info_ht = hl.split_multi(info_ht)
    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(
            split_info_annotation(info_ht.info, info_ht.a_index)
        ),
        lowqual=split_lowqual_annotation(info_ht.lowqual, info_ht.a_index),
    )
    info_ht.write(info_ht_path(data_source, freeze, split=True))


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