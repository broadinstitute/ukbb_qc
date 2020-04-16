import argparse
import logging

import hail as hl

from gnomad.utils.annotations import get_lowqual_expr
from gnomad.utils.slack import try_slack
from gnomad.utils.sparse_mt import (
    default_compute_info,
    split_info_annotation,
    split_lowqual_annotation,
)
from ukbb_qc.resources.basics import get_ukbb_data
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.variant_qc import info_ht_path, var_annotations_ht_path
from ukbb_qc.utils.utils import annotate_interval_qc_filter, remove_hard_filter_samples


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("compute_info_ht")
logger.setLevel(logging.INFO)


def main(args):
    # NOTE: This script should be run after hard filtering
    hl.init(log="/compute_info_ht.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze

    logger.info("Reading in raw MT and removing hard filtered samples...")
    mt = get_ukbb_data(
        data_source,
        freeze,
        split=False,
        raw=True,
        repartition=args.repartition,
        n_partitions=args.raw_partitions,
        key_by_locus_and_alleles=True,
    )
    mt = mt.filter_rows((hl.len(mt.alleles) > 1))
    mt = remove_hard_filter_samples(data_source, freeze, mt, non_refs_only=True)

    if args.use_vqsr:
        logger.info("Reading in VQSR unsplit HT...")
        info_ht = hl.read_table(
            var_annotations_ht_path(
                f'{"vqsr" if args.vqsr_type == "AS" else "AS_TS_vqsr"}.unsplit',
                data_source,
                freeze,
            )
        )

        logger.info("Annotating raw MT with pab max...")
        # Note: we use indel_phred_het_prior=40 to be more consistent with the filtering used by DSP/Laura for VQSR
        mt = mt.select_rows()
        ht = mt.annotate_rows(
            AS_pab_max=hl.agg.array_agg(
                lambda ai: hl.agg.filter(
                    mt.LA.contains(ai) & mt.LGT.is_het(),
                    hl.agg.max(
                        hl.binom_test(mt.LAD[1], hl.sum(mt.LAD), 0.5, "two-sided")
                    ),
                ),
                hl.range(1, hl.len(mt.alleles)),
            )
        ).rows()

        logger.info("Annotating VQSR unsplit HT with pab max from raw MT...")
        info_ht = info_ht.annotate(AS_pab_max=ht[info_ht.key].AS_pab_max)
    else:
        logger.info("Computing info HT...")
        info_ht = default_compute_info(mt, site_annotations=True)

    info_ht = info_ht.checkpoint(
        info_ht_path(data_source, freeze, split=False), overwrite=args.overwrite
    )

    logger.info("Annotating info HT with lowqual...")
    # Note: we use indel_phred_het_prior=40 to be more consistent with the filtering used by DSP/Laura for VQSR
    info_ht = info_ht.annotate(
        lowqual=get_lowqual_expr(
            info_ht.alleles, info_ht.info.QUALapprox, indel_phred_het_prior=40
        ),
        AS_lowqual=get_lowqual_expr(
            info_ht.alleles, info_ht.info.AS_QUALapprox, indel_phred_het_prior=40
        ),
    )

    logger.info("Splitting info ht...")
    info_ht = hl.split_multi(info_ht)
    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(
            **split_info_annotation(info_ht.info, info_ht.a_index),
        ),
        AS_lowqual=split_lowqual_annotation(info_ht.AS_lowqual, info_ht.a_index),
    )

    info_ht = annotate_interval_qc_filter(
        data_source,
        freeze,
        info_ht,
        cov_filter_field=args.interval_cov_filter_field,
        xy_cov_filter_field=args.xy_cov_filter_field,
        pct_samples=args.interval_filter_pct_samples,
    )
    info_ht = info_ht.annotate_globals(
        vqsr=args.use_vqsr,
        vqsr_type=args.vqsr_type if args.use_vqsr else None,
        cov_filter_field=args.interval_cov_filter_field,
        xy_cov_filter_field=args.xy_cov_filter_field,
        pct_samples=args.interval_filter_pct_samples,
    )

    info_ht.write(
        info_ht_path(data_source, freeze, split=True), overwrite=args.overwrite
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--repartition",
        help="Repartition raw MT on read. Required for 300K",
        action="store_true",
    )
    parser.add_argument(
        "--raw_partitions",
        help="Number of desired partitions for the raw MT. Necessary only for 300K. Used only if --repartition is also specified",
        default=30000,
        type=int,
    )
    parser.add_argument(
        "--use_vqsr",
        help="Will use the VQSR HT instead of the info_ht",
        action="store_true",
    )
    parser.add_argument(
        "--vqsr_type",
        help="What type of VQSR was run: allele-specific, or allele-specific with transmitted singletons",
        type=str,
        choices=["AS", "AS_TS"],
        default="AS",
    )
    parser.add_argument(
        "--interval_cov_filter_field",
        help="Field used to filter intervals on autosomes and par regions",
        default="pct_samples_20x",
    )
    parser.add_argument(
        "--xy_cov_filter_field",
        help="Field used to filter intervals on sex chromosomes",
        default="pct_samples_10x",
    )
    parser.add_argument(
        "--interval_filter_pct_samples",
        help="Percent samples at specified coverage field to filter intervals",
        default=0.85,
        type=float,
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
