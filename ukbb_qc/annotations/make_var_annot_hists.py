from gnomad_hail.utils.slack import try_slack
from ukbb_qc.resources.basics import release_ht_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
import argparse
import hail as hl
import logging
import sys
from typing import Dict


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variant_histograms")
logger.setLevel(logging.INFO)


def main(args):

    hl.init(log="/variant_histograms.log", default_reference="GRCh38")
    data_source = args.data_source
    freeze = args.freeze

    metrics = [
        "FS",
        "InbreedingCoeff",
        "MQ",
        "MQRankSum",
        "QD",
        "ReadPosRankSum",
        "SOR",
        # "BaseQRankSum" NOTE: this was removed from vcf export for ukbb and therefore release mt/ht
        "VarDP",
        "rf_probability",
        "pab_max",
        "AS_VQSLOD",
    ]

    ht = hl.read_table(release_ht_path(data_source, freeze))

    # Move necessary annotations out of nested structs
    ht = ht.transmute(
        FS=ht.rf.info_FS,
        InbreedingCoeff=ht.rf.inbreeding_coeff,
        MQ=ht.rf.info_MQ,
        MQRankSum=ht.rf.info_MQRankSum,
        QD=ht.rf.info_QD,
        ReadPosRankSum=ht.rf.info_ReadPosRankSum,
        SOR=ht.rf.info_SOR,
        VarDP=ht.rf.info_VarDP,
        rf_probability=ht.rf.rf_probability,
        pab_max=ht.rf.pab_max,
        AS_VQSLOD=ht.vqsr.AS_VQSLOD,
    )
    ht = ht.explode(ht.AS_VQSLOD)
    ht = ht.transmute(AS_VQSLOD=hl.float(ht.AS_VQSLOD))
    # NOTE: histogram aggregations are done on the entire callset (not just PASS variants), on raw data

    # NOTE: run the following code in a first pass to determine bounds for metrics
    if args.first_pass:
        logger.info("Evaluating minimum and maximum values for each metric of interest")
        minmax_dict = {}
        for metric in metrics:
            minmax_dict[metric] = hl.struct(
                min=hl.agg.min(ht[metric]),
                max=hl.cond(
                    hl.agg.max(ht[metric]) < 1e10, hl.agg.max(ht[metric]), 1e10
                ),
            )
        minmax = ht.aggregate(hl.struct(**minmax_dict))
        logger.info(f"Metrics bounds: {minmax}")
    else:
        logger.info("Aggregating hists over hand-tooled ranges")
        hists = ht.aggregate(hl.struct(**define_hist_ranges(ht)))
        hist_out = hl.array(
            [hists[f"{metric}"].annotate(metric=metric) for metric in metrics]
        )

        # Aggregate QUAL stats by bin:
        logger.info("Aggregating QUAL stats by bin")
        ht = aggregate_qual_stats_by_bin(ht)
        ht = ht.group_by("metric").aggregate(
            hist=hl.agg.hist(hl.log(ht.qual, base=10), 1, 10, 36)
        )
        hists = ht.collect()
        hist_out = hist_out.extend(
            hl.array([x.hist.annotate(metric=x.metric) for x in hists])
        )

        with hl.hadoop_open(release_var_hist_path(data_type), "w") as f:
            f.write(hl.eval(hl.json(hist_out)))


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
        "--first_pass",
        help="Determine min/max values for each variant metric (to be used in hand-tooled histogram ranges",
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
