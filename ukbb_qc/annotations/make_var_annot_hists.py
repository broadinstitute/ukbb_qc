import argparse
import logging

import hail as hl

from gnomad.utils.annotations import create_frequency_bins_expr, get_annotations_hists
from gnomad.utils.slack import slack_notifications
from ukbb_qc.resources.basics import (
    logging_path,
    release_ht_path,
    release_var_hist_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.slack_creds import slack_token


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variant_histograms")
logger.setLevel(logging.INFO)


ANNOTATIONS_HISTS = {
    "FS": (0, 50, 50),
    "InbreedingCoeff": (-0.25, 0.25, 50),
    "MQ": (0, 80, 40),
    "MQRankSum": (-15, 15, 60),
    "QD": (0, 40, 40),
    "ReadPosRankSum": (-15, 15, 60),
    "SOR": (0, 10, 50),
    #'BaseQRankSum': (-15, 15, 60),
    "VarDP": (1, 9, 32),
    "AS_VQSLOD": (-30, 30, 60),
    "rf_tp_probability": (0, 1, 50),
    "AS_pab_max": (0, 1, 50),
}


def main(args):

    hl.init(log="/variant_histograms.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)

    # NOTE: histogram aggregations on these metrics are done on the entire callset (not just PASS variants), on raw data
    metrics = list(ANNOTATIONS_HISTS.keys())

    try:
        logger.info("Reading in release HT")
        ht = hl.read_table(release_ht_path(*tranche_data)).select_globals()

        # Select annotations and move necessary annotations into info struct
        # get_annotations_hists checks within ht.info (should this be changed?)
        ht = ht.select(
            qual=ht.qual,
            freq=ht.freq,
            info=hl.struct(
                FS=ht.rf.info.FS,
                InbreedingCoeff=ht.info.InbreedingCoeff,
                MQ=ht.rf.info_MQ,
                MQRankSum=ht.info.MQRankSum,
                QD=ht.rf.info.QD,
                ReadPosRankSum=ht.info.ReadPosRankSum,
                SOR=ht.rf.info.SOR,
                VarDP=ht.rf.info.VarDP,
                rf_tp_probability=ht.rf.rf_probability,
                AS_pab_max=ht.info.AS_pab_max,
                AS_VQSLOD=hl.float(ht.vqsr.AS_VQSLOD),
            ),
        )

        logger.info("Getting info annotation histograms")
        hist_ranges_expr = get_annotations_hists(ht, ANNOTATIONS_HISTS, ["VarDP"])

        # NOTE: this is code was run on gnomAD v2 and has not been run since
        if args.first_pass:
            logger.info(
                "Evaluating minimum and maximum values for each metric of interest, and \
                caps maximum value at 1e10"
            )
            minmax_dict = {}
            for metric in metrics:
                minmax_dict[metric] = hl.struct(
                    min=hl.agg.min(ht[metric]),
                    max=hl.if_else(
                        hl.agg.max(ht[metric]) < 1e10, hl.agg.max(ht[metric]), 1e10
                    ),
                )
            minmax = ht.aggregate(hl.struct(**minmax_dict))
            logger.info(f"Metrics bounds: {minmax}")

        else:
            logger.info("Aggregating hists over hand-tooled ranges")
            hists = ht.aggregate(
                hl.array(
                    [
                        hist_expr.annotate(metric=hist_metric)
                        for hist_metric, hist_expr in hist_ranges_expr.items()
                    ]
                ).extend(
                    hl.array(
                        hl.agg.group_by(
                            create_frequency_bins_expr(
                                AC=ht.freq[1].AC, AF=ht.freq[1].AF
                            ),
                            # NOTE: used qual here for tranche 2, should we use QUALapprox for future tranches?
                            hl.agg.hist(hl.log10(ht.qual), 1, 10, 36),
                        )
                    ).map(lambda x: x[1].annotate(metric=x[0]))
                ),
                _localize=False,
            )

            logger.info("Writing output")
            with hl.hadoop_open(release_var_hist_path(*tranche_data), "w") as f:
                f.write(hl.eval(hl.json(hists)))
    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(*tranche_data))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--first_pass",
        help="Determine min/max values for each variant metric (to be used in hand-tooled histogram ranges)",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
