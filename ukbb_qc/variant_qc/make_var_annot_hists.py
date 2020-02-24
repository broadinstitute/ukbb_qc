from gnomad_hail.utils.slack import try_slack
from ukbb_qc.resources.basics import release_ht_path
import argparse
import hail as hl
import logging
import sys
from typing import Dict


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variant_histograms")
logger.setLevel(logging.INFO)


def define_hist_ranges(ht: hl.Table) -> Dict[str: hl.expr.StructExpression]:
    """
    Creates histograms for relevant metrics

    :param Table ht: Table with variant metrics
    :return: Dictionary of metrics and their histograms
    :rtype: Dict[str: hl.expr.StructExpression]
    """
    hist_dict = {
        "FS": hl.agg.hist(ht.FS, 0, 50, 50),  # NOTE: in 2.0.2 release this was on (0,20)
        "InbreedingCoeff": hl.agg.hist(ht.inbreeding_coeff, -0.25, 0.25, 50),
        "MQ": hl.agg.hist(ht.MQ, 0, 80, 40),
        "MQRankSum": hl.agg.hist(ht.MQRankSum, -15, 15, 60),
        "QD": hl.agg.hist(ht.QD, 0, 40, 40),
        "ReadPosRankSum": hl.agg.hist(ht.ReadPosRankSum, -15, 15, 60),
        "SOR": hl.agg.hist(ht.SOR, 0, 10, 50),
        #"BaseQRankSum": hl.agg.hist(ht.BaseQRankSum, -15, 15, 60),
        "VarDP": hl.agg.hist(hl.log(ht.VarDP, base=10), 1, 9, 32),  # NOTE: in 2.0.2 release this was on (0,8)
        "AS_VQSLOD": hl.agg.hist(ht.AS_VQSLOD, -30, 30, 60),  # NOTE: in 2.0.2 release this was on (-20,20)
        "rf_probability": hl.agg.hist(ht.rf_probability, 0, 1, 50),
        "pab_max": hl.agg.hist(ht.pab_max, 0, 1, 50)
    }
    return hist_dict


def aggregate_qual_stats_by_bin(ht: hl.Table) -> hl.Table:
    """
    Creates bins for frequencies in preparation for aggregating QUAL by bin

    :param Table ht: Table with variant frequencies
    :return: Table with bin annotation
    :rtype: hl.Table
    """
    # NOTE: freq[1] is raw frequencies
    ht = ht.annotate(metric=(hl.case()
                             .when(ht.freq[1].AC== 1, "binned_singleton")
                             .when(ht.freq[1].AC == 2, "binned_doubleton")
                             .when((ht.freq[1].AC > 2) & (ht.freq[1].AF < 0.00005), "binned_0.00005")
                             .when((ht.freq[1].AF >= 0.00005) & (ht.freq[1].AF < 0.0001), "binned_0.0001")
                             .when((ht.freq[1].AF >= 0.0001) & (ht.freq[1].AF < 0.0002), "binned_0.0002")
                             .when((ht.freq[1].AF >= 0.0002) & (ht.freq[1].AF < 0.0005), "binned_0.0005")
                             .when((ht.freq[1].AF >= 0.0005) & (ht.freq[1].AF < 0.001), "binned_0.001")
                             .when((ht.freq[1].AF >= 0.001) & (ht.freq[1].AF < 0.002), "binned_0.002")
                             .when((ht.freq[1].AF >= 0.002) & (ht.freq[1].AF < 0.005), "binned_0.005")
                             .when((ht.freq[1].AF >= 0.005) & (ht.freq[1].AF < 0.01), "binned_0.01")
                             .when((ht.freq[1].AF >= 0.01) & (ht.freq[1].AF < 0.02), "binned_0.02")
                             .when((ht.freq[1].AF >= 0.02) & (ht.freq[1].AF < 0.05), "binned_0.05")
                             .when((ht.freq[1].AF >= 0.05) & (ht.freq[1].AF < 0.1), "binned_0.1")
                             .when((ht.freq[1].AF >= 0.1) & (ht.freq[1].AF < 0.2), "binned_0.2")
                             .when((ht.freq[1].AF >= 0.2) & (ht.freq[1].AF < 0.5), "binned_0.5")
                             .when((ht.freq[1].AF >= 0.5) & (ht.freq[1].AF <= 1), "binned_1")
                             .default(hl.null(hl.tstr))))
    return ht


def main(args):
    hl.init(log="/variant_histograms.log", default_reference="GRCh38")

    metrics = [
                "FS", "InbreedingCoeff", "MQ", "MQRankSum", "QD", "ReadPosRankSum", "SOR", 
                # "BaseQRankSum" NOTE: this was removed from vcf export for ukbb and therefore release mt/ht
                "VarDP", "rf_probability", "pab_max", "AS_VQSLOD"
    ]

    ht = hl.read_table(release_ht_path(data_source, freeze))

    # Move necessary annotations out of nested structs
    ht = ht.transmute(
        FS=ht.rf.info_FS, InbreedingCoeff=ht.rf.inbreeding_coeff, MQ=ht.rf.info_MQ,
        MQRankSum=ht.rf.info_MQRankSum, QD=ht.rf.info_QD, ReadPosRankSum=ht.rf.info_ReadPosRankSum,
        SOR=ht.rf.info_SOR, VarDP=ht.rf.info_VarDP, rf_probability=ht.rf.rf_probability,
        pab_max=ht.rf.pab_max, AS_VQSLOD=ht.vqsr.AS_VQSLOD
    )
    
    # NOTE: histogram aggregations are done on the entire callset (not just PASS variants), on raw data

    # NOTE: run the following code in a first pass to determine bounds for metrics
    if args.first_pass:
        logger.info("Evaluating minimum and maximum values for each metric of interest")
        minmax_dict = {}
        for metric in metrics:
            minmax_dict[metric] = hl.struct(
                min=hl.agg.min(ht[metric]), 
                max=hl.cond(hl.agg.max(ht[metric])<1e10, hl.agg.max(ht[metric]), 1e10)
            )
        minmax = ht.aggregate(hl.struct(**minmax_dict))
        logger.info(f"Metrics bounds: {minmax}")
    else:
        logger.info("Aggregating hists over hand-tooled ranges")
        hists = ht.aggregate(hl.struct(**define_hist_ranges(ht)))
        hist_out = hl.array([hists[f"{metric}"].annotate(metric=metric) for metric in metrics])

        # Aggregate QUAL stats by bin:
        logger.info("Aggregating QUAL stats by bin")
        ht = aggregate_qual_stats_by_bin(ht)
        ht = ht.group_by("metric").aggregate(hist=hl.agg.hist(hl.log(ht.qual, base=10), 1, 10, 36))
        hists = ht.collect()
        hist_out = hist_out.extend(hl.array([x.hist.annotate(metric=x.metric) for x in hists]))

        with hl.hadoop_open(release_var_hist_path(data_type), "w") as f:
            f.write(hl.eval(hl.json(hist_out)))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--data_source", help="Data source", 
        choices=["regeneron", "broad"], default="broad"
    )
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--first_pass", 
        help="Determine min/max values for each variant metric (to be used in hand-tooled histogram ranges", 
        action="store_true"
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
