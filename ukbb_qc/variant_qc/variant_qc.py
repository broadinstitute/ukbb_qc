import argparse
import json
import logging
import sys
from typing import Dict, List
import uuid

import hail as hl

from gnomad.utils.file_utils import file_exists
from gnomad.utils.filtering import add_filters_expr
from gnomad.utils.slack import slack_notifications
from gnomad.variant_qc.pipeline import train_rf_model
from gnomad.variant_qc.random_forest import (
    apply_rf_model,
    load_model,
    median_impute_features,
    pretty_print_runs,
    save_model,
)
from ukbb_qc.resources.basics import get_ukbb_data
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.variant_qc import (
    info_ht_path,
    rf_run_hash_path,
    rf_path,
    rf_annotated_path,
    score_quantile_bin_path,
    var_annotations_ht_path,
)
from ukbb_qc.slack_creds import slack_token

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("run_ukbb_variant_qc")
logger.setLevel(logging.INFO)


LABEL_COL = "rf_label"
TRAIN_COL = "rf_train"
PREDICTION_COL = "rf_prediction"
INFO_FEATURES = ["AS_QD", "AS_pab_max", "AS_MQRankSum", "AS_SOR", "AS_ReadPosRankSum"]
FEATURES = [
    "InbreedingCoeff",
    "variant_type",
    "allele_type",
    "n_alt_alleles",
    "was_mixed",
    "has_star",
    "AS_QD",
    "AS_pab_max",
    "AS_MQRankSum",
    "AS_SOR",
    "AS_ReadPosRankSum",
]
TRUTH_DATA = ["hapmap", "omni", "mills", "kgp_phase1_hc", "ukbb_array_con_common"]

INBREEDING_COEFF_HARD_CUTOFF = -0.3


def create_rf_ht(
    data_source: str,
    freeze: int,
    impute_features_by_variant_type: bool = True,
    group: str = "raw",
) -> hl.Table:
    """
    Creates a Table with all necessary annotations for RF.

    Annotations that are added:

        Features for RF:
            - InbreedingCoeff
            - variant_type
            - allele_type
            - n_alt_alleles
            - was_mixed
            - has_star
            - AS_QD
            - AS_pab_max
            - AS_MQRankSum
            - AS_SOR
            - AS_ReadPosRankSum

        Training sites (bool):
            - transmitted_singleton
            - sibling_singleton
            - fail_hard_filters - (ht.QD < 2) | (ht.FS > 60) | (ht.MQ < 30)

    Numerical features are median-imputed. If impute_features_by_variant_type is set, imputation is done based on the
    median of the variant type.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param bool impute_features_by_variant_type: Whether to impute features median by variant type
    :param str group: Whether to use 'raw' or 'adj' genotypes
    :return: Hail Table ready for RF
    :rtype: Table
    """
    logger.info("Loading split hardcall MatrixTable and annotation Tables")
    ht = get_ukbb_data(data_source, freeze, meta_root="meta").rows()

    inbreeding_ht = hl.read_table(
        var_annotations_ht_path("inbreeding_coeff", data_source, freeze)
    )
    trio_stats_ht = hl.read_table(
        var_annotations_ht_path("trio_stats", data_source, freeze)
    )
    trio_stats_ht = trio_stats_ht.select(
        f"n_transmitted_{group}", f"ac_children_{group}"
    )
    sibling_stats_ht = hl.read_table(
        var_annotations_ht_path("sib_stats", data_source, freeze)
    )
    sibling_stats_ht = sibling_stats_ht.select(
        f"n_sib_shared_variants_{group}", f"ac_sibs_{group}"
    )
    truth_data_ht = hl.read_table(
        var_annotations_ht_path("truth_data", data_source, freeze)
    )
    truth_data_ht = truth_data_ht.select(*TRUTH_DATA)
    allele_data_ht = hl.read_table(
        var_annotations_ht_path("allele_data", data_source, freeze)
    )
    allele_counts_ht = hl.read_table(
        var_annotations_ht_path("allele_counts", data_source, freeze)
    )
    info_ht = hl.read_table(info_ht_path(data_source, freeze, split=True))
    info_ht = info_ht.transmute(**info_ht.info)
    info_ht = info_ht.select(
        "AS_lowqual", "interval_qc_pass", "FS", "MQ", "QD", *INFO_FEATURES
    )

    logger.info("Annotating MatrixTable with all columns from annotation Tables")
    ht = ht.annotate(
        **info_ht[ht.key],
        **inbreeding_ht[ht.key],
        **trio_stats_ht[ht.key],
        **sibling_stats_ht[ht.key],
        **truth_data_ht[ht.key],
        **allele_data_ht[ht.key].allele_data,
        **allele_counts_ht[ht.key],
    )
    # Filter to only variants found in high qual samples (included controls in count) and with no LowQual filter
    ht = ht.filter((ht[f"ac_qc_samples_{group}"] > 0) & ~ht.AS_lowqual)
    ht = ht.select(
        "a_index",
        "was_split",
        "interval_qc_pass",
        *FEATURES,
        *TRUTH_DATA,
        **{
            "transmitted_singleton": (ht[f"n_transmitted_{group}"] == 1)
            & (ht[f"ac_qc_samples_{group}"] == 2)
            & (ht[f"ac_children_{group}"] == 1),
            "sibling_singleton": (ht[f"n_sib_shared_variants_{group}"] == 1)
            & (ht[f"ac_qc_samples_{group}"] == 2)
            & (ht[f"ac_sibs_{group}"] == 2),
            "fail_hard_filters": (ht.QD < 2) | (ht.FS > 60) | (ht.MQ < 30),
        },
        singleton=ht.ac_release_samples_raw == 1,
        ac_raw=ht.ac_release_samples_raw,
        ac=ht.ac_release_samples_adj,
    )

    ht = ht.repartition(5000, shuffle=False)
    ht = ht.persist()

    if impute_features_by_variant_type:
        ht = median_impute_features(ht, {"variant_type": ht.variant_type})

    summary = ht.group_by(
        "omni",
        "mills",
        "transmitted_singleton",
        "ukbb_array_con_common",
        "sibling_singleton",
    ).aggregate(n=hl.agg.count())
    logger.info("Summary of truth data annotations:")
    summary.show(20)

    return ht


def get_rf_runs(rf_json_fp: str) -> Dict:
    """
    Loads RF run data from JSON file.

    :param rf_json_fp: File path to rf json file.
    :return: Dictionary containing the content of the JSON file, or an empty dictionary if the file wasn't found.
    """
    if file_exists(rf_json_fp):
        with hl.hadoop_open(rf_json_fp) as f:
            return json.load(f)
    else:
        logger.warning(
            f"File {rf_json_fp} could not be found. Returning empty RF run hash dict."
        )
        return {}


def get_run_data(
    data_source: str,
    freeze: int,
    interval_qc_filter: bool,
    transmitted_singletons: bool,
    sibling_singletons: bool,
    array_con_common: bool,
    adj: bool,
    vqsr_training: bool,
    test_intervals: List[str],
    features_importance: Dict[str, float],
    test_results: List[hl.tstruct],
) -> Dict:
    """
    Creates a Dict containing information about the RF input arguments and feature importance

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param bool interval_qc_filter: True if RF training variants were filtered to intervals with good coverage
    :param bool transmitted_singletons: True if transmitted singletons were used in training
    :param bool sibling_singletons: True if sibling singletons were used in training
    :param bool array_con_common: True if common variants with high concordance to the array data were used in training
    :param bool adj: True if training variants were filtered by adj
    :param List of str test_intervals: Intervals withheld from training to be used in testing
    :param Dict of float keyed by str features_importance: Feature importance returned by the RF
    :param List of struct test_results: Accuracy results from applying RF model to the test intervals
    :return: Dict of RF information
    """
    if vqsr_training:
        no_transmitted_singletons = no_sibling_singletons = no_array_con_common = None

    run_data = {
        "data": data_source,
        "freeze": freeze,
        "input_args": {
            "interval_qc_filter": interval_qc_filter,
            "transmitted_singletons": transmitted_singletons,
            "sibling_singletons": sibling_singletons,
            "array_con_common": array_con_common,
            "adj": adj,
            "vqsr_training": vqsr_training,
        },
        "features_importance": features_importance,
        "test_intervals": test_intervals,
    }

    if test_results is not None:
        tps = 0
        total = 0
        for row in test_results:
            values = list(row.values())
            # Note: values[0] is the TP/FP label and values[1] is the prediction
            if values[0] == values[1]:
                tps += values[2]
            total += values[2]
        run_data["test_results"] = [dict(x) for x in test_results]
        run_data["test_accuracy"] = tps / total

    return run_data


def generate_final_rf_ht(
    ht: hl.Table,
    ac0_filter_expr: hl.expr.BooleanExpression,
    ts_ac_filter_expr: hl.expr.BooleanExpression,
    snp_cutoff: Union[int, float],
    indel_cutoff: Union[int, float],
    determine_cutoff_from_bin: bool = False,
    aggregated_bin_ht: Optional[hl.Table] = None,
    inbreeding_coeff_cutoff: float = -0.3,
) -> hl.Table:
    """
    Prepares finalized RF model given an RF result table from `rf.apply_rf_model` and cutoffs for filtering.

    If `determine_cutoff_from_bin` is True, `aggregated_bin_ht` must be supplied to determine the SNP and indel RF
    probabilities to use as cutoffs from an aggregated quantile bin Table like one created by
    `compute_grouped_binned_ht` in combination with `score_bin_agg`.

    :param ht: RF result table from `rf.apply_rf_model` to prepare as the final RF Table
    :param ac0_filter_expr: Expression that indicates if a variant should be filtered as allele count 0 (AC0)
    :param ts_ac_filter_expr: Expression in `ht` that indicates if a variant is a transmitted singleton
    :param snp_cutoff: RF probability or bin (if `determine_cutoff_from_bin` True) to use for SNP variant QC filter
    :param indel_cutoff: RF probability or bin (if `determine_cutoff_from_bin` True) to use for indel variant QC filter
    :param determine_cutoff_from_bin: If True the RF probability will be determined using bin info in `aggregated_bin_ht`
    :param aggregated_bin_ht: File with aggregate counts of variants based on quantile bins
    :param inbreeding_coeff_cutoff: InbreedingCoeff hard filter to use for variants
    :return: Finalized random forest Table annotated with variant filters
    """
    # Determine SNP and indel RF cutoffs if given bin instead of RF probability
    if determine_cutoff_from_bin:
        snp_rf_cutoff, indel_rf_cutoff = aggregated_bin_ht.aggregate(
            [
                hl.agg.filter(
                    snv & (aggregated_bin_ht.bin == snp_cutoff),
                    hl.agg.min(aggregated_bin_ht.min_score),
                )
                for snv, cutoff in [
                    (aggregated_bin_ht.snv, snp_cutoff),
                    (~aggregated_bin_ht.snv, indel_cutoff),
                ]
            ]
        )
        snp_cutoff_global = hl.struct(bin=snp_cutoff, min_score=snp_rf_cutoff)
        indel_cutoff_global = hl.struct(bin=indel_cutoff, min_score=indel_rf_cutoff)

        logger.info(
            f"Using a SNP RF probability cutoff of {snp_rf_cutoff} and an indel RF probability cutoff of {indel_rf_cutoff}."
        )
    else:
        snp_cutoff_global = hl.struct(min_score=snp_cutoff)
        indel_cutoff_global = hl.struct(min_score=indel_cutoff)

    ht = ht.annotate_globals(
        rf_snv_cutoff=snp_cutoff_global, rf_indel_cutoff=indel_cutoff_global
    )

    # Add filters to RF HT
    filters = dict()

    if ht.any(hl.is_missing(ht.rf_probability["TP"])):
        raise ValueError("Missing RF probability!")

    filters["RF"] = (
        hl.is_snp(ht.alleles[0], ht.alleles[1])
        & (ht.rf_probability["TP"] < ht.rf_snv_cutoff.min_score)
    ) | (
        ~hl.is_snp(ht.alleles[0], ht.alleles[1])
        & (ht.rf_probability["TP"] < ht.rf_indel_cutoff.min_score)
    )

    filters["InbreedingCoeff"] = hl.or_else(
        ht.InbreedingCoeff < inbreeding_coeff_cutoff, False
    )
    filters["AC0"] = ac0_filter_expr

    ht = ht.annotate(filters=add_filters_expr(filters=filters))

    # Fix annotations for release
    annotations_expr = {
        "tp": hl.or_else(ht.tp, False),
        "transmitted_singleton": hl.or_missing(
            ts_ac_filter_expr, ht.transmitted_singleton
        ),
        "rf_probability": ht.rf_probability["TP"],
    }
    if "feature_imputed" in ht.row:
        annotations_expr.update(
            {
                x: hl.or_missing(~ht.feature_imputed[x], ht[x])
                for x in [f for f in ht.row.feature_imputed]
            }
        )

    ht = ht.transmute(**annotations_expr)

    return ht


def main(args):
    hl.init(log="/ukbb_variant_qc.log")

    data_source = "broad"
    freeze = args.freeze

    if args.list_rf_runs:
        logger.info(f"RF runs for {data_source}.freeze_{freeze}:")
        pretty_print_runs(get_rf_runs(rf_run_hash_path(data_source, freeze)))

    if args.annotate_for_rf:
        ht = create_rf_ht(
            data_source,
            freeze,
            impute_features_by_variant_type=not args.impute_features_no_variant_type,
            group="adj" if args.adj else "raw",
        )
        ht.write(
            rf_annotated_path(data_source, freeze, args.adj), overwrite=args.overwrite
        )
        logger.info(f"Completed annotation wrangling for random forests model training")

    if args.train_rf:
        run_hash = str(uuid.uuid4())[:8]
        rf_runs = get_rf_runs(rf_run_hash_path(data_source, freeze))
        while run_hash in rf_runs:
            run_hash = str(uuid.uuid4())[:8]

        ht = hl.read_table(rf_annotated_path(data_source, freeze, args.adj))

        if args.vqsr_training:
            vqsr_ht = hl.read_table(
                var_annotations_ht_path(
                    "vqsr" if args.vqsr_type == "AS" else "AS_TS_vqsr",
                    data_source,
                    freeze,
                )
            )
            ht = ht.annotate(
                vqsr_POSITIVE_TRAIN_SITE=vqsr_ht[ht.key].info.POSITIVE_TRAIN_SITE,
                vqsr_NEGATIVE_TRAIN_SITE=vqsr_ht[ht.key].info.NEGATIVE_TRAIN_SITE,
            )
            tp_expr = ht.vqsr_POSITIVE_TRAIN_SITE
            fp_expr = ht.vqsr_NEGATIVE_TRAIN_SITE
        else:
            fp_expr = ht.fail_hard_filters
            tp_expr = ht.omni | ht.mills
            if not args.no_array_con_common:
                tp_expr = tp_expr | ht.ukbb_array_con_common
            if not args.no_transmitted_singletons:
                tp_expr = tp_expr | ht.transmitted_singleton
            if not args.no_sibling_singletons:
                tp_expr = tp_expr | ht.sibling_singleton

            if isinstance(args.test_intervals, str):
                test_intervals = [args.test_intervals]
            if args.test_intervals:
                test_intervals = [hl.parse_locus_interval(x) for x in args.test_intervals]

        if args.interval_qc_filter:
            interval_filter_expr = ht.interval_qc_pass
        else:
            interval_filter_expr = True

        ht, rf_model = train_rf_model(
            ht.filter(interval_filter_expr),
            rf_features=FEATURES,
            tp_expr=tp_expr,
            fp_expr=fp_expr,
            fp_to_tp=args.fp_to_tp,
            num_trees=args.num_trees,
            max_depth=args.max_depth,
            test_expr=hl.literal(test_intervals).any(lambda interval: interval.contains(ht.locus)),
        )

        logger.info("Saving RF model")
        save_model(
            rf_model,
            rf_path(data_source, freeze, data="model", run_hash=run_hash),
            overwrite=args.overwrite,
        )

        ht = ht.checkpoint(
            rf_path(data_source, freeze, data="training", run_hash=run_hash),
            overwrite=args.overwrite,
        )

        logger.info("Adding run to RF run list")
        rf_runs[run_hash] = get_run_data(
            data_source,
            freeze,
            interval_qc_filter=args.interval_qc_filter,
            transmitted_singletons=not args.no_transmitted_singletons,
            sibling_singletons=not args.no_sibling_singletons,
            array_con_common=not args.no_array_con_common,
            adj=args.adj,
            vqsr_training=args.vqsr_training,
            test_intervals=hl.eval(ht.test_intervals),
            features_importance=hl.eval(ht.features_importance),
            test_results=hl.eval(ht.test_results),
        )

        with hl.hadoop_open(rf_run_hash_path(data_source, freeze), "w") as f:
            json.dump(rf_runs, f)

        logger.info(f"Completed training of random forests model")

    else:
        run_hash = args.run_hash

    if args.apply_rf:
        logger.info(f"Applying RF model {run_hash} to {data_source}.freeze_{freeze}.")

        rf_model = load_model(
            rf_path(data_source, freeze, data="model", run_hash=run_hash)
        )
        ht = hl.read_table(
            rf_path(data_source, freeze, data="training", run_hash=run_hash)
        )

        ht = apply_rf_model(ht, rf_model, FEATURES, label=LABEL_COL)

        logger.info("Finished applying RF model")
        ht = ht.annotate_globals(rf_hash=run_hash)
        ht = ht.checkpoint(
            rf_path(data_source, freeze, "rf_result", run_hash=run_hash),
            overwrite=args.overwrite,
        )

        ht_summary = ht.group_by("tp", TRAIN_COL, LABEL_COL, PREDICTION_COL).aggregate(
            n=hl.agg.count()
        )
        ht_summary.show(n=20)

    if args.finalize:
        ht = hl.read_table(
            rf_path(data_source, freeze, "rf_result", run_hash=args.run_hash)
        )
        freq_ht = hl.read_table(
            var_annotations_ht_path("ukb_freq", data_source, freeze)
        )

        aggregated_bin_ht = score_quantile_bin_path(
            args.run_hash, data_source, freeze, aggregated=True
        )
        if not file_exists(aggregated_bin_ht):
            sys.exit(
                f"Could not find binned HT for RF  run {args.run_hash} ({aggregated_bin_ht}). Please run create_ranked_scores.py for that hash."
            )
        ht = generate_final_rf_ht(
            ht,
            ac0_filter_expr=freq_ht[ht.key].freq[0].AC == 0,
            ts_ac_filter_expr=freq_ht[ht.key].freq[1].AC == 2,
            snp_cutoff=args.snp_cutoff,
            indel_cutoff=args.indel_cutoff,
            determine_cutoff_from_bin=args.treat_cutoff_as_prob,
            aggregated_bin_ht=aggregated_bin_ht,
            inbreeding_coeff_cutoff=INBREEDING_COEFF_HARD_CUTOFF,
        )
        # This column is added by the RF module based on a 0.5 threshold which doesn't correspond to what we use
        ht = ht.drop(ht[PREDICTION_COL])

        ht.describe()
        ht.show(20)
        ht.write(var_annotations_ht_path(data_source, freeze, "rf"), args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--run_hash",
        help="Run hash. Created by --train_rf and only needed for --apply_rf without running --train_rf",
        required=False,
    )
    parser.add_argument(
        "--snp_cutoff", help="Percentile to set RF cutoff", type=float, default=90.0
    )
    parser.add_argument(
        "--indel_cutoff", help="Percentile to set RF cutoff", type=float, default=80.0
    )
    parser.add_argument(
        "--treat_cutoff_as_prob",
        help="If set snp_cutoff and indel_cutoff will be probability rather than percentile ",
        action="store_true",
    )
    actions = parser.add_argument_group("Actions")
    actions.add_argument(
        "--list_rf_runs",
        help="Lists all previous RF runs, along with their hash,  parameters and testing results.",
        action="store_true",
    )
    actions.add_argument(
        "--annotate_for_rf",
        help="Creates an annotated ht with features for RF",
        action="store_true",
    )
    actions.add_argument("--train_rf", help="Trains RF model", action="store_true")
    actions.add_argument(
        "--apply_rf", help="Applies RF model to the data", action="store_true"
    )
    actions.add_argument("--finalize", help="Write final RF model", action="store_true")
    annotate_params = parser.add_argument_group("Annotate features params")
    annotate_params.add_argument(
        "--impute_features_no_variant_type",
        help="If set, feature imputation is NOT split by variant type.",
        action="store_true",
    )
    rf_params = parser.add_argument_group("Random Forest parameters")
    rf_params.add_argument(
        "--fp_to_tp",
        help="Ratio of FPs to TPs for creating the RF model. If set to 0, all training examples are used. (default=1.0)",
        default=1.0,
        type=float,
    )
    rf_params.add_argument(
        "--test_intervals",
        help='The specified interval(s) will be held out for testing and used for evaluation only. (default to "chr20")',
        nargs="+",
        type=str,
        default="chr20",
    )
    rf_params.add_argument(
        "--num_trees",
        help="Number of trees in the RF model. (default=500)",
        default=500,
        type=int,
    )
    rf_params.add_argument(
        "--max_depth",
        help="Maxmimum tree depth in the RF model. (default=5)",
        default=5,
        type=int,
    )
    training_params = parser.add_argument_group("Training data parameters")
    training_params.add_argument(
        "--adj", help="Use adj genotypes.", action="store_true"
    )
    training_params.add_argument(
        "--vqsr_training", help="Use VQSR training examples", action="store_true"
    )
    training_params.add_argument(
        "--vqsr_type",
        help="If a string is provided the VQSR training annotations will be used for training.",
        type=str,
    )
    training_params.add_argument(
        "--no_transmitted_singletons",
        help="Do not use transmitted singletons for training.",
        action="store_true",
    )
    training_params.add_argument(
        "--no_sibling_singletons",
        help="Do not use sibling singletons for training.",
        action="store_true",
    )
    training_params.add_argument(
        "--no_array_con_common",
        help="Do not use common concordant array variants.",
        action="store_true",
    )
    training_params.add_argument(
        "--interval_qc_filter",
        help="Should interval QC be applied before RF training.",
        action="store_true",
    )

    args = parser.parse_args()

    if not args.run_hash and not args.train_rf and args.apply_rf:
        sys.exit(
            "Error: --run_hash is required when running --apply_rf without running --train_rf too."
        )

    if args.run_hash and args.train_rf:
        sys.exit(
            "Error: --run_hash and --train_rf are mutually exclusive. --train_rf will generate a run hash."
        )

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
