import json
import logging
import argparse
import uuid
import sys
from typing import Dict, List
from gnomad_hail.utils import rf
from gnomad_hail.utils.generic import bi_allelic_site_inbreeding_expr
from gnomad_hail.utils.gnomad_functions import pretty_print_runs
from gnomad_hail.utils.slack import try_slack
from gnomad_qc.v2.variant_qc.variantqc import sample_rf_training_examples
from ukbb_qc.utils import annotate_interval_qc_filter
from ukbb_qc.resources.basics import get_ukbb_data, CURRENT_FREEZE
from ukbb_qc.resources.variant_qc import (rf_run_hash_path, var_annotations_ht_path, score_ranking_path,
                                          rf_path, rf_annotated_path)
import hail as hl

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("run_ukbb_variant_qc")
logger.setLevel(logging.INFO)


LABEL_COL = 'rf_label'
TRAIN_COL = 'rf_train'
prediction_col_name = 'rf_prediction'

SITES_FEATURES = [
    'info_MQRankSum',
    'info_SOR',
    'InbreedingCoeff',
    'info_ReadPosRankSum'
]

ALLELE_FEATURES = [
    'variant_type',
    'allele_type',
    'n_alt_alleles',
    'was_mixed',
    'has_star'
    'qd',
    'pab_max'
]

INBREEDING_COEFF_HARD_CUTOFF = -0.3


def get_rf_runs(rf_json_fp: str) -> Dict:
    """
    Loads RF run data from JSON file.

    :param rf_json_fp: file path to rf json file
    :return: Dictionary containing the content of the JSON file, or an empty dictionary if the file wasn't found.
    """
    if hl.utils.hadoop_exists(rf_json_fp):
        with hl.hadoop_open() as f:
            return json.load(f)
    else:
        logger.warning(f"File {rf_json_fp} could not be found. Returning empty RF run hash dict.")
        return {}


def get_features_list(
        sites_features: bool,
        allele_features: bool
) -> List[str]:
    """
    Returns the list of features to use based on desired arguments (currently only VQSR / alleles)

    :param bool sites_features: Whether to use site-level features
    :param bool allele_features: Whether to use Allele-specific features
    :return: List of features to use
    """

    features = []
    if sites_features:
        features.extend(SITES_FEATURES)
    if allele_features:
        features.extend(ALLELE_FEATURES)

    return features


def create_rf_ht(
        data_source: str,
        freeze: int,
        n_variants_median: int = 50000,
        impute_features_by_variant_type: bool = True,
        group: str = 'qc_samples_raw',
        vqsr_type = 'AS'
) -> hl.Table:
    """
    Creates a Table with all necessary columns for RF:
    - Features columns
    - Training criteria columns

    Numerical features are median-imputed. If impute_features_by_variant_type is set, imputation is done based on the median of the variant type.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param int n_variants_median: Number of variants to use for median computation
    :param bool impute_features_by_variant_type: Whether to impute features median by variant type
    :param str group: Whether to use 'raw' or 'adj' genotypes
    :return: Hail Table ready for RF
    :rtype: Table
    """

    def get_site_features_expr(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
        """
        Returns expressions to annotate site-level features used by both VQSR and RF.

        :param Table ht: Table to create annotation expression for.
        :return: Dict with keys containing column names and values expressions
        :rtype: Dict of str: Expression
        """
        return dict(zip(SITES_FEATURES, [
            ht.info.MQRankSum,
            ht.info.SOR,
            ht.info.InbreedingCoeff,
            ht.info.ReadPosRankSum
        ]))

    def get_allele_features_expr(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
        """
        Returns expressions to annotate allele-level features (RF only)

        :param Table ht: Table to create annotation expression for.
        :return: Dict with keys containing column names and values expressions
        :rtype: Dict of str: Expression
        """
        return dict(zip(ALLELE_FEATURES, [
            ht.allele_data.variant_type,
            ht.allele_data.allele_type,
            ht.allele_data.n_alt_alleles,
            ht.allele_data.was_mixed,
            ht.allele_data.has_star,
            ht.info.AS_QD,
            ht.info.AS_pab_max
        ]))

    def get_training_sites_expr(ht: hl.Table, family_stats_group_index: int) -> Dict[str, hl.expr.Expression]:
        """
        Returns expressions to columns to select training examples

        :param Table ht: Table to create annotation expression for.
        :param int family_stats_group_index: Index of group to get stats for
        :return: Dict with keys containing column names and values expressions
        :rtype: Dict of str: Expression
        """
        return {
            'transmitted_singleton': (ht.family_stats[family_stats_group_index].tdt.t == 1) &
                                     (ht.family_stats[family_stats_group_index].unrelated_qc_callstats.AC[1] == 1),
            'fail_hard_filters': (ht.info.QD < 2) | (ht.info.FS > 60) | (ht.info.MQ < 30),
            'info_POSITIVE_TRAIN_SITE': ht.info.POSITIVE_TRAIN_SITE,
            'info_NEGATIVE_TRAIN_SITE': ht.info.NEGATIVE_TRAIN_SITE
        }

    mt = get_ukbb_data(data_source, freeze, meta_root='meta')
    ht_family_stats = hl.read_table(var_annotations_ht_path(data_source, freeze, 'family_stats'))
    ht_call_stats = hl.read_table(var_annotations_ht_path(data_source, freeze, 'call_stats'))
    ht_truth_data = hl.read_table(var_annotations_ht_path(data_source, freeze, 'truth_data'))
    ht_allele_data = hl.read_table(var_annotations_ht_path(data_source, freeze, 'allele_data'))

    # TODO: remove hard path for info HT
    info_ht = hl.read_table(var_annotations_ht_path(data_source, freeze, "vqsr" if vqsr_type == "AS" else "AS_TS_vqsr"))
    mt = mt.annotate_rows(info=info_ht[mt.row_key].info)
    mt = mt.annotate_rows(info=mt.info.annotate(AC=hl.agg.sum(mt.GT.n_alt_alleles())))

    ht = mt.annotate_rows(
        n_nonref=ht_call_stats[mt.row_key].qc_callstats[0].AC[1],
        singleton=mt.info.AC == 1,#[mt.a_index - 1] == 1,
        info_ac=mt.info.AC, #[mt.a_index - 1],
        **ht_call_stats[mt.row_key],
        **ht_family_stats[mt.row_key],
        **ht_truth_data[mt.row_key],
        **ht_allele_data[mt.row_key]
    ).rows()

    ht = ht.annotate(AC=ht.qc_callstats.find(lambda x: x.meta.get('group') == group).AC[1])

    # Filter to only variants found in high qual samples and with no LowQual filter
    ht = ht.filter((ht.qc_callstats.find(lambda x: x.meta.get('group') == group).AC[1] > 0) & ~info_ht[ht.key].filters.contains("LowQual"))

    #What is this doing? Is it just getting the index for each of the groups in qc_stats and family stats, I think this all changes
    family_stats_groups = ht_family_stats.aggregate(hl.agg.take(ht_family_stats.family_stats.map(lambda x: x.meta['group']), 1))[0]
    family_stats_group_index = next(x[0] for x in enumerate(family_stats_groups) if x[1] == 'raw')

    ht = ht.select(
        **get_allele_features_expr(ht, qc_stats_group_index),
        **get_site_features_expr(ht),
        **get_training_sites_expr(ht, family_stats_group_index),
        omni=ht.truth_data.omni,
        mills=ht.truth_data.mills,
        #ukbb_array=ht.truth_data.ukbb_array,
        #ukbb_array_con=ht.truth_data.ukbb_array_con,
        ukbb_array_con_common=ht.truth_data.ukbb_array_con_common,
        sib_singletons=ht.truth_data.sib_singletons,
        n_nonref=ht.n_nonref,
        singleton=ht.singleton,
        was_split=ht.was_split,
        info_ac=ht.info_ac
    )

    # Annotate variants by type (e.g. for downsampling purpose)
    variants_by_type = ht.aggregate(hl.agg.counter(ht.variant_type))
    logger.info('Variants by type:\n{}'.format('\n'.join(['{}: {}'.format(k, v) for k, v in variants_by_type.items()])))

    ht = ht.annotate_globals(variants_by_type=variants_by_type)

    #ht = ht.persist()
    ht = ht.checkpoint(f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/rf_annotations_cp.ht', overwrite=True)

    ht = ht.repartition(5000, shuffle=False)

    # Compute medians
    # numerical_features = [k for k, v in ht.row.dtype.items() if annotation_type_is_numeric(v)]
    numerical_features = [k for k, v in ht.row.dtype.items() if v == hl.tint or v == hl.tfloat]

    if impute_features_by_variant_type:
        prob_sample = hl.literal({v: min([n, n_variants_median])/n for v, n in variants_by_type.items()})
        ht_by_variant = ht.group_by(ht.variant_type).partition_hint(1)
        medians = ht_by_variant.aggregate(
            **{feature: hl.median(
                hl.agg.filter(
                    hl.rand_bool(prob_sample[ht.variant_type]),
                    hl.agg.collect(ht[feature])
                )
            ) for feature in numerical_features}
        ).collect()

        ht = ht.annotate_globals(feature_medians={x.variant_type: x for x in medians})
        ht = ht.annotate(
            **{f: hl.or_else(ht[f], ht.globals.feature_medians[ht.variant_type][f]) for f in numerical_features},
            feature_imputed=hl.struct(
                **{f: hl.is_missing(ht[f]) for f in numerical_features}
            )
        )
    else:
        prob_sample = min([sum(variants_by_type.values()), n_variants_median]) / sum(variants_by_type.values())
        medians = ht.aggregate(
            hl.struct(
                **{feature: hl.median(
                    hl.agg.filter(
                        hl.rand_bool(prob_sample),
                        hl.agg.collect(ht[feature])
                    )
                ) for feature in numerical_features}
            )
        )
        ht = ht.annotate_globals(features_median=medians)
        ht = ht.annotate(
            **{f: hl.or_else(ht[f], ht.globals.features_median[f]) for f in numerical_features},
            feature_imputed=hl.struct(
                **{f: hl.is_missing(ht[f]) for f in numerical_features}
            )
        )

    return ht


def get_run_data(
        data_source: str,
        freeze: int,
        interval_qc_filter,
        vqsr_training: bool,
        no_transmitted_singletons: bool,
        array_con: bool,
        array_con_common: bool,
        adj: bool,
        test_intervals: List[str],
        features_importance: Dict[str, float],
        test_results: List[hl.tstruct] = None
) -> Dict:
    run_data = {
        'data': data_source,
        'freeze': freeze,
        'input_args': {
            'interval_qc_filter':interval_qc_filter,
            'vqsr_training': vqsr_training,
            'no_transmitted_singletons': no_transmitted_singletons,
            'array_con': array_con,
            'array_con_common': array_con_common,
            'adj': adj
        },
        'features_importance': features_importance,
        'test_intervals': test_intervals
    }

    if test_results:
        tps = 0
        total = 0
        for row in test_results:
            values = list(row.values())
            if values[0] == values[1]:
                tps += values[2]
            total += values[2]
        run_data['test_results'] = [dict(x) for x in test_results]
        run_data['test_accuracy'] = tps / total

    return run_data


def train_rf(data_source: str, freeze: int, args):

    # Get unique hash for run and load previous runs
    run_hash = str(uuid.uuid4())[:8]
    rf_runs = get_rf_runs(rf_run_hash_path(data_source, freeze))
    while run_hash in rf_runs:
        run_hash = str(uuid.uuid4())[:8]

    ht = hl.read_table(rf_annotated_path(data_source, freeze, args.adj))


    #interval_qc_ht = hl.read_table(interval_qc_path(data_source, freeze))
    #good_intervals_ht = interval_qc_ht.filter(interval_qc_ht.pct_samples_20x > args.pct_samples_20x).key_by('interval')
    #ht = ht.annotate(interval_qc_pass=hl.is_defined(good_intervals_ht[ht.locus]))
    ht = annotate_interval_qc_filter('broad', 5, ht)


    # TODO: Revisit this
    #summary = ht.group_by('omni', 'mills', 'transmitted_singleton', 'ukbb_array', 'ukbb_array_con', 'ukbb_array_con_common', 'sib_singletons').aggregate(n=hl.agg.count())
    summary = ht.group_by('omni', 'mills', 'transmitted_singleton', 'ukbb_array_con_common', 'sib_singletons').aggregate(n=hl.agg.count())
    logger.info("Summary of truth data annotations:")
    summary.show(20)

    # ht = ht.repartition(500, shuffle=False)

    if not args.vqsr_training:
        tp_expr = ht.omni | ht.mills
        if args.array_con_common:
            tp_expr = tp_expr | ht.ukbb_array_con_common
        # TODO: Revisit this
        #elif args.array_con:
        #    tp_expr = tp_expr | ht.ukbb_array_con
        #elif args.array:
        #    tp_expr = tp_expr | ht.ukbb_array

        if not args.no_transmitted_singletons:
            tp_expr = tp_expr | ht.transmitted_singleton  # | ht.info_POSITIVE_TRAIN_SITE

        if not args.no_sibling_singletons:
            tp_expr = tp_expr | ht.sib_singletons  # | ht.info_POSITIVE_TRAIN_SITE

        ht = ht.annotate(
            tp=tp_expr
        )

    test_intervals_str = [] if not args.test_intervals else [args.test_intervals] if isinstance(args.test_intervals, str) else args.test_intervals
    test_intervals_locus = [hl.parse_locus_interval(x, reference_genome='GRCh38') for x in test_intervals_str]

    if test_intervals_locus:
        ht = ht.annotate_globals(
            test_intervals=test_intervals_locus
        )

    ht = sample_rf_training_examples(ht,
                                     tp_col='info_POSITIVE_TRAIN_SITE' if args.vqsr_training else 'tp',  # TODO: remove VQSR option
                                     fp_col='info_NEGATIVE_TRAIN_SITE' if args.vqsr_training else 'fail_hard_filters',  # TODO: remove VQSR option
                                     fp_to_tp=args.fp_to_tp)
    ht = ht.persist()

    summary = ht.group_by('tp', 'fail_hard_filters', 'rf_train', 'rf_label').aggregate(n=hl.agg.count())
    logger.info("Summary of TP/FP and RF training labels:")
    summary.show(n=20)

    rf_ht = ht.filter(ht[TRAIN_COL])
    if args.interval_qc_filter:
        rf_ht = rf_ht.filter(rf_ht.interval_qc_pass)

    rf_features = get_features_list(True,
                                    not (args.vqsr_features or args.median_features),
                                    args.vqsr_features,
                                    args.median_features)

    logger.info("Training RF model:\nfeatures: {}\nnum_tree: {}\nmax_depth:{}\nTest intervals: {}".format(
        ",".join(rf_features),
        args.num_trees,
        args.max_depth,
        ",".join(test_intervals_str)))

    rf_model = rf.train_rf(rf_ht,
                           features=rf_features,
                           label=LABEL_COL,
                           num_trees=args.num_trees,
                           max_depth=args.max_depth)

    logger.info("Saving RF model")
    rf.save_model(rf_model,
                  rf_path(data_source, freeze, data='model', run_hash=run_hash),
                  overwrite=args.overwrite
                  )

    test_results = None
    if args.test_intervals:
        logger.info("Testing model {} on intervals {}".format(run_hash, ",".join(test_intervals_str)))
        test_ht = hl.filter_intervals(ht, test_intervals_locus, keep=True)
        test_ht = test_ht.checkpoint(f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/test_rf.ht', overwrite=True)  # TODO: revisit temp placement?
        test_ht = test_ht.filter(hl.is_defined(test_ht[LABEL_COL]))
        test_results = rf.test_model(test_ht,
                                     rf_model,
                                     features=get_features_list(True, not (args.vqsr_features or args.median_features), args.vqsr_features, args.median_features),
                                     label=LABEL_COL)

        ht = ht.annotate_globals(
            test_results=test_results
        )

    logger.info("Writing RF training HT")
    features_importance = rf.get_features_importance(rf_model)
    ht = ht.annotate_globals(
        features_importance=features_importance,
        features=get_features_list(True, not (args.vqsr_features or args.median_features), args.vqsr_features, args.median_features),
        vqsr_training=args.vqsr_training,
        no_transmitted_singletons=args.no_transmitted_singletons,
        adj=args.adj
    )
    ht.write(rf_path(data_source, freeze, data='training', run_hash=run_hash), overwrite=args.overwrite)

    logger.info("Adding run to RF run list")
    rf_runs[run_hash] = get_run_data(
        data_source,
        freeze,
        args.interval_qc_filter,
        args.vqsr_training,
        args.no_transmitted_singletons,
        args.array_con,
        args.array_con_common,
        args.adj,
        test_intervals_str,
        features_importance,
        test_results
    )
    with hl.hadoop_open(rf_run_hash_path(data_source, freeze), 'w') as f:
        json.dump(rf_runs, f)

    return run_hash


def prepare_final_ht(data_source: str, freeze: int, run_hash: str, snp_cutoff: int, indel_cutoff: int, treat_cutoff_as_prob: bool) -> hl.Table:

    # Get snv and indel RF cutoffs based on bin
    binned_ht_path = score_ranking_path(data_source, freeze, run_hash, binned=True)
    if not hl.hadoop_exists(score_ranking_path(data_source, freeze, run_hash, binned=True)):
        sys.exit(f"Could not find binned HT for RF  run {run_hash} ({binned_ht_path}). Please run create_ranked_scores.py for that hash.")
    binned_ht = hl.read_table(binned_ht_path)
    if treat_cutoff_as_prob:
        snp_cutoff_global = hl.struct(min_score=snp_cutoff)
        indel_cutoff_global = hl.struct(min_score=indel_cutoff)
    else:
        snp_rf_cutoff, indel_rf_cutoff = binned_ht.aggregate([hl.agg.filter(binned_ht.snv & (binned_ht.bin == snp_cutoff), hl.agg.min(binned_ht.min_score)),
                                                          hl.agg.filter(~binned_ht.snv & (binned_ht.bin == indel_cutoff), hl.agg.min(binned_ht.min_score))])
        snp_cutoff_global = hl.struct(bin=snp_cutoff, min_score=snp_rf_cutoff)
        indel_cutoff_global = hl.struct(bin=indel_cutoff, min_score=indel_rf_cutoff)

    # Add filters to RF HT
    ht = hl.read_table(rf_path(data_source, freeze, 'rf_result', run_hash=run_hash))
    ht = ht.annotate_globals(rf_hash=run_hash,
                             rf_snv_cutoff=snp_cutoff_global,
                             rf_indel_cutoff=indel_cutoff_global)
    rf_filter_criteria = (hl.is_snp(ht.alleles[0], ht.alleles[1]) & (ht.rf_probability['TP'] < ht.rf_snv_cutoff.min_score)) | (
            ~hl.is_snp(ht.alleles[0], ht.alleles[1]) & (ht.rf_probability['TP'] < ht.rf_indel_cutoff.min_score))
    ht = ht.annotate(filters=hl.case()
                     .when(rf_filter_criteria, {'RF'})
                     .when(~rf_filter_criteria, hl.empty_set(hl.tstr))
                     .or_error('Missing RF probability!'))

    inbreeding_coeff_filter_criteria = hl.is_defined(ht.inbreeding_coeff) & (
            ht.inbreeding_coeff < INBREEDING_COEFF_HARD_CUTOFF)
    ht = ht.annotate(filters=hl.cond(inbreeding_coeff_filter_criteria,
                                     ht.filters.add('InbreedingCoeff'), ht.filters))

    freq_ht = hl.read_table(var_annotations_ht_path(data_source, freeze, 'join_freq'))
    ac0_filter_criteria = freq_ht[ht.key].freq[0].AC == 0

    ht = ht.annotate(
        filters=hl.cond(ac0_filter_criteria, ht.filters.add('AC0'), ht.filters)
    )

    # Fix annotations for release
    annotations_expr = {
        'tp': hl.or_else(ht.tp, False),
        'transmitted_singleton': hl.or_missing(freq_ht[ht.key].freq[1].AC == 2, ht.transmitted_singleton),
        'rf_probability': ht.rf_probability["TP"]
    }
    if 'feature_imputed' in ht.row:
        annotations_expr.update(
            {x: hl.or_missing(~ht.feature_imputed[x], ht[x]) for x in [f for f in ht.row.feature_imputed]}
        )

    ht = ht.transmute(
        **annotations_expr
    )

    #This column is added by the RF module based on a 0.5 threshold which doesn't correspond to what we use
    ht = ht.drop(ht.rf_prediction)

    return ht


def main(args):
    hl.init(log='/ukbb_variant_qc.log')

    data_source = args.data_source
    freeze = args.freeze

    if args.debug:
        logger.setLevel(logging.DEBUG)

    if args.list_rf_runs:
        logger.info(f"RF runs for {data_source}.freeze_{freeze}:")
        pretty_print_runs(get_rf_runs(rf_run_hash_path(data_source, freeze)))

    if args.annotate_for_rf:
        ht = create_rf_ht(data_source, freeze,
                          n_variants_median=args.n_variants_median,
                          impute_features_by_variant_type=not args.impute_features_no_variant_type,
                          group='qc_samples_raw')  # group='adj' if args.adj else 'raw'
        ht.write(rf_annotated_path(data_source, freeze, args.adj), overwrite=args.overwrite)
        logger.info(f"Completed annotation wrangling for random forests model training")

    run_hash = train_rf(data_source, freeze, args) if args.train_rf else args.run_hash
    logger.info(f"Completed training of random forests model")

    if args.apply_rf:
        logger.info(f"Applying RF model {run_hash} to {data_source}.freeze_{freeze}.")

        rf_model = rf.load_model(rf_path(data_source, freeze, data='model', run_hash=run_hash))
        ht = hl.read_table(rf_path(data_source, freeze, data='training', run_hash=run_hash))  #.persist()

        ht = rf.apply_rf_model(ht, rf_model, get_features_list(True, not (args.vqsr_features or args.median_features), args.vqsr_features, args.median_features), label=LABEL_COL)

        logger.info('Finished applying RF model')
        ht = ht.checkpoint(rf_path(data_source, freeze, 'rf_result', run_hash=run_hash), overwrite=args.overwrite)

        ht_summary = ht.group_by('tp', 'rf_train', 'rf_label', 'rf_prediction').aggregate(n=hl.agg.count())
        ht_summary.show(n=20)


    if args.finalize:
        ht = prepare_final_ht(data_source, freeze, args.run_hash, args.snp_cutoff, args.indel_cutoff, args.treat_cutoff_as_prob)
        ht = annotate_interval_qc_filter(data_source, freeze, ht, pct_samples=args.pct_samples_20x)
        ht.describe()
        ht.show(20)
        ht.write(var_annotations_ht_path(data_source, freeze, 'rf'), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)

    parser.add_argument('--run_hash', help='Run hash. Created by --train_rf and only needed for --apply_rf without running --train_rf', required=False)
    parser.add_argument('--snp_cutoff', help='Percentile to set RF cutoff', type=float, default=90.0)
    parser.add_argument('--indel_cutoff', help='Percentile to set RF cutoff', type=float, default=80.0)
    parser.add_argument('--treat_cutoff_as_prob', help='If set snp_cutoff and indel_cutoff will be probability rather than percentile ', action='store_true')

    parser.add_argument('--debug', help='Prints debug statements', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    actions = parser.add_argument_group('Actions')
    actions.add_argument('--list_rf_runs', help='Lists all previous RF runs, along with their hash,  parameters and testing results.', action='store_true')
    actions.add_argument('--annotate_for_rf', help='Creates an annotated ht with features for RF', action='store_true')
    actions.add_argument('--train_rf', help='Trains RF model', action='store_true')
    actions.add_argument('--apply_rf', help='Applies RF model to the data', action='store_true')
    actions.add_argument('--finalize', help='Write final RF model', action='store_true')

    annotate_params = parser.add_argument_group('Annotate features params')
    annotate_params.add_argument('--n_variants_median', help='Number of variants to use for median computation.', type=int, default=20000)
    annotate_params.add_argument('--impute_features_no_variant_type', help='If set, feature imputation is NOT split by variant type.', action='store_true')

    rf_params = parser.add_argument_group('Random Forest parameters')
    rf_params.add_argument('--fp_to_tp', help='Ratio of FPs to TPs for creating the RF model. If set to 0, all training examples are used. (default=1.0)', default=1.0,
                           type=float)
    rf_params.add_argument('--test_intervals', help='The specified interval(s) will be held out for testing and used for evaluation only. (default to "chr20")', nargs='+', type=str, default='chr20')
    rf_params.add_argument('--num_trees', help='Number of trees in the RF model. (default=500)', default=500, type=int)
    rf_params.add_argument('--max_depth', help='Maxmimum tree depth in the RF model. (default=5)', default=5, type=int)

    training_params = parser.add_argument_group('Training data parameters')
    training_params.add_argument('--adj', help='Use adj genotypes.', action='store_true')
    training_params.add_argument('--vqsr_features', help='Use VQSR features only (+ snv/indel variant type)',
                                 action='store_true')
    training_params.add_argument('--median_features', help='Use gnomAD 2.0.2 features',
                                 action='store_true')
    training_params.add_argument('--vqsr_training', help='Use VQSR training examples', action='store_true')
    training_params.add_argument('--no_transmitted_singletons', help='Do not use transmitted singletons for training.',
                                 action='store_true')
    training_params.add_argument('--no_sibling_singletons', help='Do not use sibling singletons for training.',
                                 action='store_true')
    training_params.add_argument('--array_con_common', help='Use only array variants >0.1% AF.',
                                 action='store_true')
    training_params.add_argument('--array_con', help='Use only highly concordant array variants.',
                                 action='store_true')
    training_params.add_argument('--array', help='Use all array variants.',
                                 action='store_true')
    training_params.add_argument('--interval_qc_filter', help='Should interval QC be applied', action='store_true')
    training_params.add_argument('--pct_samples_20x',
                        help='Percent samples at 20X to filter intervals',
                        default=0.85,
                        type=float)


    args = parser.parse_args()

    if not args.run_hash and not args.train_rf and args.apply_rf:
        sys.exit('Error: --run_hash is required when running --apply_rf without running --train_rf too.')

    if args.run_hash and args.train_rf:
        sys.exit('Error: --run_hash and --train_rf are mutually exclusive. --train_rf will generate a run hash.')

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
