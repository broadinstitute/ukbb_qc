from gnomad_hail import *
from ukbb_qc.resources import *
from ukbb_qc.variant_qc.variant_qc_defaults import *

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variant_qc_functions")
logger.setLevel(logging.INFO)


def create_rf_ht(
        n_variants_median: int = 50000,
        impute_features_by_variant_type: bool = True,
        group: str = 'qc_samples_raw'
) -> hl.Table:
    """
    Creates a Table with all necessary columns for RF:
    - Features columns (both VQSR and RF)
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

    mt = get_ukbb_data(data_source, freeze, meta_root='meta')
    ht_family_stats = hl.read_table(var_annotations_ht_path(data_source, freeze, 'family_stats'))
    ht_qc_stats = hl.read_table(var_annotations_ht_path(data_source, freeze, 'qc_stats'))
    ht_call_stats = hl.read_table(var_annotations_ht_path(data_source, freeze, 'call_stats'))
    ht_truth_data = hl.read_table(var_annotations_ht_path(data_source, freeze, 'truth_data'))
    ht_allele_data = hl.read_table(var_annotations_ht_path(data_source, freeze, 'allele_data'))

    # TODO: remove hard path for info HT
    info_ht = hl.read_table(var_annotations_ht_path(data_source, freeze, "vqsr" if vqsr_type == "AS" else "AS_TS_vqsr"))
    mt = mt.annotate_rows(info=info_ht[mt.row_key].info)
    mt = mt.annotate_rows(info=mt.info.annotate(AC=hl.agg.sum(mt.GT.n_alt_alleles())))

    # Add inbreeding coefficient
    ic_mt = mt.filter_cols(mt.meta.related_filter, keep=False)
    ic_mt = ic_mt.filter_rows(hl.agg.sum(ic_mt.GT.n_alt_alleles()) > 0)
    ic_mt = ic_mt.annotate_rows(inbreeding_coeff=bi_allelic_site_inbreeding_expr(ic_mt.GT))
    ic_ht = ic_mt.rows().select('inbreeding_coeff')
    mt = mt.annotate_rows(**ic_ht[mt.row_key])

    ht = mt.annotate_rows(
        n_nonref=ht_call_stats[mt.row_key].qc_callstats[0].AC[1],
        singleton=mt.info.AC == 1,#[mt.a_index - 1] == 1,
        info_ac=mt.info.AC, #[mt.a_index - 1],
        **ht_call_stats[mt.row_key],
        **ht_family_stats[mt.row_key],
        **ht_qc_stats[mt.row_key],
        **ht_truth_data[mt.row_key],
        **ht_allele_data[mt.row_key]
    ).rows()

    ht = ht.annotate(AC=ht.qc_callstats.find(lambda x: x.meta.get('group') == group).AC[1])

    # Filter to only variants found in high qual samples and with no LowQual filter
    ht = ht.filter((ht.qc_callstats.find(lambda x: x.meta.get('group') == group).AC[1] > 0) & ~info_ht[ht.key].filters.contains("LowQual"))

    qc_stats_groups = ht_qc_stats.aggregate(hl.agg.take(ht_qc_stats.qc_stats.map(lambda x: x.meta['group']), 1))[0]
    qc_stats_group_index = next(x[0] for x in enumerate(qc_stats_groups) if x[1] == group)
    family_stats_groups = ht_family_stats.aggregate(hl.agg.take(ht_family_stats.family_stats.map(lambda x: x.meta['group']), 1))[0]
    family_stats_group_index = next(x[0] for x in enumerate(family_stats_groups) if x[1] == 'raw')

    ht = ht.select(
        **get_allele_features_expr(ht, qc_stats_group_index),
        **get_site_features_expr(ht),
        **get_vqsr_features_expr(ht),
        **get_training_sites_expr(ht, family_stats_group_index),
        **get_median_features_expr(ht, qc_stats_group_index),
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

