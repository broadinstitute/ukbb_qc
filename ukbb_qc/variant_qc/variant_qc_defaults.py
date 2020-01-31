

LABEL_COL = 'rf_label'
TRAIN_COL = 'rf_train'
prediction_col_name = 'rf_prediction'

SITES_FEATURES = [
    'info_MQRankSum',
    'info_SOR',
    'inbreeding_coeff',
    'info_ReadPosRankSum'
]

VQSR_FEATURES = [
    'info_FS',
    'info_QD',
    'info_MQ',
    'info_VarDP'
]

ALLELE_FEATURES = [
    'variant_type',
    'allele_type',
    'n_alt_alleles',
    'was_mixed',
    'has_star',
    'qd',
    'pab_max'
]

MEDIAN_FEATURES = [
    'gq_median',
    'dp_median',
    'nrq_median',
    'ab_median'
]

INBREEDING_COEFF_HARD_CUTOFF = -0.3


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
        ht.inbreeding_coeff,
        ht.info.ReadPosRankSum
    ]))


def get_vqsr_features_expr(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
    """
    Returns expressions to annotate site-level features only used by VQSRs

    :param Table ht: Table to create annotation expression for.
    :return: Dict with keys containing column names and values expressions
    :rtype: Dict of str: Expression
    """
    return dict(zip(VQSR_FEATURES, [
        ht.info.FS,
        ht.info.QD,
        ht.info.MQ,
        ht.info.VarDP
    ]))


def get_allele_features_expr(ht: hl.Table, qc_stats_group_index: int) -> Dict[str, hl.expr.Expression]:
    """
    Returns expressions to annotate allele-level features (RF only)

    :param Table ht: Table to create annotation expression for.
    :param int qc_stats_group_index: Index of group to get stats for
    :return: Dict with keys containing column names and values expressions
    :rtype: Dict of str: Expression
    """
    return dict(zip(ALLELE_FEATURES, [
        ht.allele_data.variant_type,
        ht.allele_data.allele_type,
        ht.allele_data.n_alt_alleles,
        ht.allele_data.was_mixed,
        ht.allele_data.has_star,
        ht.qc_stats[qc_stats_group_index].qd,
        ht.qc_stats[qc_stats_group_index].pab.max
    ]))


def get_median_features_expr(ht: hl.Table, qc_stats_group_index: int) -> Dict[str, hl.expr.Expression]:
    """
    Returns expressions to annotate 2.0.2 allele-specific features (RF only)

    :param Table ht: Table to create annotation expression for.
    :param int qc_stats_group_index: Index of group to get stats for
    :return: Dict with keys containing column names and values expressions
    :rtype: Dict of str: Expression
    """
    return dict(zip(MEDIAN_FEATURES, [
        ht.qc_stats[qc_stats_group_index].gq_median,
        ht.qc_stats[qc_stats_group_index].dp_median,
        ht.qc_stats[qc_stats_group_index].nrq_median,
        ht.qc_stats[qc_stats_group_index].ab_median
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
        'fail_hard_filters': (ht.info.QD < 2) | (ht.info.FS > 60) | (ht.info.MQ < 30)
        'info_POSITIVE_TRAIN_SITE': ht.info.POSITIVE_TRAIN_SITE,
        'info_NEGATIVE_TRAIN_SITE': ht.info.NEGATIVE_TRAIN_SITE
    }
