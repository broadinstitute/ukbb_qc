from .basics import DataException

# Sample QC files
def sample_qc_prefix(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    if data_source not in DATA_SOURCES:
        raise DataException("This data_source is currently not present")
    if freeze not in FREEZES:
        raise DataException("This freeze is currently not present")
    return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/sample_qc'


def sample_list_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path of sample list for sample check
    :param str data_source: Will be regeneron or broad
    :param int freeze: One of data freezes
    :return: Path to sample list file
    :rtype: str
    """
    return f'{sample_qc_prefix(data_source, freeze)}/samples.list'


def meta_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path of meta ht for plotting
    :param str data_source: Will be regeneron or broad
    :param int freeze: One of data freezes
    :return: Path to meta ht
    :rtype: str
    """
    return f'{sample_qc_prefix(data_source, freeze)}/meta.ht'


def hard_filters_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/hard_filters_flagged.ht'


def array_sample_concordance_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/array_concordance/sample_concordance.ht'


def array_variant_concordance_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/array_concordance/variant_concordance.ht'


def ploidy_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    if data_source == 'broad' and freeze >= 5:
        return f'{sample_qc_prefix(data_source, freeze)}/sex_check/ploidy.ht'
    else:
        raise DataException("No ploidy file specified for this data_source and freeze yet")


def sex_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/sex_check/sex.ht'


def sex_tsv_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/sex_check/sex.tsv'


def qc_mt_path(data_source: str, freeze: int = CURRENT_FREEZE, ld_pruned: bool = False) -> str:
    """
    Returns path of MatrixTable for sample QC purposes
    :param bool ld_pruned: Should the qc matrix be LD pruned
    :return: Path MatrixTable for sample QC purposes
    :rtype: str
    """
    ld_pruned = '.pruned' if ld_pruned else ''
    return f'{sample_qc_prefix(data_source, freeze)}/qc_data/high_callrate_common_biallelic_snps{ld_pruned}.mt'


def qc_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/qc_data/high_callrate_common_biallelic_snps.ht'


def callrate_mt_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/platform_pca/callrate.mt'


def platform_pca_scores_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/platform_pca/platform_pca_scores.ht'


def platform_pca_loadings_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/platform_pca/platform_pca_loadings.ht'


def platform_pca_results_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/platform_pca/platform_pca_results.ht'


def relatedness_pca_scores_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/relatedness/pruned.pca_scores.ht'


def relatedness_ht_path(data_source: str, freeze: int = CURRENT_FREEZE, method: str = 'pc_relate') -> str:
    method = f'.{method}' if method != 'pc_relate' else ''
    return f'{sample_qc_prefix(data_source, freeze)}/relatedness/relatedness{method}.ht'


def duplicates_ht_path(data_source: str, freeze: int = CURRENT_FREEZE, dup_sets: bool = False, method: str = 'pc_relate') -> str:
    dup_sets = f'_sets' if dup_sets else ''
    method = f'.{method}' if method != 'pc_relate' else ''
    return f'{sample_qc_prefix(data_source, freeze)}/relatedness/duplicate{dup_sets}{method}.ht'


def inferred_ped_path(data_source: str, freeze: int = CURRENT_FREEZE, method: str = 'pc_relate') -> str:
    method = f'.{method}' if method != 'pc_relate' else ''
    return f'{sample_qc_prefix(data_source, freeze)}/relatedness/ped{method}.txt'


def related_drop_path(data_source: str, freeze: int = CURRENT_FREEZE, method: str = 'pc_relate') -> str:
    method = f'.{method}' if method != 'pc_relate' else ''
    return f'{sample_qc_prefix(data_source, freeze)}/relatedness/related_samples_to_drop{method}.ht'


def gnomad_ancestry_loadings_liftover_path(checkpoint: bool = False):
    if checkpoint:
        return 'gs://broad-ukbb/temp/gnomad_joint_unrelated_pca_loadings.ht'
    else:
        return 'gs://broad-ukbb/resources/gnomad.joint.unrelated.pca_loadings_lift.ht'


def ancestry_pca_scores_ht_path(data_source: str, freeze: int = CURRENT_FREEZE, population: str = None) -> str:
    pop = f'.{population}' if population else ''
    return f'{sample_qc_prefix(data_source, freeze)}/population_pca/pca_scores{pop}.ht'


def ancestry_pca_loadings_ht_path(data_source: str, freeze: int = CURRENT_FREEZE, population: str = None) -> str:
    pop = f'.{population}' if population else ''
    return f'{sample_qc_prefix(data_source, freeze)}/population_pca/pca_loadings{pop}.ht'


def ancestry_cluster_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/population_pca/cluster_assignments.ht'


def ancestry_cluster_array_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/population_pca/array_cluster_assignments.ht'


def ancestry_cluster_joint_scratch_array_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/population_pca/joint_scratch_array_cluster_assignments.ht'


def ancestry_hybrid_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/population_pca/hybrid_pop_assignments.ht'


def ancestry_pc_project_scores_ht_path(data_source: str, freeze: int = CURRENT_FREEZE, data_type: str = None) -> str:
    """
    Returns path of Table for scores and pop assignments from pc_project on gnomAD PCs
    :param str data_type: either None for UKBB only or joint for merged UKBB and gnomAD
    :return: Path to Table
    :rtype: str
    """
    data_type = f'.{data_type}' if data_type else ''
    return f'{sample_qc_prefix(data_source, freeze)}/population_pca/pc_project_scores_pop_assign{data_type}.ht'


def platform_pop_outlier_ht_path(data_source: str, freeze: int = CURRENT_FREEZE, pop_assignment_method: str = None) -> str:
    pop_assignment_method = f'.{pop_assignment_method}' if pop_assignment_method else ''
    return f'{sample_qc_prefix(data_source, freeze)}/outlier_detection/outlier_detection{pop_assignment_method}.ht'


def qc_temp_data_prefix(data_source: str, freeze: int = CURRENT_FREEZE):
    return f'{sample_qc_prefix(data_source, freeze)}/temp/'


def get_regeneron_relatedness_path(freeze: int = CURRENT_FREEZE, relationship: str = None) -> str:
    """
    Returns path to regeneron relatedness inference files
    :param int freeze: One of data freezes
    :param str relationship: relatedness relationship, need to be 2nd-degree, full-sibling, or parent-child
    :return: Output relatedness file path
    :rtype: str
    """

    if freeze not in FREEZES:
        raise DataException("This freeze is currently not present")
    if relationship not in ['2nd-degree', 'full-sibling', 'parent-child']:
        raise DataException("This regeneron relationship file not present")
    if freeze == 4:
        freeze_str = "Four"
    return f'gs://broad-ukbb/regeneron.freeze_{freeze}/data/pharma_relatedness_analysis/UKB_Freeze_{freeze_str}.NF.pVCF_{relationship}_relationships.genome'


def get_regeneron_broad_relatedness_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/relatedness/regeneron_joint_relatedness.ht'


def get_ukbb_self_reported_ancestry_path(freeze: int = CURRENT_FREEZE) -> str:
    return f'gs://broad-ukbb/resources/ukb24295.phenotypes.freeze_{freeze}.ht'


def get_regeneron_ancestry_path(freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to regeneron relatedness inference files
    :param int freeze: One of data freezes
    :param str relationship: relatedness relationship, need to be 2nd-degree, full-sibling, or parent-child
    :return: Output relatedness file path
    :rtype: str
    """

    if freeze not in FREEZES:
        raise DataException("This freeze is currently not present")
    if freeze == 4:
        freeze_str = "Four"
    return f'gs://broad-ukbb/regeneron.freeze_{freeze}/data/pharma_relatedness_analysis/UKB_Freeze_{freeze_str}.NF.splitmulti.commonsnps_samples_ancestries.txt'


def get_joint_regeneron_ancestry_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/population_pca/regeneron_ukb_joint_ancestry.ht'


def interval_qc_path(data_source: str, freeze: int = CURRENT_FREEZE, chrom: str = None, ht: bool = True) -> str:
    if chrom is None:
        chrom = ""
    else:
        chrom = f'.{chrom}'
    prefix = f'{sample_qc_prefix(data_source, freeze)}/interval_qc/coverage_by_target{chrom}'
    if ht:
        return prefix + '.ht'
    else:
        return prefix + '.txt'
