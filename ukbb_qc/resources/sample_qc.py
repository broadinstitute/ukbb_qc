from gnomad_hail.resources.resource_utils import DataException
from .resource_utils import CURRENT_FREEZE, DATA_SOURCES, FREEZES


def sample_qc_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to sample QC for specified data source and freeze

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to sample QC bucket for that data source and freeze
    :rtype: str
    """
    if data_source not in DATA_SOURCES:
        raise DataException("This data_source is currently not present")
    if freeze not in FREEZES:
        raise DataException("This freeze is currently not present")
    return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/sample_qc'


def qc_temp_data_prefix(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to temporary sample QC data. 
    Used to store densified MT data (freeze 6 and newer) or outlier detection data.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to bucket with temporary sample QC data
    :rtype: str
    """
    return f'{sample_qc_path(data_source, freeze)}/temp/'


# Interval resources
def interval_qc_path(data_source: str, freeze: int = CURRENT_FREEZE, chrom: str = None, ht: bool = True) -> str:
    """
    Returns path to interval QC results (either as a Table or text file).

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param str chrom: Chromosome of interest
    :param bool ht: Whether to return path to Table
    :return: Path to interval QC results
    :rtype: str
    """
    if chrom is None:
        chrom = ""
    else:
        chrom = f'.{chrom}'
    return f'{sample_qc_path(data_source, freeze)}/interval_qc/coverage_by_target{chrom}{".ht" if ht else ".txt"}'


def densified_interval_qc_path(
    data_source: str, freeze: int = CURRENT_FREEZE, repartitioned: bool = False
    ) -> str:
    """
    Returns path to MatrixTable densified only to interval QC sites

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool repartitioned: Whether MatrixTable has been repartitioned
    :return: Path to densified interval qc sites only MatrixTable
    :rtype: str
    """
    if freeze < 6:
        raise DataException("This resource doesn't exist for anything older than tranche 3/freeze 6")
    name = 'dense.repartitioned' if repartitioned else 'dense'
    return f'{qc_temp_data_prefix(data_source, freeze)}{name}.mt'


def f_stat_sites_path() -> str:
    """
    Returns path to Table with high callrate, common, pass interval QC, biallelic SNP positions on chromosome X used in sex imputation.
    NOTE: The sites came from freeze 5 (the last dataset with AF).

    :return: Path to Table with sites for sex imputation
    :rtype: str
    """
    data_source = 'broad'
    freeze = 5
    return f'{qc_temp_data_prefix(data_source, freeze)}f_stat_sites.ht'


def qc_sites_path() -> str:
    """
    Returns path to Table with sites to use in QC MatrixTable generation.
    NOTE: Sites were generated using QC MT from tranche 2/freeze 5.

    :return: Path to QC MT sites Table
    :rtype: str
    """
    data_source = 'broad'
    freeze = 5
    return f'{sample_qc_path(data_source, freeze)}/qc_data/qc_sites.ht'


# QC resources (meta ht, qc mt, qc ht, qc sites ht)
def meta_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to metadata Table

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to metadata Table
    :rtype: str
    """
    return f'{sample_qc_path(data_source, freeze)}/meta.ht'


def qc_mt_path(data_source: str, freeze: int = CURRENT_FREEZE, ld_pruned: bool = False) -> str:
    """
    Returns path to MatrixTable filtered to high callrate, common, biallelc snps for sample QC purposes

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool ld_pruned: Whether to return the LD pruned version of the MatrixTable
    :return: Path MatrixTable for sample QC purposes
    :rtype: str
    """
    if ld_pruned and freeze > 5:
        raise DataException("Only one version of QC MT stored for tranches after tranche 2 (created using tranche 2 sites)")
    ld_pruned = '.pruned' if ld_pruned else ''
    return f'{sample_qc_path(data_source, freeze)}/qc_data/high_callrate_common_biallelic_snps{ld_pruned}.mt'


def qc_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:  
    """
    Returns path to Table for sample QC purposes (specifically, relatedness checks)

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to Table for sample QC purposes
    :rtype: str
    """
    return f'{sample_qc_path(data_source, freeze)}/qc_data/high_callrate_common_biallelic_snps.ht'


# Hard filters resources (including sex imputation resources)
def ploidy_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with sex chromosome ploidies

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to ploidy Table
    :rtype: str
    """
    if data_source == 'broad' and freeze >= 5:
        return f'{sample_qc_path(data_source, freeze)}/sex_check/ploidy.ht'
    else:
        raise DataException("No ploidy file specified for this data source and freeze yet")


def sex_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with inferred sex karyotypes

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to sex Table
    :rtype: str
    """
    return f'{sample_qc_path(data_source, freeze)}/sex_check/sex.ht'


def hard_filters_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with hard filter status

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to hard filters Table
    :rtype: str
    """
    return f'{sample_qc_path(data_source, freeze)}/hard_filters_flagged.ht'


# Array resources
def array_sample_concordance_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with array sample concordance status

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to array sample concordance Table
    :rtype: str
    """
    return f'{sample_qc_path(data_source, freeze)}/array_concordance/sample_concordance.ht'


def array_variant_concordance_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with array variant concordance status

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to array variant concordance Table
    :rtype: str
    """
    return f'{sample_qc_path(data_source, freeze)}/array_concordance/variant_concordance.ht'


# Platform PCA resources
def callrate_mt_path(
    data_source: str, freeze: int = CURRENT_FREEZE, 
    interval_filtered: bool = False
    ) -> str:
    """
    Returns path to callrate MatrixTable 

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool interval_filtered: Whether data was filtered to high coverage intervals
    :return: Path to array variant concordance Table
    :rtype: str
    """
    filtered = '.interval_filtered.' if interval_filtered else ''
    return f'{sample_qc_path(data_source, freeze)}/platform_pca/callrate{filtered}.mt'


def platform_pca_scores_ht_path(
    data_source: str, freeze: int = CURRENT_FREEZE, 
    interval_filtered: bool = False
    ) -> str:
    """
    Returns path to Table with platform PCA scores

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool interval_filtered: Whether data was filtered to high coverage intervals
    :return: Path to platform PCA scores Table
    :rtype: str
    """
    filtered = '.interval_filtered.' if interval_filtered else ''
    return f'{sample_qc_path(data_source, freeze)}/platform_pca/platform_pca_scores{filtered}.ht'


def platform_pca_loadings_ht_path(
    data_source: str, freeze: int = CURRENT_FREEZE,
    interval_filtered: bool = False
    ) -> str:
    """
    Returns path to Table with platform PCA loadings

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool interval_filtered: Whether data was filtered to high coverage intervals
    :return: Path to platform PCA loadings Table
    :rtype: str
    """
    filtered = '.interval_filtered.' if interval_filtered else ''
    return f'{sample_qc_path(data_source, freeze)}/platform_pca/platform_pca_loadings{filtered}.ht'


def platform_pca_results_ht_path(
    data_source: str, freeze: int = CURRENT_FREEZE,
    interval_filtered: bool = False
    ) -> str:
    """
    Returns path to Table with platform PCA results

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool interval_filtered: Whether data was filtered to high coverage intervals
    :return: Path to platform PCA results Table
    :rtype: str
    """
    filtered = '.interval_filtered.' if interval_filtered else ''
    return f'{sample_qc_path(data_source, freeze)}/platform_pca/platform_pca_results{filtered}.ht'


# Relatedness inference resources
def relatedness_pca_scores_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with relatedness PCA scores (run prior to PC-relate)

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to relatedness PCA scores Table
    :rtype: str
    """
    return f'{sample_qc_path(data_source, freeze)}/relatedness/pruned.pca_scores.ht'


def relatedness_ht_path(data_source: str, freeze: int = CURRENT_FREEZE, method: str = 'pc_relate') -> str:
    """
    Returns path to Table output by PC-relate with inferred relationships

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param str method: Method of inferring relatedness. Default: hail's pc_relate
    :return: Path to relatedness Table
    :rtype: str
    """

    method = f'.{method}' if method != 'pc_relate' else ''
    return f'{sample_qc_path(data_source, freeze)}/relatedness/relatedness{method}.ht'


def duplicates_ht_path(data_source: str, freeze: int = CURRENT_FREEZE, dup_sets: bool = False, method: str = 'pc_relate') -> str:
    """
    Returns path to Table with inferred duplicates

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param str method: Method of inferring relatedness. Default: hail's pc_relate
    :return: Path to duplicates Table
    :rtype: str
    """
    dup_sets = f'_sets' if dup_sets else ''
    method = f'.{method}' if method != 'pc_relate' else ''
    return f'{sample_qc_path(data_source, freeze)}/relatedness/duplicate{dup_sets}{method}.ht'


def inferred_ped_path(data_source: str, freeze: int = CURRENT_FREEZE, method: str = 'pc_relate') -> str:
    """
    Returns path to Table with inferred families

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param str method: Method of inferring relatedness. Default: hail's pc_relate
    :return: Path to inferred family Table
    :rtype: str
    """

    method = f'.{method}' if method != 'pc_relate' else ''
    return f'{sample_qc_path(data_source, freeze)}/relatedness/ped{method}.txt'


def related_drop_path(data_source: str, freeze: int = CURRENT_FREEZE, method: str = 'pc_relate') -> str:
    """
    Returns path to Table with related samples to drop

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param str method: Method of inferring relatedness. Default: hail's pc_relate
    :return: Path to related samples to drop Table
    :rtype: str
    """
    method = f'.{method}' if method != 'pc_relate' else ''
    return f'{sample_qc_path(data_source, freeze)}/relatedness/related_samples_to_drop{method}.ht'


# Broad-Regneron relatedness comparison resources (only relevant for freeze 4)
def get_regeneron_relatedness_path(freeze: int = CURRENT_FREEZE, relationship: str = None) -> str:
    """
    Returns path to regeneron relatedness inference files

    :param int freeze: One of data freezes
    :param str relationship: relatedness relationship, need to be 2nd-degree, full-sibling, or parent-child
    :return: Output relatedness file path
    :rtype: str
    """
    if freeze != 4:
        raise DataException("Regeneron relatedness only exists for tranche 1/freeze 4")
    if relationship not in ['2nd-degree', 'full-sibling', 'parent-child']:
        raise DataException("This regeneron relationship file not present")
    freeze_str = "Four"
    return f'gs://broad-ukbb/regeneron.freeze_{freeze}/data/pharma_relatedness_analysis/UKB_Freeze_{freeze_str}.NF.pVCF_{relationship}_relationships.genome'


def get_regeneron_broad_relatedness_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with both Broad and Regeneron inferred relationships

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to Table with Broad and Regeneron inferred relationships
    :rtype: str
    """
    if freeze != 4:
        raise DataException("Broad-Regeneron relatedness comparisons only exists for tranche 1/freeze 4")
    return f'{sample_qc_path(data_source, freeze)}/relatedness/regeneron_joint_relatedness.ht'


# Population inference resources
def gnomad_ancestry_loadings_liftover_path(checkpoint: bool = False):
    """
    Returns path to gnomAD ancestry loadings lifted to build 38

    :param bool checkpoint: Whether to return path to temporary Table
    :return: Path to gnomAD ancestry loadings Tablee
    :rtype: str
    """
    if checkpoint:
        return 'gs://broad-ukbb/temp/gnomad_joint_unrelated_pca_loadings.ht'
    else:
        return 'gs://broad-ukbb/resources/gnomad.joint.unrelated.pca_loadings_lift.ht'


def ancestry_pca_scores_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with ancestry PCA scores from scratch PCA 

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to ancestry PCA scores Table
    :rtype: str
    """
    return f'{sample_qc_path(data_source, freeze)}/population_pca/pca_scores.ht'


def ancestry_pca_loadings_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with ancestry PCA loadings from scratch PCA 

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to ancestry PCA loadings Table
    :rtype: str
    """
    return f'{sample_qc_path(data_source, freeze)}/population_pca/pca_loadings.ht'


def ancestry_cluster_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with ancestry PCA cluster assignments using exome data

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to ancestry cluster assignments Table
    :rtype: str
    """
    return f'{sample_qc_path(data_source, freeze)}/population_pca/cluster_assignments.ht'


def ancestry_cluster_array_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with ancestry PCA cluster assignments using array data

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to ancestry cluster assingments Table using array data
    :rtype: str
    """
    return f'{sample_qc_path(data_source, freeze)}/population_pca/array_cluster_assignments.ht'


def ancestry_cluster_joint_scratch_array_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with ancestry PCA cluster assignments using exome and array data

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to ancestry cluster assingments Table using exome and array data
    :rtype: str
    """
    return f'{sample_qc_path(data_source, freeze)}/population_pca/joint_scratch_array_cluster_assignments.ht'


def ancestry_pc_project_scores_ht_path(data_source: str, freeze: int = CURRENT_FREEZE, data_type: str = None) -> str:
    """
    Returns path of Table for scores and pop assignments from pc_project on gnomAD PCs

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param str data_type: either None for UKBB only or joint for merged UKBB and gnomAD
    :return: Path to Table
    :rtype: str
    """
    data_type = f'.{data_type}' if data_type else ''
    return f'{sample_qc_path(data_source, freeze)}/population_pca/pc_project_scores_pop_assign{data_type}.ht'


def ancestry_hybrid_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with hybrid scratch and gnomAD PC project ancestry assignments

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to Table with hybrid ancestry assignments
    :rtype: str
    """
    return f'{sample_qc_path(data_source, freeze)}/population_pca/hybrid_pop_assignments.ht'


# UKBB self reported ancestry
def get_ukbb_self_reported_ancestry_path(freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with UKBB self reported ancestries

    :param int freeze: One of data freezes
    :return: Path to Table with self reported ancestry
    :rtype: str
    """
    if freeze not in FREEZES:
        raise DataException("This freeze doesn't exist yet")
    return f'gs://broad-ukbb/resources/ukb24295.phenotypes.freeze_{freeze}.ht'


# Broad-Regeneron ancestry comparison resources (relevant for freeze 4 only)
def get_regeneron_ancestry_path(freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Regeneron relatedness inference files

    :param int freeze: One of data freezes
    :param str relationship: relatedness relationship, need to be 2nd-degree, full-sibling, or parent-child
    :return: Output relatedness file path
    :rtype: str
    """

    if freeze != 4:
        raise DataException("Regeneron ancestry only exists for tranche 1/freeze 4")
    freeze_str = "Four"
    return f'gs://broad-ukbb/regeneron.freeze_{freeze}/data/pharma_relatedness_analysis/UKB_Freeze_{freeze_str}.NF.splitmulti.commonsnps_samples_ancestries.txt'


def get_joint_regeneron_ancestry_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with both Broad and Regeneron inferred ancestries

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to Table with Broad and Regeneron inferred ancestries
    :rtype: str
    """
    if freeze != 4:
        raise DataException("Regeneron ancestry only exists for tranche 1/freeze 4")
    return f'{sample_qc_path(data_source, freeze)}/population_pca/regeneron_ukb_joint_ancestry.ht'


# Outlier detection resources
def platform_pop_outlier_ht_path(data_source: str, freeze: int = CURRENT_FREEZE, pop_assignment_method: str = None) -> str:
    """
    Returns path to Table containing samples flagged for outlier sample QC metrics

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param str pop_assignment_method: Method used to infer populations (gnomad_qc_project_pop, HDBSCAN_pop_cluster, hybrid_pop)
    :return: Path to Table with outlier samples flagged
    :rtype: str
    """
    pop_assignment_method = f'.{pop_assignment_method}' if pop_assignment_method else ''
    return f'{sample_qc_path(data_source, freeze)}/outlier_detection/outlier_detection{pop_assignment_method}.ht'
