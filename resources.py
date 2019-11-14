import hail as hl
from typing import *

CURRENT_FREEZE = 4
DATA_SOURCES = ['regeneron', 'broad']
FREEZES = [4, 5]
CURRENT_HAIL_VERSION = "0.2"

def get_ukbb_data(data_source: str, freeze: int = CURRENT_FREEZE, adj: bool = False, split: bool = True,
                  raw: bool = False, non_refs_only: bool = False, meta_root: Optional[str] = None) -> hl.MatrixTable:
    """
    Wrapper function to get UKBB data as MT. By default, returns split hardcalls (with adj annotated but not filtered)
    :param str data_source: Will be regeneron or broad
    :param int freeze: One of data freezes
    :param bool adj: Whether the returned data should be filtered to adj genotypes
    :param bool split: Whether the dataset should be split (only applies to raw=False)
    :param bool raw: Whether to return the raw data (not recommended: unsplit, and no special consideration on sex chromosomes)
    :param bool non_refs_only: Whether to return the non-ref-genotype only MT (warning: no special consideration on sex chromosomes)
    :param str meta_root: Root annotation name for metadata (e.g., 'meta')
    :return: hardcalls dataset
    :rtype: MatrixTable
    """
    from gnomad_hail.utils import filter_to_adj

    if raw and split:
        raise DataException('No split raw data. Use of hardcalls is recommended.')

    if non_refs_only:
        mt = hl.read_matrix_table(get_ukbb_data_path(data_source, freeze, split=split, non_refs_only=non_refs_only))
    else:
        mt = hl.read_matrix_table(get_ukbb_data_path(data_source, freeze, hardcalls=not raw, split=split))

    if adj:
        mt = filter_to_adj(mt)

    if meta_root:
        meta_ht = hl.read_table(meta_ht_path(data_source, freeze))
        mt = mt.annotate_cols(**{meta_root: meta_ht[mt.s]})

    return mt


def get_ukbb_data_path(data_source: str, freeze: int = CURRENT_FREEZE, hardcalls: bool = False, split: bool = True,
                         non_refs_only: bool = False) -> str:
    """
    Wrapper function to get paths to UKBB data
    :param str data_source: Will be regeneron or broad
    :param int freeze: One of data freezes
    :param bool hardcalls: Whether hardcalls should be returned
    :param bool split: Whether the dataset should be split (applies to hardcalls and non_refs_only)
    :param bool non_refs_only: Whether non-ref-genotype only MT should be returned
    :return: Path to chosen MT
    :rtype: str
    """
    if hardcalls and non_refs_only:
        raise DataException('No dataset with hardcalls and non_refs_only')
    if data_source not in DATA_SOURCES:
        raise DataException("This data_source is currently not present")
    if hardcalls:
        return hardcalls_mt_path(data_source, freeze, split)
    elif non_refs_only:
        return non_refs_only_mt_path(data_source, freeze, split)
    else:
        return raw_mt_path(data_source, freeze)


def array_mt_path(liftover: bool = False) -> str:
    """
    :param bool liftover: One of data freezes
    :return: Path to array MT
    :rtype: str
    """
    return f'gs://broad-ukbb/resources/array/ukbb_array{"_liftover_GRCh38" if liftover else ""}.mt'


def raw_mt_path(data_source: str, freeze: int = CURRENT_FREEZE, is_temp = False) -> str:
    """
    Warning: unsplit and no special consideration on sex chromosomes
    :param str data_source: Will be regeneron or broad
    :param int freeze: One of data freezes
    :return: Path to chosen raw MT
    :rtype: str
    """
    tempstr = '.temp' if is_temp else ''
    if data_source == 'regeneron':
        return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/data/{data_source}.freeze_{freeze}.nf.mt'
    elif data_source == 'broad':
        return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/data/{data_source}.freeze_{freeze}{tempstr}.mt'


def hardcalls_mt_path(data_source: str, freeze: int = CURRENT_FREEZE, split: bool = True) -> str:
    """
    Warning: unsplit and no special consideration on sex chromosomes
    :param str data_source: Will be regeneron or broad
    :param int freeze: One of data freezes
    :param bool split: Whether the dataset should be split
    :return: Path to chosen hardcalls MT
    :rtype: str
    """
    return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/hardcalls/hardcalls{".split" if split else ""}.mt'


def non_refs_only_mt_path(data_source: str, freeze: int = CURRENT_FREEZE, split: bool = True) -> str:
    return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/non_refs_only/non_refs_only{".split" if split else ""}.mt'


def get_array_data_path(extension: str, chrom: str) -> str:
    """
    Get path to UKBB array data for a specific chromosome in plink format.
    :param str extension: file type bim, bam, or fam
    :param str chrom: Chromosome to return
    :return: Path to array data for given chromosome and file type
    :rtype: str
    """
    if extension == 'bed':
        return f'gs://broad-ukbb/resources/array/ukb_cal_{chrom}_v2.{extension}'
    elif extension == 'bim':
        return f'gs://broad-ukbb/resources/array/ukb_snp_{chrom}_v2.{extension}'
    elif extension == 'fam':
        return 'gs://broad-ukbb/resources/array/ukb26041_cal_chr22_v2_s488292.fam'


array_sample_map = 'gs://broad-ukbb/resources/array/Project_26041_bridge.csv'
ukbb_calling_intervals_path = 'gs://broad-ukbb/resources/ukbb_exome_calling.interval_list'
broad_calling_intervals_path = 'gs://broad-ukbb/resources/broad_exome_calling.interval_list'
lcr_intervals_path = 'gs://broad-ukbb/resources/LCRFromHengH38_chr1-22_XY.txt'
hg38_selfchain_path = 'gs://broad-ukbb/resources/hg38_self_chain_nosamepos_withalts_gt10k.bed.gz'
hg38_segdup_path = 'gs://broad-ukbb/resources/hg38.segdups_sorted_merged_gt10kb.bed.gz'


ukbb_calling_intervals_summary = 'gs://broad-ukbb/resources/ukbb_exome_calling_intervals.summary.txt'


def gnomad_ancestry_loadings_liftover_path(checkpoint: bool = False):
    if checkpoint:
        return 'gs://broad-ukbb/temp/gnomad_joint_unrelated_pca_loadings.ht'
    else:
        return 'gs://broad-ukbb/resources/gnomad.joint.unrelated.pca_loadings_lift.ht'


def array_sample_map_ht(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/array_sample_map.ht'


def get_lcr_intervals() -> hl.Table:
    return hl.import_locus_intervals(lcr_intervals_path, reference_genome='GRCh38', skip_invalid_intervals=True)

def get_selfchain_intervals() -> hl.Table:
    return hl.import_locus_intervals(hg38_selfchain_path, reference_genome='GRCh38', skip_invalid_intervals=True)

def get_segdup_intervals() -> hl.Table:
    return hl.import_locus_intervals(hg38_segdup_path, reference_genome='GRCh38', skip_invalid_intervals=True)


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


def get_ht_checkpoint_path(data_source: str, freeze: int = CURRENT_FREEZE, name: str = None) -> str:
    """
    Creates a checkpoint path for Table
    :param str data_type: Will be regeneron or broad
    :param int freeze: One of data freezes
    :param str name: Name of intermediate Table
    :return: Output checkpoint path
    :rtype: str
    """

    return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/{name}.ht'


def capture_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    if data_source == 'broad' and freeze == 4:
        return 'gs://broad-ukbb/resources/xgen_plus_spikein.Homo_sapiens_assembly38.targets.pad50.intervals.ht'
    elif data_source == 'regeneron':
        return 'gs://broad-ukbb/resources/ukbb_exome_calling_intervals.summary.ht'
    else:
        raise DataException("No interval file specified for this data_source and freeze yet")


def get_mt_checkpoint_path(data_source: str, freeze: int = CURRENT_FREEZE, name: str = None) -> str:
    """
    Creates a checkpoint path for Table
    :param str data_type: Will be regeneron or broad
    :param int freeze: One of data freezes
    :param str name: Name of intermediate Table
    :return: Output checkpoint path
    :rtype: str
    """

    return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/{name}.mt'


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


ukbb_phenotype_path = "gs://broad-ukbb/resources/ukb24295.phenotypes.txt"


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


def interval_qc_path(data_source: str, freeze: int = CURRENT_FREEZE, chrom: int = None, ht: bool = True) -> str:
    if chrom is None:
        chrom = ""
    else:
        chrom = f'.chr{chrom}'
    prefix = f'{sample_qc_prefix(data_source, freeze)}/interval_qc/coverage_by_target{chrom}'
    if ht:
        return prefix + '.ht'
    else:
        return prefix + '.txt'


# Variant QC annotations
def variant_qc_prefix(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    if data_source not in DATA_SOURCES:
        raise DataException("This data_source is currently not present")
    if freeze not in FREEZES:
        raise DataException("This freeze is currently not present")

    return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/variant_qc'


def var_annotations_ht_path(data_source: str, freeze: int, annotation_type: str) -> str:
    """
    Get variant-level annotations

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param str annotation_type: Tag describing variant-level annotations in HT, e.g., "vep", "qc_stats", "call_stats", "vqsr"
    :return: Path to annotations Table
    :rtype: str
    """
    return f'{variant_qc_prefix(data_source, freeze)}/variant_annotations/{annotation_type}.ht'


def sample_annotations_table_path(data_source: str, freeze: int, annotation_type: str) -> str:
    """
    Get sample-level annotations

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param str annotation_type: Tag describing variant-level annotations in HT, e.g., "family_stats"
    :return: Path to annotations HT
    :rtype: str
    """
    return f'{variant_qc_prefix(data_source, freeze)}/sample_annotations/{annotation_type}.ht'


def rf_run_hash_path(data_source: str, freeze: int = CURRENT_FREEZE):
    """
    Returns the path to the json file containing the RF runs list.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Path to json file
    :rtype: str
    """

    return f'{variant_qc_prefix(data_source, freeze)}/rf/runs.json'


def rf_annotated_path(
        data_source: str,
        freeze: int = CURRENT_FREEZE,
        adj: bool = False) -> str:
    """
    Returns the path to the RF-ready annotated HT

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param bool adj: Whether to load 'adj' or 'raw'
    :return:
    """

    return f'{variant_qc_prefix(data_source, freeze)}/rf/rf_annotated.{"adj" if adj else "raw"}.ht'


def rf_path(data_source: str, freeze: int = CURRENT_FREEZE,
            data: str = 'rf_result',
            run_hash: str = None
            ) -> str:
    """

    Gets the path to the desired RF data.
    Data can take the following values:
        - 'training': path to the training data for a given run
        - 'model': path to pyspark pipeline RF model
        - 'rf_result' (default): path to HT containing result of RF filtering

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param str data: One of 'pre_rf', 'training', 'model' or 'rf_result' (default)
    :param str run_hash: Hash of RF run to load
    :return:
    """

    extension = 'model' if data == 'model' else 'ht'
    return f'{variant_qc_prefix(data_source, freeze)}/rf/{run_hash}/{data}.{extension}'


def score_ranking_path(data_source: str, freeze: int,
            data: str,
            binned: bool = False,
            ) -> str:
    """
    Returns the path to non-RF metrics score rankings Tables, e.g.:
    * vqsr

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param data: The score data to return
    :param binned: Whether to get the binned data
    :return: Path to desired hail Table
    :rtype: str
    """

    return f'{variant_qc_prefix(data_source, freeze)}/score_rankings/{data}{"_binned" if binned else ""}.ht'

def binned_concordance_path(data_source: str, freeze: int, truth_sample: str, metric: str):
    '''
    :param str data_source: 'broad' or 'regeneron'
    :param int freeze: UKBB tranche version
    :param str truth_sample: Which truth sample concordance to analyze (e.g., "NA12878" or "syndip")
    :param str metric: One of the evaluation metrics (or a RF hash)
    :return: Path to Hail Table
    :rtype: str
    '''

    return f'{variant_qc_prefix(data_source, freeze)}/rf/{truth_sample}.{metric}.binned_concordance.ht'


def omni_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/1000G_omni2.5.hg38.mt'.format(hail_version)


def mills_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/Mills_and_1000G_gold_standard.indels.hg38.mt'.format(hail_version)


def hapmap_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/hapmap_3.3.hg38.mt'.format(hail_version)


def hapmap_ht_path():
    return 'gs://broad-ukbb/resources//hapmap_3.3.hg38.ht'

def dbsnp_ht_path():
    return 'gs://gnomad-public/resources/grch38/dbsnp_b151_grch38_all_20180418.ht'


def kgp_high_conf_snvs_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/1000G_phase1.snps.high_confidence.hg38.mt'.format(hail_version)


def release_prefix(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    if data_source not in DATA_SOURCES:
        raise DataException("This data_source is currently not present")
    if freeze not in FREEZES:
        raise DataException("This freeze is currently not present")
    return  f'gs://broad-ukbb/{data_source}.freeze_{freeze}/release'


def release_mt_path(data_source: str, freeze: int, nested=True, temp=False) -> str:
    '''
    Fetch filepath for release Hail MatrixTables

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param bool nested: If True, fetch Table in which variant annotations (e.g., freq, popmax, faf, and age histograms)
        are in array format ("nested"); if False, fetch Table in which nested variant annotations are unfurled
    :param bool temp: If True, fetch Table in which nested variant annotations are unfurled but listed under 'info' rather
        than at the top level; used for sanity-checking sites
    :return: Filepath for desired Hail Table
    :rtype: str
    '''
    tag = 'nested' if nested else 'flat'
    tag = tag + '.temp' if temp else tag
    return f'{release_prefix(data_source, freeze)}/mt/{tag}.mt'


def release_ht_path(data_source: str, freeze: int, nested=True, temp=False) -> str:
    '''
    Fetch filepath for release Hail Tables

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param bool nested: If True, fetch Table in which variant annotations (e.g., freq, popmax, faf, and age histograms)
        are in array format ("nested"); if False, fetch Table in which nested variant annotations are unfurled
    :param bool temp: If True, fetch Table in which nested variant annotations are unfurled but listed under 'info' rather
        than at the top level; used for sanity-checking sites
    :return: Filepath for desired Hail Table
    :rtype: str
    '''
    tag = 'nested' if nested else 'flat'
    tag = tag + '.temp' if temp else tag
    return f'{release_prefix(data_source, freeze)}/ht/{tag}.sites.ht'


def release_vcf_path(data_source: str, freeze: int, contig=None) -> str:
    '''
    Fetch filepath for release (variant-only) VCFs

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param str contig: String containing the name of the desired reference contig
    :return: Filepath for the desired VCF
    :rtype: str
    '''
    if contig:
        return f'{release_prefix(data_source, freeze)}/vcf/{contig}.vcf.bgz'
    else:
        return f'gs://gnomad-public/release/{release}/vcf/{data_type}/gnomad.{data_type}.{release_tag}.vcf.bgz'


class DataException(Exception):
    pass
