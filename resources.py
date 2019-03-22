import hail as hl

CURRENT_FREEZE = 4
DATA_SOURCES = ['regeneron']
FREEZES = [4]


def get_ukbb_data(data_source: str, freeze: int = CURRENT_FREEZE, adj: bool = False, split: bool = True,
                  raw: bool = False, non_refs_only: bool = False) -> hl.MatrixTable:
    """
    Wrapper function to get UKBB data as MT. By default, returns split hardcalls (with adj annotated but not filtered)
    :param str data_source: Will be regeneron or broad
    :param int freeze: One of data freezes
    :param bool adj: Whether the returned data should be filtered to adj genotypes
    :param bool split: Whether the dataset should be split (only applies to raw=False)
    :param bool raw: Whether to return the raw data (not recommended: unsplit, and no special consideration on sex chromosomes)
    :param bool non_refs_only: Whether to return the non-ref-genotype only MT (warning: no special consideration on sex chromosomes)
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


def raw_mt_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Warning: unsplit and no special consideration on sex chromosomes
    :param str data_source: Will be regeneron or broad
    :param int freeze: One of data freezes
    :return: Path to chosen raw MT
    :rtype: str
    """
    # Note: nf wont apply to broad so should we just remove nf or add as a parameter?
    return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/data/{data_source}.freeze_{freeze}.nf.mt'


def hardcalls_mt_path(data_source: str, freeze: int = CURRENT_FREEZE, split: bool = True) -> str:
    """
    Warning: unsplit and no special consideration on sex chromosomes
    :param str data_source: Will be regeneron or broad
    :param int freeze: One of data freezes
    :param bool split: Whether the dataset should be split
    :return: Path to chosen hardcalls MT
    :rtype: str
    """
    return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/data/{data_source}.freeze_{freeze}.hardcalls{".split" if split else ""}.mt'


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


def array_sample_concordance_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/array_concordance/sample_concordance.ht'

def array_variant_concordance_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/array_concordance/variant_concordance.ht'

def sex_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/sex_check/sex.ht'

def qc_mt_path(data_source: str, freeze: int = CURRENT_FREEZE, ld_pruned: bool = False) -> str:
    """
    Returns path of MatrixTable for sample QC purposes
    :param bool ld_pruned: Should the qc matrix be LD pruned
    :return: Path MatrixTable for sample QC purposes
    :rtype: str
    """
    ld_pruned = '.pruned' if ld_pruned else ''
    return f'{sample_qc_prefix(data_source, freeze)}/high_callrate_common_biallelic_snps{ld_pruned}.mt'

def qc_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/high_callrate_common_biallelic_snps.ht'
    return f'{sample_qc_prefix(data_source, freeze)}/high_callrate_common_biallelic_snps.ht'

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

def relatedness_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/relatedness/relatedness.ht'

def duplicates_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/relatedness/dups.ht'

def inferred_ped_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/relatedness/ped.txt'

def related_drop_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    return f'{sample_qc_prefix(data_source, freeze)}/relatedness/related_samples_to_drop.ht'

def ancestry_pca_scores_ht_path(data_source: str, freeze: int = CURRENT_FREEZE, population: str = None) -> str:
    pop = f'.{population}' if population else ''
    return f'{sample_qc_prefix(data_source, freeze)}/ancestry/unrelated.pca_scores{pop}.ht'

def ancestry_pca_loadings_ht_path(data_source: str, freeze: int = CURRENT_FREEZE, population: str = None) -> str:
    pop = f'.{population}' if population else ''
    return f'{sample_qc_prefix(data_source, freeze)}/ancestry/unrelated.pca_loadings{pop}.ht'


class DataException(Exception):
    pass
