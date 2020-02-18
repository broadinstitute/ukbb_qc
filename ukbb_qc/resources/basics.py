import hail as hl
from typing import Optional
from gnomad_hail.utils.generic import file_exists
from gnomad_hail.resources.resource_utils import DataException
from .resource_utils import CURRENT_FREEZE, DATA_SOURCES, FREEZES, CURRENT_HAIL_VERSION
from .sample_qc import meta_ht_path


broad_calling_intervals_path = 'gs://broad-ukbb/resources/broad_exome_calling.interval_list'
ukbb_calling_intervals_summary = 'gs://broad-ukbb/resources/ukbb_exome_calling_intervals.summary.txt'
ukbb_phenotype_path = "gs://broad-ukbb/resources/ukb24295.phenotypes.txt"


# UKBB data resources
def excluded_samples_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to list of samples to exclude from QC due to withdrawn consents

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to excluded samples list
    :rtype: str
    """
    return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/data/samples_to_exclude.txt'


def get_ukbb_data(data_source: str, freeze: int = CURRENT_FREEZE, key_by_locus_and_alleles: bool = False,
                    adj: bool = False, split: bool = True, raw: bool = False, 
                    non_refs_only: bool = False, meta_root: Optional[str] = None
    ) -> hl.MatrixTable:
    """
    Wrapper function to get UKBB data as MatrixTable. By default, returns split hardcalls (with adj annotated but not filtered).

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool key_by_locus_and_alleles: Whether to key the MatrixTable by locus and alleles
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
    if not raw and not split:
        raise DataException('No unsplit hardcalls. Use of split hardcalls is recommended.')

    if non_refs_only:
        mt = hl.read_matrix_table(get_ukbb_data_path(data_source, freeze, split=split, non_refs_only=non_refs_only))
    else:
        mt = hl.read_matrix_table(get_ukbb_data_path(data_source, freeze, hardcalls=not raw, split=split))

    if adj:
        mt = filter_to_adj(mt)

    if meta_root:
        meta_ht = hl.read_table(meta_ht_path(data_source, freeze))
        mt = mt.annotate_cols(**{meta_root: meta_ht[mt.s]})

    if key_by_locus_and_alleles:
        mt = hl.MatrixTable(hl.ir.MatrixKeyRowsBy(mt._mir, ['locus', 'alleles'], is_sorted=True)) # taken from v3 qc

    if file_exists(excluded_samples_path):
        ht = hl.import_table(excluded_samples_path, no_header=True, impute=True)
        excluded_samples = hl.literal(ht.aggregate(hl.agg.collect_as_set(ht.f0)))
        mt = mt.filter_cols(~excluded_samples.contains(mt.s))
        mt = mt.filter_rows(hl.agg.any(mt.LGT.is_non_ref()) & hl.agg.any(hl.is_defined(mt.END)))

    return mt


def get_ukbb_data_path(data_source: str, freeze: int = CURRENT_FREEZE, hardcalls: bool = False, split: bool = True,
                         non_refs_only: bool = False) -> str:
    """
    Wrapper function to get paths to UKBB data.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool hardcalls: Whether to return hardcalls
    :param bool split: Whether the dataset should be split (applies to hardcalls and non_refs_only)
    :param bool non_refs_only: Whether non-ref-genotype only MatrixTable should be returned
    :return: Path to chosen MatrixTable
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


def array_mt_path(liftover: bool = False, checkpoint: bool = False) -> str:
    """
    Function to return path to UKBB array MatrixTable

    :param bool liftover: Whether to return the path to the version of the data lifted over to build 38
    :param bool checkpoint: Whether to return the temporary checkpoint path 
    :return: Path to array MatrixTable
    :rtype: str
    """
    if checkpoint:
        return f'gs://broad-ukbb/temp/ukbb_array{"_liftover_GRCh38" if liftover else ""}.mt'
    else:
        return f'gs://broad-ukbb/resources/array/ukbb_array{"_liftover_GRCh38" if liftover else ""}.mt'


def raw_mt_path(data_source: str, freeze: int = CURRENT_FREEZE, is_temp = False, densified: bool = False) -> str:
    """
    Returns path to raw UKBB MatrixTable. WARNING: unsplit and no special consideration on sex chromosomes

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool is_temp: Whether to return path to temp raw MatrixTable
    :param bool densified: Whether to return the densified sparse MatrixTable. Only relevant for freeze 5.
    :return: Path to raw MatrixTable
    :rtype: str
    """
    tempstr = '.temp' if is_temp else ''
    if data_source == 'regeneron':
        return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/data/{data_source}.freeze_{freeze}.nf.mt'
    elif data_source == 'broad':
        # create dict for raw mts for freeze 5 and later
        dsp_prefix = 'gs://broad-pharma5-ukbb-outputs'
        raw_mt_names = {
            5: 'hail_dataproc_20191108115937',
            6: 'hail_dataproc_20200130092005.mt'
        }
        if freeze == 4 or (freeze == 5 and densified):
            return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/data/{data_source}.freeze_{freeze}{tempstr}.mt'
        else:
            return f'{dsp_prefix}/{raw_mt_names[freeze]}'


def hardcalls_mt_path(data_source: str, freeze: int = CURRENT_FREEZE, split: bool = True) -> str:
    """
    Returns path to hardcalls MatrixTable.
    
    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool split: Whether the dataset should be split
    :return: Path to chosen hardcalls MatrixTable
    :rtype: str
    """
    return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/hardcalls/hardcalls{".split" if split else ""}.mt'


def non_refs_only_mt_path(data_source: str, freeze: int = CURRENT_FREEZE, split: bool = True) -> str:
    """
    Returns path to non-reference genotypes only MatrixTable

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool split: Whether the dataset should be split
    :return: Path to non-ref only MatrixTable
    :rtype: str
    """
    return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/non_refs_only/non_refs_only{".split" if split else ""}.mt'


def get_checkpoint_path(data_source: str, freeze: int = CURRENT_FREEZE, name: str = None, mt: bool = False) -> str:
    """
    Creates a checkpoint path for Table or MatrixTable

    :param str data_type: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param str name: Name of intermediate Table/MatrixTable
    :param bool mt: Whether path is for a MatrixTable, default is False
    :return: Output checkpoint path
    :rtype: str
    """
    return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/{name}.{"mt" if mt else "ht"}'


# Sparse MatrixTable resources
def last_END_positions_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Retrieves path to Table with last END positions annotation

    :param str data_type: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to Table with last_END_position annotation
    :rtype: str
    """
    return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/sparse_resources/last_END_positions.ht'


# Array resources
def get_array_data_path(extension: str, chrom: str) -> str:
    """
    Get path to UKBB array data for a specific chromosome in plink format.

    :param str extension: File type: bim, bam, or fam
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


def array_sample_map_path(freeze: int = CURRENT_FREEZE) -> str:
    """
    Get path to UKBB array sample mapping csv file.

    :param int freeze: One of data freezes
    :return: Path to array sample map csv
    :rtype: str
    """
    if freeze == 4:
        return 'gs://broad-ukbb/resources/array/Project_26041_bridge.csv'
    elif freeze == 5:
        return 'gs://broad-ukbb/resources/array/linking_file_200K_withbatch.csv'


def array_sample_map_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> hl.Table:
    """
    Returns array sample map Table.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: array sample map Table
    :rtype: hl.Table
    """
    return f'gs://broad-ukbb/{data_source}.freeze_{freeze}/array_sample_map.ht'


# Capture intervals 
def capture_ht_path(data_source: str) -> str:
    """
    Returns path to capture intervals Table.

    :param str data_source: One of 'regeneron' or 'broad'
    :return: Path to capture intervals Table
    :rtype: str
    """
    if data_source == 'broad':
        return 'gs://broad-ukbb/resources/xgen_plus_spikein.Homo_sapiens_assembly38.targets.pad50.merged.ht'
    elif data_source == 'regeneron':
        return 'gs://broad-ukbb/resources/ukbb_exome_calling_intervals.summary.ht'
    else:
        raise DataException("This data_source is currently not present")


# Release resources
def get_release_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to release data.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to release data bucket
    :rtype: str
    """
    if data_source not in DATA_SOURCES:
        raise DataException("This data_source is currently not present")
    if freeze not in FREEZES:
        raise DataException("This freeze is currently not present")
    return  f'gs://broad-ukbb/{data_source}.freeze_{freeze}/release'


def release_mt_path(data_source: str, freeze: int, temp=False) -> str:
    '''
    Fetch filepath for release Hail MatrixTables

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param bool temp: If True, fetch Table in which nested variant annotations are unfurled but listed under 'info' rather
        than at the top level; created as a checkpoint during VCF generation
    :return: Filepath for release MatrixTable
    :rtype: str
    '''
    tag = '.temp' if temp else ''
    return f'{get_release_path(data_source, freeze)}/mt/nested{tag}.mt'


def release_ht_path(data_source: str, freeze: int) -> str:
    '''
    Fetch filepath for release Hail Tables

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Filepath for release Table
    :rtype: str
    '''
    return f'{get_release_path(data_source, freeze)}/ht/release.sites.ht'


def release_vcf_path(data_source: str, freeze: int, contig: str = None) -> str:
    '''
    Fetch bucket for release (variant-only) VCFs

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param str contig: String containing the name of the desired reference contig
    :return: Filepath for the desired VCF
    :rtype: str
    '''
    if contig:
        return f'{get_release_path(data_source, freeze)}/vcf/{contig}.vcf.bgz'
    else:
        return f'{get_release_path(data_source, freeze)}/vcf/sharded_vcf/ukbb_freeze_{freeze}_release.vcf.bgz'
