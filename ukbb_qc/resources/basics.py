import hail as hl
from typing import Optional
from gnomad.utils.generic import file_exists
from gnomad.resources.resource_utils import DataException
from .resource_utils import CURRENT_FREEZE, DATA_SOURCES, FREEZES, CURRENT_HAIL_VERSION
from .sample_qc import meta_ht_path
from ukbb_qc.assessment.sanity_checks import sample_check


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("basics")
logger.setLevel(logging.INFO)


broad_calling_intervals_path = (
    "gs://broad-ukbb/resources/broad_exome_calling.interval_list"
)
ukbb_calling_intervals_summary = (
    "gs://broad-ukbb/resources/ukbb_exome_calling_intervals.summary.txt"
)


# UKBB data resources
def excluded_samples_path() -> str:
    """
    Returns path to list of samples to exclude from QC due to withdrawn consents

    :return: Path to excluded samples list
    :rtype: str
    """
    return f"gs://broad-ukbb/resources/withdrawn_consents/w26041_20200204.csv"


def get_ukbb_data(
    data_source: str,
    freeze: int = CURRENT_FREEZE,
    key_by_locus_and_alleles: bool = False,
    adj: bool = False,
    split: bool = True,
    raw: bool = False,
    ukbb_samples_only: bool = True,
    meta_root: Optional[str] = None,
) -> hl.MatrixTable:
    """
    Wrapper function to get UKBB data as MatrixTable. By default, returns split hardcalls (with adj annotated but not filtered).

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool key_by_locus_and_alleles: Whether to key the MatrixTable by locus and alleles
    :param bool adj: Whether the returned data should be filtered to adj genotypes
    :param bool split: Whether the dataset should be split (only applies to raw=False)
    :param bool raw: Whether to return the raw data (not recommended: unsplit, and no special consideration on sex chromosomes)
    :param bool ukbb_samples_only: Whether to return only UKBB samples (exclude control samples). Default is True.
    :param str meta_root: Root annotation name for metadata (e.g., 'meta')
    :return: hardcalls dataset
    :rtype: MatrixTable
    """
    from gnomad.utils import filter_to_adj

    if raw and split:
        raise DataException("No split raw data. Use of hardcalls is recommended.")
    if not raw and not split:
        raise DataException(
            "No unsplit hardcalls. Use of split hardcalls is recommended."
        )

    mt = hl.read_matrix_table(
        get_ukbb_data_path(data_source, freeze, hardcalls=not raw, split=split)
    )

    if adj:
        mt = filter_to_adj(mt)

    if meta_root:
        meta_ht = hl.read_table(meta_ht_path(data_source, freeze))
        mt = mt.annotate_cols(**{meta_root: meta_ht[mt.s]})

    if key_by_locus_and_alleles:
        mt = hl.MatrixTable(
            hl.ir.MatrixKeyRowsBy(mt._mir, ["locus", "alleles"], is_sorted=True)
        )  # taken from v3 qc

    # Remove any samples with withdrawn consents if file with excluded samples exists
    if file_exists(excluded_samples_path(freeze)):
        if not file_exists(f"{array_sample_map_ht_path(data_source, freeze)}_SUCCESS"):
            raise DataException(
                f"Need to import array sample map ht for freeze {freeze}!"
            )
        else:
            sample_map_ht = hl.read_table(array_sample_map_ht_path(data_source, freeze))
        mt = mt.annotate_cols(
            withdrawn_consent=sample_map_ht[mt.s.split("_")[1]].withdrawn_consent
        )
        mt = mt.filter_cols(~mt.withdrawn_consent)
        mt = mt.drop("withdrawn_consent")

        # Double check all withdrawn samples were actually excluded
        withdrawn_ht = hl.import_table(
            excluded_samples_path(freeze), no_header=True, impute=True,
        )
        mt_samples = mt.cols()
        mt_samples = mt_samples.key_by(s_id=mt.s.split("_")[1])

        # Call sample_check function
        # If withdrawn_ids_missing is False, then some samples with withdrawn consents are still present in MT
        withdrawn_ids_missing, extra_mt_ids = sample_check(
            mt_samples, withdrawn_ht, show_mismatch=False
        )
        if not withdrawn_ids_missing:
            raise DataException(
                "Withdrawn samples present in MT. Double check sample filtration"
            )
        else:
            logger.info("No withdrawn samples found in MT")

    if ukbb_samples_only:
        mt = mt.filter_cols(mt.s.contains("UKB"))

    gt_expr = mt.GT if split else mt.LGT
    mt = mt.filter_rows(
        hl.agg.any(gt_expr.is_non_ref()) | hl.agg.any(hl.is_defined(mt.END))
    )

    return mt


def get_ukbb_data_path(
    data_source: str,
    freeze: int = CURRENT_FREEZE,
    hardcalls: bool = False,
    split: bool = True,
) -> str:
    """
    Wrapper function to get paths to UKBB data.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool hardcalls: Whether to return hardcalls
    :param bool split: Whether the dataset should be split (applies to hardcalls only)
    :return: Path to chosen MatrixTable
    :rtype: str
    """
    if data_source not in DATA_SOURCES:
        raise DataException("This data_source is currently not present")
    if hardcalls:
        return hardcalls_mt_path(data_source, freeze)
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


def raw_mt_path(
    data_source: str,
    freeze: int = CURRENT_FREEZE,
    is_temp=False,
    densified: bool = False,
) -> str:
    """
    Returns path to raw UKBB MatrixTable. WARNING: unsplit and no special consideration on sex chromosomes

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool is_temp: Whether to return path to temp raw MatrixTable
    :param bool densified: Whether to return the densified sparse MatrixTable. Only relevant for freeze 5.
    :return: Path to raw MatrixTable
    :rtype: str
    """
    tempstr = ".temp" if is_temp else ""
    if data_source == "regeneron":
        return f"gs://broad-ukbb/{data_source}.freeze_{freeze}/data/{data_source}.freeze_{freeze}.nf.mt"
    elif data_source == "broad":
        # create dict for raw mts for freeze 5 and later
        dsp_prefix = "gs://broad-pharma5-ukbb-outputs"
        raw_mt_names = {
            5: "hail_dataproc_20191108115937",
            6: "hail_dataproc_20200130092005.mt",
        }
        if freeze == 4 or (freeze == 5 and densified):
            return f"gs://broad-ukbb/{data_source}.freeze_{freeze}/data/{data_source}.freeze_{freeze}{tempstr}.mt"
        else:
            return f"{dsp_prefix}/{raw_mt_names[freeze]}"


def hardcalls_mt_path(data_source: str, freeze: int = CURRENT_FREEZE,) -> str:
    """
    Returns path to hardcalls MatrixTable.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool split: Whether the dataset should be split
    :return: Path to chosen hardcalls MatrixTable
    :rtype: str
    """
    return f"gs://broad-ukbb/{data_source}.freeze_{freeze}/hardcalls/hardcalls.mt"


def get_checkpoint_path(
    data_source: str, freeze: int = CURRENT_FREEZE, name: str = None, mt: bool = False
) -> str:
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
def last_END_positions_ht_path(
    data_source: str = "broad", freeze: int = CURRENT_FREEZE
) -> str:
    """
    Retrieves path to Table with last END positions annotation

    :param str data_type: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to Table with last_END_position annotation
    :rtype: str
    """
    return f"gs://broad-ukbb/{data_source}.freeze_{freeze}/sparse_resources/last_END_positions.ht"


# Array resources
def get_array_data_path(extension: str, chrom: str) -> str:
    """
    Get path to UKBB array data for a specific chromosome in plink format.

    :param str extension: File type: bim, bam, or fam
    :param str chrom: Chromosome to return
    :return: Path to array data for given chromosome and file type
    :rtype: str
    """
    if extension == "bed":
        return f"gs://broad-ukbb/resources/array/ukb_cal_{chrom}_v2.{extension}"
    elif extension == "bim":
        return f"gs://broad-ukbb/resources/array/ukb_snp_{chrom}_v2.{extension}"
    elif extension == "fam":
        return "gs://broad-ukbb/resources/array/ukb26041_cal_chr22_v2_s488292.fam"


def array_sample_map_path(freeze: int = CURRENT_FREEZE) -> str:
    """
    Get path to UKBB array sample mapping csv file.

    :param int freeze: One of data freezes
    :return: Path to array sample map csv
    :rtype: str
    """
    array_map_names = {
        4: "Project_26041_bridge.csv",
        5: "linking_file_200K_withbatch.csv",
        6: "linking_file_300K_withbatch.csv",
    }
    return f"gs://broad-ukbb/resources/array/{array_map_names[freeze]}"


def array_sample_map_ht_path(
    data_source: str, freeze: int = CURRENT_FREEZE,
) -> hl.Table:
    """
    Returns array sample map Table.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: array sample map Table
    :rtype: hl.Table
    """
    if freeze <= 5:
        f"gs://broad-ukbb/{data_source}.freeze_{freeze}/array_sample_map.ht"
    return f"gs://broad-ukbb/{data_source}.freeze_{freeze}/array/array_sample_map.ht"


# Capture intervals
def capture_ht_path(data_source: str) -> str:
    """
    Returns path to capture intervals Table.

    :param str data_source: One of 'regeneron' or 'broad'
    :return: Path to capture intervals Table
    :rtype: str
    """
    if data_source == "broad":
        return "gs://broad-ukbb/resources/xgen_plus_spikein.Homo_sapiens_assembly38.targets.pad50.merged.ht"
    elif data_source == "regeneron":
        return "gs://broad-ukbb/resources/ukbb_exome_calling_intervals.summary.ht"
    else:
        raise DataException("This data_source is currently not present")


# Phenotype resources
ukbb_phenotype_path = "gs://broad-ukbb/resources/ukb24295.phenotypes.txt"


def phenotype_ht_path() -> str:
    """
    Returns path to phenotype Table.

    :return: Path to phenotype Table
    :rtype: str
    """
    return "gs://broad-ukbb/resources/ukb24295.phenotypes.ht"


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
    return f"gs://broad-ukbb/{data_source}.freeze_{freeze}/release"


def release_mt_path(data_source: str, freeze: int, temp=False) -> str:
    """
    Fetch filepath for release Hail MatrixTables

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param bool temp: If True, fetch Table in which nested variant annotations are unfurled but listed under 'info' rather
        than at the top level; created as a checkpoint during VCF generation
    :return: Filepath for release MatrixTable
    :rtype: str
    """
    tag = ".temp" if temp else ""
    return f"{get_release_path(data_source, freeze)}/mt/{data_source}.freeze_{freeze}.nested{tag}.mt"


def release_ht_path(data_source: str, freeze: int) -> str:
    """
    Fetch filepath for release Hail Tables

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Filepath for release Table
    :rtype: str
    """
    return f"{get_release_path(data_source, freeze)}/ht/{data_source}.freeze_{freeze}.release.sites.ht"


def release_var_hist_path(data_source: str, freeze: int) -> str:
    """
    Fetch bucket for release variant histograms (json files)

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Filepath for release jsons
    :rtype: str
    """
    return f"{get_release_path(data_source, freeze)}/json/{data_source}.freeze_{freeze}.json"


def release_vcf_path(data_source: str, freeze: int, contig: str = None) -> str:
    """
    Fetch bucket for release (variant-only) VCFs

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param str contig: String containing the name of the desired reference contig
    :return: Filepath for the desired VCF
    :rtype: str
    """
    if contig:
        return f"{get_release_path(data_source, freeze)}/vcf/{data_source}.freeze_{freeze}.{contig}.vcf.bgz"
    else:
        # if contig is None, return path to sharded vcf bucket
        return f"{get_release_path(data_source, freeze)}/vcf/sharded_vcf/{data_source}.freeze_{freeze}"
