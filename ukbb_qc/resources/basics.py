import logging
from typing import Optional

import hail as hl

from gnomad.utils.file_utils import file_exists
from gnomad.resources.resource_utils import DataException

from .resource_utils import CURRENT_FREEZE, DATA_SOURCES, FREEZES
from .sample_qc import meta_ht_path


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("basics")
logger.setLevel(logging.INFO)


# UKBB data resources
def pan_ancestry_txt_path() -> str:
    """
    Returns path to text file that contains pan-ancestry label per sample.

    :return: Path to file that contains pan-ancestry labels.
    :rtype: str
    """
    return "gs://broad-ukbb/resources/pan_ancestry.txt"


def pan_ancestry_ht_path() -> str:
    """
    Returns path to Table that contains pan-ancestry label per sample.

    :return: Path to HT that contains pan-ancestry labels.
    :rtype: str
    """
    return "gs://broad-ukbb/resources/pan_ancestry.ht"


def pan_ancestry_bridge_path() -> str:
    """
    Returns path to text file that contains sample ID mappings for pan-ancestry labels.

    :return: Path to text file that contains pan-ancestry sample ID mappings.
    :rtype: str
    """
    return "gs://broad-ukbb/resources/bridge_26041_31063.csv"


def excluded_samples_path(freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to list of samples to exclude from QC due to withdrawn consents

    :param int freeze: One of data freezes.
    :return: Path to excluded samples list
    :rtype: str
    """
    # NOTE: we did not have files with withdrawn sample IDs for freezes 4 and 5
    excluded_file_names = {
        6: "w26041_20200204.csv",
        7: "w26041_20200820.csv",
    }
    return f"gs://broad-ukbb/resources/withdrawn_consents/{excluded_file_names[freeze]}"


def dup_resolution_path(freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to bucket containing files necessary to resolve duplicate samples.

    Currently only exists for freeze 7/450k callset.
    :param int freeze: One of data freezes
    :return: Path to bucket with duplicate sample resolution files
    :rtype: str
    """
    if freeze != 7:
        raise DataException(
            "Duplicate resolution bucket only exists for freeze 7/450k!"
        )
    return f"gs://broad-ukbb/broad.freeze_{freeze}/duplicate_resolution"


def dup_gvcf_path(freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to file containing duplicate samples and their most recent gVCFs.

    Currently only exists for the 450k MatrixTable.

    :param int freeze: One of data freezes
    :return: Path to duplicates sample gVCF path tsv file
    :rtype: str
    """
    if freeze != 7:
        raise DataException("Duplicate gVCF path only exists for freeze 7/450k!")
    return f"{dup_resolution_path(freeze)}/duplicate_sample_map_no_UU.tsv"


def dup_mt_path(freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to MT created using most recent gVCF version for duplicate samples.

    Currently only exists for the 450k MatrixTable.

    MT was created using the hail gVCF combiner on 12/16/2020

    :param int freeze: One of data freezes
    :return: Path to duplicate sample MT
    :rtype: str
    """
    if freeze != 7:
        raise DataException("Duplicate MT only exists for freeze 7/450k!")
    return f"{dup_resolution_path(freeze)}/most_recent_dup.mt"


def dup_map_path(freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to TSV file containing duplicate sample IDs and undesired column index.

    Currently only exists for the 450k MatrixTable.

    TSV has two columns: sample ID of duplicate sample and column index to remove from the 450 MT

    :param int freeze: One of data freezes.
    :return: Path to duplicate sample mapping TSV
    :rtype: str
    """
    if freeze != 7:
        raise DataException("Duplicate map file only exists for freeze 7/450K!")
    return f"{dup_resolution_path(freeze)}/dup_remove_idx.tsv"


def known_dups_ht_path(freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table containing samples known to be on the pharma partners' remove lists.

    All 44 samples are marked as "unresolved duplicates" by the pharma partners.

    :param int freeze: One of data freezes.
    :return: Path to known duplicates HT.
    :rtype: str
    """
    if freeze != 7:
        raise DataException("Known duplicates file only exists for freeze 7/450K!")
    return f"gs://broad-ukbb/broad.freeze_{freeze}/duplicate_resolution/pharma_known_dups.ht"


def check_dups_to_remove(
    remove_ids: hl.expr.ArrayExpression, samples: hl.Table, column_name: str = "s",
) -> int:
    """
    This function checks whether the number of duplicate samples in the input Table matches expected counts.

    Used to check numbers prior to duplicate sample removal. 
    Function was written specifically to resolve duplicates in freeze 7/the 450k callset.

    :param hl.expr.ArrayExpression remove_list: ArrayExpression containing list of sample IDs to remove.
    :param hl.Table samples: Table containing all sample IDs.
    :param str column_name: Name of column containing sample IDs in Table. Default is 's'.
    :return: Number of samples to remove from Table.
    :rtype: int
    """
    # Using an HT here because aggregate_cols has been slow/memory intensive in the past
    n_samples_to_drop = samples.aggregate(
        hl.agg.count_where(remove_ids.contains(samples[column_name]))
    )
    if n_samples_to_drop != hl.eval(hl.len(remove_ids)):
        raise DataException(
            f"Expecting to remove {hl.eval(hl.len(remove_ids))} duplicate samples but found {n_samples_to_drop}. Double check samples in MT!"
        )
    return n_samples_to_drop


def get_ukbb_data(
    data_source: str,
    freeze: int = CURRENT_FREEZE,
    key_by_locus_and_alleles: bool = False,
    adj: bool = False,
    split: bool = True,
    ukbb_samples_only: bool = False,
    raw: bool = False,
    repartition: bool = False,
    n_partitions: int = 30000,
    meta_root: Optional[str] = None,
) -> hl.MatrixTable:
    """
    Wrapper function to get UKBB data as MatrixTable. By default, returns split hardcalls (with adj annotated but not filtered).

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool key_by_locus_and_alleles: Whether to key the MatrixTable by locus and alleles
    :param bool adj: Whether the returned data should be filtered to adj genotypes
    :param bool split: Whether the dataset should be split (only applies to raw=False)
    :param bool ukbb_samples_only: Whether to return only UKBB samples (exclude control samples). Default is False.
        Relevant only when running release code.
    :param bool raw: Whether to return the raw data (not recommended: unsplit, and no special consideration on sex chromosomes)
    :param bool repartition: Whether to repartition the MatrixTable. 
        Required if raw is True for tranche 3/freeze 6. Default is False
    :param int n_partitions: Number of desired partitions for MatrixTable. 
        Required if raw is True for tranche 3/freeze 6. Default is 30000
    :param str meta_root: Root annotation name for metadata (e.g., 'meta')
    :return: hardcalls dataset
    :rtype: MatrixTable
    """
    from gnomad.utils.filtering import filter_to_adj

    if raw and split:
        raise DataException("No split raw data. Use of hardcalls is recommended.")
    if not raw and not split:
        raise DataException(
            "No unsplit hardcalls. Use of split hardcalls is recommended."
        )

    # Check for array sample map -- this is loaded first in the pipeline, before the MT is loaded/processed
    if not file_exists(f"{array_sample_map_ht_path(freeze)}"):
        raise DataException(f"Need to import array sample map ht for freeze {freeze}!")

    if repartition:
        mt = hl.read_matrix_table(
            path=get_ukbb_data_path(data_source, freeze, hardcalls=not raw),
            _n_partitions=n_partitions,
        )
    else:
        mt = hl.read_matrix_table(
            get_ukbb_data_path(data_source, freeze, hardcalls=not raw)
        )
    # Add column index to samples
    mt = mt.add_col_index()

    if adj:
        mt = filter_to_adj(mt)

    if meta_root:
        meta_ht = hl.read_table(meta_ht_path(data_source, freeze))
        mt = mt.annotate_cols(**{meta_root: meta_ht[mt.s]})

    if key_by_locus_and_alleles:
        mt = hl.MatrixTable(
            hl.ir.MatrixKeyRowsBy(mt._mir, ["locus", "alleles"], is_sorted=True)
        )  # taken from v3 qc
    logger.info(f"Number of samples in MT: {mt.count_cols()}")

    # Add warning that no samples will be removed if excluded samples file doesn't exist
    if not file_exists(excluded_samples_path(freeze)):
        logger.warning(
            "No excluded samples file found. No samples will be removed from MT"
        )

    if file_exists(excluded_samples_path(freeze)) or ukbb_samples_only:

        # Read in array sample map HT
        sample_map_ht = hl.read_table(array_sample_map_ht_path(freeze))

        # Filter to UKBB samples if specified
        if ukbb_samples_only:
            mt = mt.filter_cols(hl.is_defined(sample_map_ht[mt.col_key]))

        # Remove any samples with withdrawn consents
        # NOTE: Keeping samples with missing consents to avoid filtering any samples present in MT but not in sample map HT
        if file_exists(excluded_samples_path(freeze)):
            mt = mt.filter_cols(
                (~sample_map_ht[mt.col_key].withdrawn_consent)
                | (hl.is_missing(sample_map_ht[mt.col_key].withdrawn_consent))
            )

            # Double check all withdrawn samples were actually excluded
            withdrawn_ht = hl.import_table(
                excluded_samples_path(freeze), no_header=True,
            ).key_by("f0")
            mt_samples = mt.annotate_cols(
                ukbb_app_26041_id=sample_map_ht[mt.col_key].ukbb_app_26041_id
            ).cols()
            mt_samples = mt_samples.key_by("ukbb_app_26041_id")
            withdrawn_samples_in_mt = mt_samples.filter(
                hl.is_defined(withdrawn_ht[mt_samples.ukbb_app_26041_id])
            ).count()

            if withdrawn_samples_in_mt > 0:
                raise DataException(
                    "Withdrawn samples present in MT. Double check sample filtration"
                )
            else:
                logger.info("No withdrawn samples found in MT")

    # Code to resolve duplicate samples specifically in freeze 7/the 450k callset
    if freeze == 7:
        # Remove fully duplicated IDs when reading in raw MT only
        # These 27 duplicates are discussed in this ticket:
        # https://broadinstitute.atlassian.net/browse/PO-28339
        if raw:
            logger.info("Resolving duplicate sample IDs in the 450k MT...")

            # Create list of sample IDs to remove
            remove_ids = []
            with hl.hadoop_open(dup_map_path(freeze), "r") as d:
                for line in d:
                    line = line.strip().split("\t")
                    remove_ids.append(f"{line[0]}_{line[1]}")

            # Remove sample IDs that are present in remove_ids list
            remove_ids = hl.literal(remove_ids)
            mt = mt.annotate_cols(new_s=hl.format("%s_%s", mt.s, mt.col_idx))
            # Check number of samples to remove -- should be 27 here
            samples_to_drop = check_dups_to_remove(remove_ids, mt.cols(), "new_s")
            logger.info(f"Removing {samples_to_drop} samples...")
            mt = mt.filter_cols(~remove_ids.contains(mt.new_s)).drop("new_s", "col_idx")

        # Remove samples that are known to be on the pharma's sample remove list
        # These samples were discovered when running array concordance
        # and were resolved in this notebook:
        # gs://broad-ukbb/broad.freeze_7/notebooks/array_dups.ipynb
        # All 44 of these samples have the note to "remove ID1"
        # Example:
        # "UKB_1534238_0230755066"  "UKB_1534238_236168062-DUP" "Known duplication; blacklist ID1"
        # NOTE: UKB_4316697_236168099-DUP is in the 450K MT
        # but the othersample ID known to be on the remove list, UKB_3489701_236190071,
        # is not present, so this pair is NOT in the known dups HT
        # Therefore, the number of samples to remove in this check is 43
        # "UKB_3489701_236190071"   "UKB_4316697_236168099-DUP" "Known duplication; blacklist ID1"
        dups_ht = hl.read_table(known_dups_ht_path(freeze))
        ids_to_remove = dups_ht.aggregate(
            hl.agg.collect(dups_ht["Sample Name - ID1"]), _localize=False
        )

        # Check number of samples to remove -- should be 44 here
        num_ids = check_dups_to_remove(ids_to_remove, mt.cols())
        logger.info(f"Removing {num_ids} samples...")
        mt = mt.filter_cols(~ids_to_remove.contains(mt.s))

        logger.info(f"Sample count post-filtration: {mt.count_cols()}")

    gt_expr = mt.GT if split else mt.LGT
    mt = mt.filter_rows(hl.agg.any(gt_expr.is_non_ref() | hl.is_defined(mt.END)))
    return mt


def get_ukbb_data_path(
    data_source: str, freeze: int = CURRENT_FREEZE, hardcalls: bool = False,
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
    is_temp: bool = False,
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

    if densified and freeze != 5:
        raise DataException("Densified parameter only applies to freeze 5")

    if data_source == "regeneron":
        return f"gs://broad-ukbb/{data_source}.freeze_{freeze}/data/{data_source}.freeze_{freeze}.nf.mt"
    elif data_source == "broad":
        # create dict for raw mts for freeze 5 and later
        dsp_prefix = "gs://broad-pharma5-ukbb-outputs"
        raw_mt_names = {
            5: "hail_dataproc_20191108115937",
            6: "hail_dataproc_20200130092005.mt",
            7: "hail_450k_dataproc_20201207164032.mt",
        }
        if freeze == 4 or (freeze == 5 and densified):
            return f"gs://broad-ukbb/{data_source}.freeze_{freeze}/data/{data_source}.freeze_{freeze}{tempstr}.mt"
        else:
            return f"{dsp_prefix}/{raw_mt_names[freeze]}"


def hardcalls_mt_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to hardcalls MatrixTable.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
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
def last_END_positions_ht_path(freeze: int = CURRENT_FREEZE) -> str:
    """
    Retrieves path to Table with last END positions annotation

    :param str data_type: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to Table with last_END_position annotation
    :rtype: str
    """
    data_source = "broad"
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
        7: "linking.file.455Kexomes.Broad.csv",
    }
    return f"gs://broad-ukbb/resources/array/{array_map_names[freeze]}"


def array_sample_map_ht_path(freeze: int = CURRENT_FREEZE) -> hl.Table:
    """
    Returns array sample map Table.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: array sample map Table
    :rtype: hl.Table
    """
    return f"gs://broad-ukbb/resources/array/array_sample_map_freeze_{freeze}.ht"


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
def ukbb_phenotype_path() -> str:
    """
    Returns path to phenotype text file.

    :return: Path to phenotypes.txt
    :rtype: str
    """
    return "gs://broad-ukbb/resources/ukb24295.phenotypes.txt"


def phenotype_ht_path() -> str:
    """
    Returns path to phenotype Table.

    :return: Path to phenotype Table
    :rtype: str
    """
    return "gs://broad-ukbb/resources/ukb24295.phenotypes.ht"


def geographical_ht_path() -> str:
    """
    Returns path to location Table.

    :return: Path to location Table
    :rtype: str
    """
    return "gs://broad-ukbb/resources/geographical_data.ht"


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


def release_mt_path(data_source: str, freeze: int) -> str:
    """
    Fetch filepath for release Hail MatrixTable.

    .. note::
        Annotations in this MT are not nested into individual structs.
        Instead, variant annotations are unfurled and listed under the info struct (in preparation for VCF export).

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Filepath for release MatrixTable
    :rtype: str
    """
    return f"gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/{data_source}.freeze_{freeze}.mt"


def release_ht_path(data_source: str, freeze: int) -> str:
    """
    Fetch filepath for release Hail Tables

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Filepath for release Table
    :rtype: str
    """
    return f"{get_release_path(data_source, freeze)}/ht/{data_source}.freeze_{freeze}.release.sites.ht"


def release_vcf_ht_path(data_source: str, freeze: int) -> str:
    """
    Fetch filepath for release Hail Table with annotations reformatted for VCF export

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes. Resource only exists for freeze 7.
    :return: Filepath for release VCF Table
    :rtype: str
    """
    if freeze != 7:
        raise DataException("Release VCF HT only exists for freeze 7/455K!")
    return f"{get_release_path(data_source, freeze)}/ht/{data_source}.freeze_{freeze}.release.sites.vcf.ht"


def release_var_hist_path(data_source: str, freeze: int) -> str:
    """
    Fetch bucket for release variant histograms (json files)

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Filepath for release jsons
    :rtype: str
    """
    return f"{get_release_path(data_source, freeze)}/json/{data_source}.freeze_{freeze}.json"


def release_header_path(data_source: str, freeze: int) -> str:
    """
    Fetch path to pickle file containing VCF header dictionary.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Filepath for header dictionary pickle
    :rtype: str
    """
    return f"gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/header_dict.pickle"


def release_vcf_path(data_source: str, freeze: int, contig: str = None) -> str:
    """
    Fetch bucket for release (variant-only) VCFs

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param str contig: String containing the name of the desired reference contig
    :return: Filepath for the desired VCF
    :rtype: str
    """
    if contig:
        return f"{get_release_path(data_source, freeze)}/vcf/{data_source}.freeze_{freeze}.{contig}.vcf.bgz"
    else:
        # if contig is None, return path to sharded vcf bucket
        # NOTE: need to add .bgz or else hail will not bgzip shards
        return f"{get_release_path(data_source, freeze)}/vcf/sharded_vcf/{data_source}.freeze_{freeze}.bgz"


def append_to_vcf_header_path(data_source: str, freeze: int) -> str:
    """
    Fetch path to TSV file containing extra fields to append to VCF header.

    Extra fields are VEP and dbSNP versions.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Filepath for extra fields TSV file
    :rtype: str
    """
    if freeze != 7:
        raise DataException(
            "Extra fields to append to VCF header TSV only exists for freeze 7/455K!"
        )
    return (
        f"gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/append_to_vcf_header.tsv"
    )


def annotation_hists_path(data_source: str, freeze: int) -> str:
    """
    Returns path to file containing ANNOTATIONS_HISTS dictionary.

    Dictionary contains histogram values for each metric. 
    For example, "InbreedingCoeff": [-0.25, 0.25, 50].

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Path to file with annotations histograms
    :rtype: str
    """
    return f"{get_release_path(data_source, freeze)}/annotation_hists.json"


def vqsr_sites_path(
    data_source: str, freeze: int = CURRENT_FREEZE, ht: bool = True
) -> str:
    """
    Fetch bucket for sites HT or VCF (input to VQSR)

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param bool ht: Will return path to Table if True, otherwise will return path to sharded VCF directory. Default is True
        Sharded VCF directory needs to have '.bgz' extension or the VCF shards it contains will NOT be bgzipped.
    :return: Path to VQSR sites Table or VCFs
    :rtype: str
    """
    suffix = "ht" if ht else "vcf.bgz"
    return f"{get_release_path(data_source, freeze)}/vqsr/{data_source}.freeze_{freeze}.sites_for_vqsr.{suffix}"


def release_summary_ht_path(
    data_source: str, freeze: int = CURRENT_FREEZE, intervals: bool = False
) -> str:
    """
    Returns path to HT with summary counts per variant category.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param bool intervals: Whether HT has been filtered to interval QC pass regions only. Default is False.
    :return: Path to summary counts HT
    :rtype: str
    """
    if freeze not in (6, 7):
        raise DataException(
            "Release summary HT path only exists for freezes 6/300K and 7/450K!"
        )

    return f"{get_release_path(data_source, freeze)}/summary/summary_per_variant{'_intervals' if intervals else ''}.ht"


def release_lof_mt_path(
    data_source: str, freeze: int = CURRENT_FREEZE, intervals: bool = False
) -> str:
    """
    Returns path to summary MT containing gene LoF matrix.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param bool intervals: Whether MT has been filtered to interval QC pass regions only. Default is False.
    :return: Path to summary MT with gene LoF matrix
    :rtype: str
    """
    if freeze not in (6, 7):
        raise DataException(
            "Release LoF MT path only exists for freezes 6/300K and 7/450K!"
        )

    return f"{get_release_path(data_source, freeze)}/summary/gene_lof_matrix{'_intervals' if intervals else ''}.mt"


def release_lof_ht_path(
    data_source: str, freeze: int = CURRENT_FREEZE, intervals: bool = False
) -> str:
    """
    Returns path to summary HT containing gene LoF matrix summary.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param bool intervals: Whether MT has been filtered to interval QC pass regions only. Default is False.
    :return: Path to summary HT with gene LoF matrix
    :rtype: str
    """
    if freeze not in (6, 7):
        raise DataException(
            "Release LoF HT path only exists for freezes 6/300K and 7/450K!"
        )

    return f"{get_release_path(data_source, freeze)}/summary/gene_lof_matrix_summary{'_intervals' if intervals else ''}.ht"


def get_doubleton_ht_path(
    data_source: str, freeze: int = CURRENT_FREEZE, unrelated_only: bool = False
) -> str:
    """
    Returns path to HT containing doubletons and relevant sample IDs.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param bool unrelated_only: Whether the HT contains only unrelated samples. Default is False.
    :return: Path to doubletons HT
    :rtype: str
    """
    if freeze not in (6, 7):
        raise DataException(
            "Doubleton HT path only exists for freezes 6/300K and 7/450K!"
        )

    return f"{get_release_path(data_source, freeze)}/summary/doubleton/doubletons{'_unrelated' if unrelated_only else ''}.ht"


def get_pair_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to HT containing randomly selected sample pairs.

    Used to compare to pairs in doubleton HT.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Path to comparison pairs HT
    :rtype: str
    """
    if freeze not in (6, 7):
        raise DataException(
            "Comparison pairs HT path only exists for freezes 6/300K and 7/450K!"
        )

    return (
        f"{get_release_path(data_source, freeze)}/summary/doubleton/comparison_pairs.ht"
    )


# Readviz resources
def cram_map_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Return path to TSV containing UKBB samples and their cram paths.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: Data freeze. Must be 7.
    :return: Path to TSV containing UKBB samples and their cram paths
    :rtype: str  
    """
    if freeze != 7:
        raise DataException("UKBB cram map path only exists for freeze 7/450K!")
    return "gs://gnomad-readviz/ukbb/sample_cram_map.tsv"


def non_gnomad_var_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Return path to HT containing variants unique to the UKBB (not present in gnomAD v2.1.1 or v3.1).

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: Data freeze. Must be 7.
    :return: Path to UKBB unique variants HT
    :rtype: str   
    """
    if freeze != 7:
        raise DataException(
            "UKBB unique variants (not in gnomAD) HT path only exists for freeze 7/450K!"
        )
    return "gs://gnomad-readviz/ukbb/not_in_gnomad_variants/ukbb_unique_variants.ht"


def readviz_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to HT containing variants unique to the UKBB (not present in gnomAD v2.1.1 or v3.1).

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: Data freeze. Must be 7.
    :return: Path to UKBB unique variants HT
    :rtype: str   
    """
    if freeze != 7:
        raise DataException(
            "UKBB unique variants (not in gnomAD) HT path only exists for freeze 7/450K!"
        )
    return "gs://gnomad-readviz/ukbb/not_in_gnomad_variants/samples_for_readviz.ht"


# logging path
def logging_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to bucket that stores hail logs

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Path to bucket with hail logs
    :rtype: str
    """
    return f"gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/logs/"
