import logging
import hail as hl

from gnomad.resources.resource_utils import DataException
from .resource_utils import CURRENT_FREEZE, DATA_SOURCES, FREEZES


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sample_qc_resources")
logger.setLevel(logging.INFO)


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
    return f"gs://broad-ukbb/{data_source}.freeze_{freeze}/sample_qc"


# Interval resources
def interval_qc_path(
    data_source: str, freeze: int = CURRENT_FREEZE, chrom: str = None, ht: bool = True
) -> str:
    """
    Returns path to interval QC results (either as a Table or text file).

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param str chrom: Chromosome of interest, 'autosomes', or 'sex_chr'
    :param bool ht: Whether to return path to Table. If False, will return a text file. Default is True
    :return: Path to interval QC results
    :rtype: str
    """
    if chrom is None:
        chrom = ""
    else:
        chrom = f".{chrom}"

    return f'{sample_qc_path(data_source, freeze)}/interval_qc/coverage_by_target{chrom}{".ht" if ht else ".txt"}'


def f_stat_sites_path() -> str:
    """
    Returns path to Table with high callrate, common, pass interval QC, biallelic SNP positions on chromosome X used in sex imputation.

    NOTE: The sites came from tranche 2/freeze 5 (first 200K samples and the last dataset with AF).

    :return: Path to Table with sites for sex imputation
    :rtype: str
    """
    data_source = "broad"
    freeze = 5
    return f"{sample_qc_path(data_source, freeze)}/f_stat_sites.ht"


def qc_sites_path() -> str:
    """
    Returns path to Table with sites to use in QC MatrixTable generation.

    NOTE: Sites were generated using QC MT from tranche 2/freeze 5. 
    Chose high callrate (>0.99), high coverage (85% samples > 20x), and common (AF > 0.001) bi-allelic SNPs outside of LCR intervals.

    :return: Path to QC MT sites Table
    :rtype: str
    """
    data_source = "broad"
    freeze = 5
    return f"{sample_qc_path(data_source, freeze)}/qc_data/qc_sites.ht"


# QC resources (meta ht, qc mt, qc ht, qc sites ht)
def meta_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to metadata Table

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to metadata Table
    :rtype: str
    """
    return f"{sample_qc_path(data_source, freeze)}/meta.ht"


def qc_mt_path(
    data_source: str, freeze: int = CURRENT_FREEZE, ld_pruned: bool = True
) -> str:
    """
    Returns path to MatrixTable filtered to high callrate, common, biallelc snps for sample QC purposes

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool ld_pruned: Whether to return the LD pruned version of the MatrixTable
    :return: Path MatrixTable for sample QC purposes
    :rtype: str
    """
    if not ld_pruned and freeze > 5:
        raise DataException(
            "Only one version of QC MT stored for tranches after tranche 2 (created using ld pruned tranche 2 sites)"
        )
    ld_pruned = ".pruned" if ld_pruned else ""
    return f"{sample_qc_path(data_source, freeze)}/qc_data/high_callrate_common_biallelic_snps{ld_pruned}.mt"


def get_qc_mt(
    data_source: str,
    freeze: int = CURRENT_FREEZE,
    ld_pruned: bool = True,
    filter_lowqual: bool = True,
) -> hl.MatrixTable:
    """
    Returns MatrixTable filtered to high callrate, common, biallelc snps for sample QC purposes

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool ld_pruned: Whether to return the LD pruned version of the MatrixTable. Default is True
    :param bool filter_lowqual: Whether to filter low QUAL variants immediately after reading in MatrixTable. Default is True
    :return: MatrixTable for sample QC purposes
    :rtype: hl.MatrixTable
    """
    mt = hl.read_matrix_table(qc_mt_path(data_source, freeze, ld_pruned))

    if filter_lowqual:
        logger.info("Removing low QUAL variants...")
        mt = mt.filter_rows(~mt.lowqual)

    return mt


def qc_ht_path(
    data_source: str, freeze: int = CURRENT_FREEZE, ld_pruned: bool = True
) -> str:
    """
    Returns path to Table for sample QC purposes (specifically, relatedness checks).

    This is a Table keyed by sample with metrics from hail's sample QC.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to Table for sample QC purposes
    :rtype: str
    """
    ld_pruned = ".pruned" if ld_pruned else ""
    return f"{sample_qc_path(data_source, freeze)}/qc_data/high_callrate_common_biallelic_snps{ld_pruned}.ht"


# Hard filters resources (including sex imputation resources)
def ploidy_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with sex chromosome ploidies

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to ploidy Table
    :rtype: str
    """
    if data_source == "broad" and freeze >= 5:
        return f"{sample_qc_path(data_source, freeze)}/sex_check/ploidy.ht"
    else:
        raise DataException("No ploidy file specified for this data source and freeze")


def sex_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with inferred sex karyotypes

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to sex Table
    :rtype: str
    """
    return f"{sample_qc_path(data_source, freeze)}/sex_check/sex.ht"


def hard_filters_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with hard filter status

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to hard filters Table
    :rtype: str
    """
    return f"{sample_qc_path(data_source, freeze)}/hard_filters_flagged.ht"


# Array resources
def array_concordance_sites_path() -> str:
    """
    Returns path to Table with sites to use when calculating array concordance

    All sites are autosomal sites from 200K (tranche 2/freeze 5)
    Chose high callrate (>0.95), high coverage (>85% samples at 20x), common (AF > 0.0001) sites

    :return: Path to sites Table
    :rtype: Table
    """
    return "gs://broad-ukbb/broad.freeze_5/temp/sites_for_array_concordance.ht"


def array_concordance_results_path(
    data_source: str, freeze: int = CURRENT_FREEZE, sample: bool = True,
) -> str:
    """
    Returns path to Table with array sample concordance status

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool sample: Whether to return sample concordance results. If False, returns variant results. Default is True.
    :return: Path to Table with array concordance results (sample or variant)
    :rtype: str
    """
    return f"{sample_qc_path(data_source, freeze)}/array_concordance/{'sample' if sample else 'variant'}_concordance.ht"


def get_ukbb_array_pcs_path() -> str:
    """
    Returns path to text file with UKBB array PCs

    :return: Path to text file with PCs
    :rtype: str
    """
    return "gs://broad-ukbb/resources/array/ukb_sqc_v2_ukb26041.txt"


def get_ukbb_array_pcs_ht_path() -> str:
    """
    Returns path to Table with UKBB array PCs

    :return: Path to Table with PCs
    :rtype: str
    """
    return "gs://broad-ukbb/resources/array/ukb24295.pcs.ht"


# Platform PCA resources
def callrate_mt_path(
    data_source: str, freeze: int = CURRENT_FREEZE, interval_filtered: bool = False
) -> str:
    """
    Returns path to callrate MatrixTable 

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool interval_filtered: Whether data was filtered to high coverage intervals. Default is False
    :return: Path to callrate MatrixTable
    :rtype: str
    """
    filtered = ".interval_filtered" if interval_filtered else ""
    return f"{sample_qc_path(data_source, freeze)}/platform_pca/callrate{filtered}.mt"


def platform_pca_scores_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with platform PCA scores

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to platform PCA scores Table
    :rtype: str
    """
    return f"{sample_qc_path(data_source, freeze)}/platform_pca/platform_pca_scores.ht"


def platform_pca_loadings_ht_path(
    data_source: str, freeze: int = CURRENT_FREEZE
) -> str:
    """
    Returns path to Table with platform PCA loadings

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to platform PCA loadings Table
    :rtype: str
    """
    return (
        f"{sample_qc_path(data_source, freeze)}/platform_pca/platform_pca_loadings.ht"
    )


def platform_pca_assignments_ht_path(
    data_source: str, freeze: int = CURRENT_FREEZE
) -> str:
    """
    Returns path to Table with platform PCA results

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to platform PCA results Table
    :rtype: str
    """
    # NOTE: for tranche 1 and 2 (freeze 4 and 5; 100K and 200K), this table was named "platform_pca_results"
    return f"{sample_qc_path(data_source, freeze)}/platform_pca/platform_pca_assignments.ht"


# Relatedness inference resources
def relatedness_pca_scores_ht_path(
    data_source: str, freeze: int = CURRENT_FREEZE
) -> str:
    """
    Returns path to Table with relatedness PCA scores (run prior to PC-relate)

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to relatedness PCA scores Table
    :rtype: str
    """
    return f"{sample_qc_path(data_source, freeze)}/relatedness/pruned.pca_scores.ht"


def relatedness_ht_path(
    data_source: str, freeze: int = CURRENT_FREEZE, method: str = "pc_relate"
) -> str:
    """
    Returns path to Table output by PC-relate with inferred relationships

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param str method: Method of inferring relatedness. Default: hail's pc_relate
    :return: Path to relatedness Table
    :rtype: str
    """

    method = f".{method}" if method != "pc_relate" else ""
    return f"{sample_qc_path(data_source, freeze)}/relatedness/relatedness{method}.ht"


def duplicates_ht_path(
    data_source: str,
    freeze: int = CURRENT_FREEZE,
    dup_sets: bool = False,
    method: str = "pc_relate",
) -> str:
    """
    Returns path to Table with inferred duplicates.

    Each row in duplicate sample sets (dup_sets) HT is indexed by the sample that is kept and contains the set of duplicate samples that should be filtered.
    The duplicate HT is created by exploding the dup_sets Table. 
    Each row in the duplicate HT is indexed by sample and has an annotation (`dup_filtered`) indicating which sample from the duplicate pair was kept

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param bool dup_sets: Whether to return path to duplicate set Table. Default is False
    :param str method: Method of inferring relatedness. Default: hail's pc_relate
    :return: Path to duplicates Table
    :rtype: str
    """
    dup_sets = f"_sets" if dup_sets else ""
    method = f".{method}" if method != "pc_relate" else ""
    return f"{sample_qc_path(data_source, freeze)}/relatedness/duplicate{dup_sets}{method}.ht"


def inferred_ped_path(
    data_source: str, freeze: int = CURRENT_FREEZE, method: str = "pc_relate"
) -> str:
    """
    Returns path to Table with inferred families

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param str method: Method of inferring relatedness. Default: hail's pc_relate
    :return: Path to inferred family Table
    :rtype: str
    """

    method = f".{method}" if method != "pc_relate" else ""
    return f"{sample_qc_path(data_source, freeze)}/relatedness/ped{method}.txt"


def related_drop_path(
    data_source: str, freeze: int = CURRENT_FREEZE, method: str = "pc_relate"
) -> str:
    """
    Returns path to Table with related samples to drop

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param str method: Method of inferring relatedness. Default: hail's pc_relate
    :return: Path to related samples to drop Table
    :rtype: str
    """
    method = f".{method}" if method != "pc_relate" else ""
    return f"{sample_qc_path(data_source, freeze)}/relatedness/related_samples_to_drop{method}.ht"


# Population inference resources
def gnomad_ancestry_loadings_liftover_path(checkpoint: bool = False):
    """
    Returns path to gnomAD ancestry loadings lifted to build 38

    :param bool checkpoint: Whether to return path to temporary Table
    :return: Path to gnomAD ancestry loadings Table
    :rtype: str
    """
    if checkpoint:
        return "gs://broad-ukbb/temp/gnomad_joint_unrelated_pca_loadings.ht"
    else:
        return "gs://broad-ukbb/resources/gnomad.joint.unrelated.pca_loadings_lift.ht"


def ancestry_pca_scores_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with ancestry PCA scores from scratch PCA 

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to ancestry PCA scores Table
    :rtype: str
    """
    return f"{sample_qc_path(data_source, freeze)}/population_pca/pca_scores.ht"


def ancestry_pca_loadings_ht_path(
    data_source: str, freeze: int = CURRENT_FREEZE
) -> str:
    """
    Returns path to Table with ancestry PCA loadings from scratch PCA 

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to ancestry PCA loadings Table
    :rtype: str
    """
    return f"{sample_qc_path(data_source, freeze)}/population_pca/pca_loadings.ht"


def ancestry_cluster_ht_path(
    data_source: str, freeze: int = CURRENT_FREEZE, data_type: str = "exome"
) -> str:
    """
    Returns path to Table with ancestry PCA cluster assignments using exome data

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param str data_type: Data type used in ancestry PCA. One of 'exome', 'array', or 'joint'. Default is exome
    :return: Path to ancestry cluster assignments Table
    :rtype: str
    """
    if data_type == "exome":
        name = ""
    elif data_type == "array":
        name = "array_"
    elif data_type == "joint":
        name = "joint_scratch_array_"
    return f"{sample_qc_path(data_source, freeze)}/population_pca/{name}cluster_assignments.ht"


def ancestry_pc_project_scores_ht_path(
    data_source: str, freeze: int = CURRENT_FREEZE, data_type: str = None
) -> str:
    """
    Returns path of Table for scores and pop assignments from pc_project on gnomAD PCs

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param str data_type: either None for UKBB only or joint for merged UKBB and gnomAD
    :return: Path to Table with ancestry PC project scores
    :rtype: str
    """
    data_type = f".{data_type}" if data_type else ""
    return f"{sample_qc_path(data_source, freeze)}/population_pca/pc_project_scores_pop_assign{data_type}.ht"


def ancestry_hybrid_ht_path(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to Table with hybrid scratch and gnomAD PC project ancestry assignments

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to Table with hybrid ancestry assignments
    :rtype: str
    """
    return f"{sample_qc_path(data_source, freeze)}/population_pca/hybrid_pop_assignments.ht"


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
    return f"gs://broad-ukbb/resources/ukb24295.phenotypes.freeze_{freeze}.ht"


# Broad-Regeneron ancestry comparison resources (relevant for freeze 4 only)
def get_regeneron_ancestry_path() -> str:
    """
    Returns path to Regeneron relatedness inference files

    :param int freeze: One of data freezes
    :return: Output relatedness file path
    :rtype: str
    """
    freeze = 4
    freeze_str = "Four"
    return f"gs://broad-ukbb/regeneron.freeze_{freeze}/data/pharma_relatedness_analysis/UKB_Freeze_{freeze_str}.NF.splitmulti.commonsnps_samples_ancestries.txt"


def get_joint_regeneron_ancestry_path(data_source: str) -> str:
    """
    Returns path to Table with both Broad and Regeneron inferred ancestries

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to Table with Broad and Regeneron inferred ancestries
    :rtype: str
    """
    freeze = 4
    return f"{sample_qc_path(data_source, freeze)}/population_pca/regeneron_ukb_joint_ancestry.ht"


# Outlier detection resources
def qc_temp_data_prefix(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to temporary sample QC data. 

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :return: Path to bucket with temporary sample QC data
    :rtype: str
    """
    return f"{sample_qc_path(data_source, freeze)}/temp/"


def platform_pop_outlier_ht_path(
    data_source: str,
    freeze: int = CURRENT_FREEZE,
    pop_assignment_method: str = "hybrid_pop_data",
    platform_assignment_method: str = "batch",
) -> str:
    """
    Returns path to Table containing samples flagged for outlier sample QC metrics

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of data freezes
    :param str pop_assignment_method: Method used to infer populations (gnomad_qc_project_pop_data, hybrid_pop_data). Default is hybrid_pop_data
    :param str platform_assignment_method: Method used to infer platforms (qc_platform or batch). Default is batch
    :return: Path to Table with outlier samples flagged
    :rtype: str
    """
    return f"{sample_qc_path(data_source, freeze)}/outlier_detection/outlier_detection.{pop_assignment_method}_{platform_assignment_method}.ht"
