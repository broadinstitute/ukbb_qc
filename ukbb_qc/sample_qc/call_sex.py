import logging
from gnomad_hail.utils.sample_qc import get_ploidy_cutoffs, get_sex_expr
from gnomad_hail.utils.generic import get_reference_genome
from ukbb_qc.resources.basics import capture_ht_path
from ukbb_qc.resources.resource_utils import *
from ukbb_qc.resources.sample_qc import ploidy_ht_path, sex_ht_path


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("call_sex")
logger.setLevel(logging.INFO)


def impute_sex(mt: hl.MatrixTable, ploidy_ht: hl.Table, aaf_threshold = 0.001, male_threshold: float = 0.75, female_threshold: float = 0.5) -> hl.Table:
    """
    Imputes sex and annotates matrix table with results, outputs a histogram of fstat values
    NOTE: decide on thresholds for male/female by looking at fstat plots

    :param MatrixTable mt: MatrixTable containing samples to be ascertained for sex
    :param Table ploidy_ht: Table with sex chromosome ploidies
    :param str aaf_threshold: Minimum alternate allele frequency
    :param str male_threshold: Threshold above which a sample will be called male
    :param str female_threshold: Threshold below which a sample will be called female
    :return: Table with imputed sex annotations stashed in column annotation 'sex_karyotype'
    :rtype: Table
    """
    logger.info("Filtering to chrX...")
    mt = hl.filter_intervals(mt, 
                            [hl.parse_locus_interval(x_contig, reference_genome=get_reference_genome(mt.locus).name) 
                            for x_contig in get_reference_genome(mt.locus).x_contigs]) 

    # TODO update this for sparse
    sex_ht = hl.impute_sex(mt.GT, aaf_threshold=aaf_threshold, male_threshold=male_threshold, female_threshold=female_threshold)
    sex_ht = sex_ht.annotate(**ploidy_ht[sex_ht.s])
    x_ploidy_cutoffs, y_ploidy_cutoffs = get_ploidy_cutoffs(sex_ht, 0.5)
    return sex_ht.annotate(
            **get_sex_expr(
                sex_ht.chrX_ploidy,
                sex_ht.chrY_ploidy,
                x_ploidy_cutoffs,
                y_ploidy_cutoffs
        )
    )


def run_impute_ploidy(mt: hl.MatrixTable, data_source: str, freeze: int) -> hl.Table:
    """
    Detects sex chromosome ploidy for samples in given MatrixTable and writes ploidy Table

    :param MatrixTable mt: Input MatrixTable 
    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Table with imputed sex chromosome ploidies
    :rtype: Table
    """
    logger.info("Imputing sex ploidy...")
    ploidy_ht = impute_sex_ploidy(mt, included_intervals=capture_ht_path(data_source, freeze))
    ploidy_ht = ploidy_ht.checkpoint(ploidy_ht_path(data_source, freeze), overwrite=True)
    return ploidy_ht


def run_impute_sex(mt: hl.MatrixTable, data_source: str, freeze: int) -> hl.Table:
    """
    Imputes sex karyotypes for samples in given MatrixTable and writes Table

    :param MatrixTable mt: Input MatrixTable 
    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Table with imputed sex karyotypes
    :rtype: Table
    """
    # NOTE: Broad callset does not have filter annotations
    if data_source == "regeneron":
        mt = mt.filter_rows(hl.is_missing(mt.filters))

    logger.info("Filtering to biallelic SNPs...")
    mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]))

    logger.info("Imputing sex...")
    sex_ht = impute_sex(mt, ploidy_ht)
    sex_ht = sex_ht.checkpoint(sex_ht_path(data_source, freeze), overwrite=True)
