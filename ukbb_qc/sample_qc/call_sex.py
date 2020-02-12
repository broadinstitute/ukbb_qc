import hail as hl
import logging
from gnomad_hail.utils.sample_qc import default_annotate_sex
from ukbb_qc.resources.basics import capture_ht_path
from ukbb_qc.resources.sample_qc import f_stat_sites_path, ploidy_ht_path, sex_ht_path


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("call_sex")
logger.setLevel(logging.INFO)


def run_impute_ploidy(mt: hl.MatrixTable, data_source: str, freeze: int) -> None:
    """
    Detects sex chromosome ploidy for samples in given MatrixTable and writes ploidy Table

    :param MatrixTable mt: Input MatrixTable 
    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: None
    :rtype: None
    """
    logger.info("Imputing sex ploidy...")
    ploidy_ht = impute_sex_ploidy(mt, included_intervals=capture_ht_path(data_source, freeze))
    ploidy_ht = ploidy_ht.checkpoint(ploidy_ht_path(data_source, freeze), overwrite=True)


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
    sites_ht = hl.read_table(f_stats_sites_path())
    sex_ht = default_annotate_sex(mt, ploidy_ht, sites_ht, gt_expr='LGT', aaf_expr='AF')
    sex_ht = sex_ht.checkpoint(sex_ht_path(data_source, freeze), overwrite=True)
