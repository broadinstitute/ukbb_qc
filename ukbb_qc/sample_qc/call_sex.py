import logging
from gnomad_hail.utils.sample_qc import get_ploidy_cutoffs,get_sex_expr
from gnomad_hail.utils.generic import get_reference_genome
from ukbb_qc.resources import *


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

    logger.info('Filtering to common, high-callrate, biallelic SNPs on chrX and imputing sex...')
    mt = hl.filter_intervals(mt, 
                            [hl.parse_locus_interval(x_contig, reference_genome=get_reference_genome(mt.locus).name) 
                            for x_contig in get_reference_genome(mt.locus).x_contigs]) 
    mt = mt.filter_rows((hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99) & (mt.variant_qc.AF[1] > aaf_threshold))

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


def run_impute_sex(mt: hl.MatrixTable, data_source: str, freeze: str) -> hl.Table:
    """
    Calls sex for the samples in a given mt and writes ht and output sex.txt to specified directory

    :param MatrixTable mt: MatrixTable with samples for sexcheck
    :param str data_source: One of regeneron or broad
    :return: Table with sex annotations
    :rtype: Table
    """

    if data_source == 'regeneron':
        mt = mt.filter_rows(hl.is_missing(mt.filters))
    # NOTE: Broad callset does not have filter annotations
    mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]))

    logger.info('Imputing sex ploidy...')
    ploidy_ht = impute_sex_ploidy(mt, included_intervals=capture_ht_path(data_source, freeze))
    ploidy_ht = ploidy_ht.checkpoint(ploidy_ht_path(data_source, freeze), overwrite=True)
    # NOTE: inspect ploidy_ht output before continuing

    logger.info('Imputing sex...')
    sex_ht = impute_sex(mt, ploidy_ht)
    sex_ht = sex_ht.checkpoint(sex_ht_path(data_source, freeze), overwrite=True)
