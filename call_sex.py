import logging
from gnomad_hail.utils.generic import get_reference_genome
from ukbb_qc.resources import *


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("call_sex")
logger.setLevel(logging.INFO)


def get_chrY_cov(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Calculates mean chrY coverage, mean chr20 coverage, and normalized chrY coverage for each sample

    :param MatrixTable mt: MatrixTable containing samples to be checked
    :return: MatrixTable with coverage annotations
    :rtype: MatrixTable
    """

    mt = mt.filter_rows((mt.locus.contig == 'chr20') | (mt.locus.in_y_nonpar()), keep=True)
    mt = mt.annotate_cols(chrY_frac_cov=hl.agg.filter(mt.locus.contig == 'chrY', hl.agg.fraction(hl.is_defined(mt.DP))),
                          chr20_frac_cov=hl.agg.filter(mt.locus.contig == 'chr20', hl.agg.fraction(hl.is_defined(mt.DP)))
                          )
    # NOTE: switching from mean_cov to frac_cov for Broad callset -- will need to test on the next Regeneron VCF
    mt = mt.annotate_cols(normalized_Y_frac=hl.cond(mt.chr20_frac_cov > 0,
                                                    hl.cond(mt.chrY_frac_cov/mt.chr20_frac_cov > 0,
                                                            mt.chrY_frac_cov/mt.chr20_frac_cov, 0), -99))
    return mt.cols()


def impute_sex(mt: hl.MatrixTable, male_threshold: float = 0.75, female_threshold: float = 0.5, aaf_threshold = 0.05) -> hl.Table:
    """
    Imputes sex and annotates matrix table with results, outputs a histogram of fstat values
    NOTE: decide on thresholds for male/female by looking at fstat plotse

    :param MatrixTable mt: MatrixTable containing samples to be ascertained for sex
    :param str male_threshold: Threshold above which a sample will be called male
    :param str female_threshold: Threshold below which a sample will be called female
    :param str aaf_threshold: Minimum alternate allele frequency
    :return: Table with imputed sex annotations stashed in column annotation 'sex'
    :rtype: Table
    """

    logger.info('Filtering to common, high-callrate, biallelic SNPs on chrX and imputing sex...')
    mt = hl.filter_intervals(mt, 
                            [hl.parse_locus_interval(x_contig, reference_genome=get_reference_genome(mt.locus).name) 
                            for x_contig in get_reference_genome(mt.locus).x_contigs]) 
    mt = mt.filter_rows((hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99) & (mt.variant_qc.AF[0] > 0.05))

    sex_ht = hl.impute_sex(mt.GT, aaf_threshold=aaf_threshold, male_threshold=male_threshold, female_threshold=female_threshold)
    mt = mt.annotate_cols(**sex_ht[mt.col_key])
    sex_ht = mt.cols()  # NOTE: sex_ht now has coverage annotations from MT

    logger.info('Updating imputed sex with normalized chrY coverage information')
    sex_ht = sex_ht.annotate(ambiguous_sex=hl.is_missing(sex_ht.is_female),
                             sex_aneuploidy=(
                                    (sex_ht.is_female) & hl.is_defined(sex_ht.normalized_Y_frac)
                                    & (sex_ht.normalized_Y_frac > 0.25)
                                    | (~sex_ht.is_female) & hl.is_defined(sex_ht.normalized_Y_frac)
                                    & (sex_ht.normalized_Y_frac <= 0.25)
                                    ))
    
    sex_expr = (hl.case()
                        .when(sex_ht.ambiguous_sex, "ambiguous_sex")
                        .when(sex_ht.sex_aneuploidy, "sex_aneuploidy")
                        .when(sex_ht.is_female, "female")
                        .default("male"))

    sex_ht = sex_ht.annotate(sex=sex_expr)
    return sex_ht


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

    logger.info('Computing coverage on chromosome Y...')
    coverage_ht = get_chrY_cov(mt)
    mt = mt.annotate_cols(**coverage_ht[mt.col_key])

    sex_ht = impute_sex(mt)
    sex_ht = sex_ht.checkpoint(sex_ht_path(data_source, freeze), overwrite=True)
    sex_ht.export(sex_tsv_path(data_source, freeze))
