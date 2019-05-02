import logging
from gnomad_hail.utils.generic import get_reference_genome
from ukbb_qc.resources import *


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("call_sex")
logger.setLevel(logging.INFO)


def get_chrY_cov(mt: hl.MatrixTable, build: str) -> hl.MatrixTable:
    """
    Calculates mean chrY coverage, mean chr20 coverage, and normalized chrY coverage for each sample

    :param MatrixTable mt: MatrixTable containing samples to be checked
    :param str build: Reference genome build (GRCh37, GRCh38)
    :return: MatrixTable with coverage annotations
    :rtype: MatrixTable
    """

    logger.info('Filtering to chromosome 20 and non-PAR regions on chrY')
    if build == 'GRCh37': 
        chr20 = '20'
        chrY = 'Y'
    else:
        chr20 = 'chr20'
        chrY = 'chrY'

    mt = hl.filter_intervals(
                mt, 
                [hl.parse_locus_interval(chr20, reference_genome=build), 
                hl.parse_locus_interval(chrY, reference_genome=build)])
    mt = mt.filter_rows((mt.locus.contig == chr20) | (mt.locus.in_y_nonpar()),keep = True)

    logger.info('Calculating coverage on chromosome 20 and chromosome Y')
    mt = mt.annotate_cols(chrY_mean_cov=hl.agg.filter(mt.locus.contig == chrY, hl.agg.mean(mt.DP)))
    mt = mt.annotate_cols(chr20_mean_cov=hl.agg.filter(mt.locus.contig == chr20, hl.agg.mean(mt.DP)))

    logger.info('Using chr20 coverage to normalize chrY coverage')
    mt = mt.annotate_cols(normalized_Y_cov=hl.cond(mt.chr20_mean_cov > 0, mt.chrY_mean_cov/mt.chr20_mean_cov, -99))

    return mt.cols()


def annotate_sex(
                mt: hl.MatrixTable, 
                male_threshold: float = 0.75, 
                female_threshold: float = 0.5, 
                aaf_threshold = 0.05
        ) -> hl.Table:
    """
    Imputes sex and annotates matrix table with results, outputs a histogram of fstat values
    NOTE: decide on thresholds for male/female by looking at fstat plots

    :param MatrixTable mt: MatrixTable containing samples to be ascertained for sex
    :param str male_threshold: Threshold above which a sample will be called male
    :param str female_threshold: Threshold below which a sample will be called female
    :param str aaf_threshold: Minimum alternate allele frequency
    :return: Table with imputed sex annotations stashed in column annotation 'sex'
    :rtype: Table
    """

    logger.info('Filtering to chrX and imputing sex')
    mt = hl.filter_intervals(mt, 
                            [hl.parse_locus_interval(x_contig, reference_genome=get_reference_genome(mt.locus).name) 
                            for x_contig in get_reference_genome(mt.locus).x_contigs]) 
    sex_ht = hl.impute_sex(mt.GT, aaf_threshold = aaf_threshold, male_threshold = male_threshold, female_threshold = female_threshold)

    # annotate mt with imputed sex
    mt = mt.annotate_cols(**sex_ht[mt.col_key])
    sex_ht = mt.cols()

    logger.info('Updating imputed sex with normalized chrY coverage information')
    sex_ht = sex_ht.annotate(
                            ambiguous_sex = hl.is_missing(sex_ht.is_female),
                            sex_aneuploidy = (
                                    (sex_ht.is_female) & hl.is_defined(sex_ht.normalized_Y_cov) 
                                    & (sex_ht.normalized_Y_cov > 0.1)
                                    | (~sex_ht.is_female) & hl.is_defined(sex_ht.normalized_Y_cov) 
                                    & (sex_ht.normalized_Y_cov < 0.1)
                                    )
                            )
    
    sex_expr = (hl.case()
                        .when(sex_ht.ambiguous_sex, "ambiguous_sex")
                        .when(sex_ht.sex_aneuploidy, "sex_aneuploidy")
                        .when(sex_ht.is_female, "female")
                        .default("male"))

    sex_ht = sex_ht.annotate(sex = sex_expr)
    return sex_ht


def impute_sex(mt: hl.MatrixTable, build: str, data_source:str, min_callrate: float = 0.99) -> hl.Table:
    """
    Calls sex for the samples in a given mt and writes ht and output sex.txt to specified directory

    :param MatrixTable mt: MatrixTable with samples for sexcheck
    :param str build: Genome build (GRCh37, GRCh38)
    :param str data_source: One of regeneron or broad
    :param float min_callrate: Minimum callrate, default 0.99
    :return: Table with sex annotations
    :rtype: Table
    """
    logger.info('Filtering to common PASS variants (PASS: empty filters field for broad, missing filters value for regeneron)')
    if data_source == 'regeneron':
        mt = mt.filter_rows(hl.is_missing(mt.filters))
        mt = mt.filter_rows(mt.variant_qc.AF[0] > 0.05)

    else:
        mt = mt.filter_rows(mt.filters.length() == 0, keep=True)
        mt = mt.filter_rows(mt.info.AF[0] > 0.05)

    logger.info('Filtering to biallelic SNPs')
    mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]))

    logger.info('Filtering to variants with high callrate')
    mt = mt.filter_rows(hl.agg.fraction(hl.is_defined(mt.GT)) > min_callrate)

    logger.info('Starting sexcheck')
    logger.info('Getting coverage on chromosome Y')
    coverage_ht = get_chrY_cov(mt, build)

    logger.info('Annotating mt with coverage information')
    mt = mt.annotate_cols(**coverage_ht[mt.col_key]) 
   
    logger.info('Imputing sex') 
    sex_ht = annotate_sex(mt)
    sex_ht = sex_ht.select(
                        sex_ht.is_female,
                        sex_ht.f_stat,
                        sex_ht.n_called,
                        sex_ht.expected_homs,
                        sex_ht.observed_homs,
                        sex_ht.sex, 
                        sex_ht.chrY_mean_cov,
                        sex_ht.chr20_mean_cov,
                        sex_ht.normalized_Y_cov
                        )

    # NOTE: code crashes here trying to write out sex ht
    sex_ht = sex_ht.checkpoint(sex_ht_path(data_source, CURRENT_FREEZE), overwrite=True)
    return sex_ht
