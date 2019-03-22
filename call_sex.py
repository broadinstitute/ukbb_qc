# stolen from Kristen (kml_call_sex notebook) with minor edits for UKBB
import hail as hl
import hail.expr.aggregators as agg
import logging


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("call_sex")
logger.setLevel(logging.INFO)


def impute_sex(mt: hl.MatrixTable, build: str, outdir: str) -> hl.Table:
    """
    Calls sex for the samples in a given mt and writes ht and output sex.txt to specified directory
    :param MatrixTable mt: MatrixTable with samples for sexcheck
    :param str build: Genome build (GRCh37, GRCh38)
    :param str outdir: Path to output directory
    :return: Table with sex annotations
    :rtype: Table
    """

    def get_y_cov(mt: hl.MatrixTable) -> hl.MatrixTable:
        """
        Calculates mean chrY coverage, mean chr20 coverage, and normalized chrY coverage for each sample
        :param MatrixTable mt: MatrixTable containing samples to be checked
        :return: MatrixTable with coverage annotations
        :rtype: MatrixTable
        """
        
        logger.info('Filtering to chromosome 20 and non-PAR regions on chrY')
        mt = mt.filter_rows((mt.locus.contig == "chr20") | (mt.locus.in_y_nonpar()),keep = True)
        
        logger.info('Filtering to PASS variants')
        mt = mt.filter_rows(hl.is_missing(mt.filters))
        
        logger.info('Filtering to common SNVs')
        mt = mt.filter_rows(mt.variant_qc.AF[0] > 0.05)
       
        logger.info('Calculating coverage on chromosome 20 and chromosome Y')
        mt = mt.annotate_cols(
        chry_cov = hl.agg.filter(mt.locus.contig == "chrY",hl.agg.mean(mt.DP)))
        mt = mt.annotate_cols(
        chr20_cov = hl.agg.filter(mt.locus.contig == "chr20",hl.agg.mean(mt.DP)))

        logger.info('Using chr20 coverage to normalize chrY coverage')
        mt = mt.annotate_cols(normalized_y_coverage = mt.chry_cov/mt.chr20_cov)
        return(mt)


    def annotate_sex(
                mt: hl.MatrixTable, 
                build: str, 
                outdir: str, 
                male_threshold: float = 0.75, 
                female_threshold: float = 0.5, 
                aaf_threshold = 0.05
        ) -> hl.Table:
        """
        Imputes sex and annotates matrix table with results, outputs a histogram of fstat values
        NOTE: decide on thresholds for male/female by looking at fstat plots
        :param MatrixTable mt: MatrixTable containing samples to be ascertained for sex
        :param str build: Genome build (GRCh37, GRCh38)
        :param str outdir: Directory to output results
        :param str male_threshold: Threshold above which a sample will be called male
        :param str female_threshold: Threshold below which a sample will be called female
        :param str aaf_threshold: Minimum alternate allele frequency
        :return: Table with imputed sex annotations stashed in column annotation 'sex'
        :rtype: Table
        """

        # filter to the X chromosome and impute sex
        if build == "GRCh37":
            mt = hl.filter_intervals(mt, [hl.parse_locus_interval('X', reference_genome="GRCh37")])
        else:
            mt = hl.filter_intervals(mt, [hl.parse_locus_interval('chrX', reference_genome="GRCh38")])
        
        sex_ht = hl.impute_sex(mt.GT, aaf_threshold = aaf_threshold, male_threshold = male_threshold, female_threshold = female_threshold)

        # annotate mt with imputed sex
        mt = mt.annotate_cols(**sex_ht[mt.col_key])
        sex_ht = mt.cols()

        # update imputed sex with normalized chrY coverage info
        sex_ht = sex_ht.annotate(
                                ambiguous_sex = hl.is_missing(sex_ht.is_female),
                                sex_aneuploidy = (
                                        (sex_ht.is_female) & hl.is_defined(sex_ht.normalized_y_coverage) 
                                        & (sex_ht.normalized_y_coverage > 0.1)
                                        | (~sex_ht.is_female) & hl.is_defined(sex_ht.normalized_y_coverage) 
                                        & (sex_ht.normalized_y_coverage < 0.1)
                                        )
                                )
        
        sex_expr = (hl.case()
                            .when(sex_ht.ambiguous_sex, "ambiguous_sex")
                            .when(sex_ht.sex_aneuploidy, "sex_aneuploidy")
                            .when(sex_ht.is_female, "female")
                            .default("male"))

        sex_ht = sex_ht.annotate(sex = sex_expr)
        return sex_ht

    logger.info('Starting sexcheck')
    logger.info('Getting coverage on chromosome Y')
    mt = get_y_cov(mt)
   
    logger.info('Imputing sex') 
    sex_ht = annotate_sex(mt, build, outdir)
    sex_ht = sex_ht.select(
                        sex_ht.is_female,
                        sex_ht.f_stat,
                        sex_ht.n_called,
                        sex_ht.expected_homs,
                        sex_ht.observed_homs,
                        sex_ht.sex, 
                        sex_ht.chry_cov,
                        sex_ht.chr20_cov,
                        sex_ht.normalized_y_coverage
                        )

    logger.info('Writing out sex ht')
    sex_ht_path = outdir + '/sex.ht'
    sex_ht.write(sex_ht_path, overwrite=True)
    sex_ht = hl.read_table(sex_ht_path)

    logger.info('Writing out sex.txt')
    outfile = outdir + '/sex.txt'
    sex_ht.export(outfile)
    return sex_ht
