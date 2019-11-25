import logging
from gnomad_hail.utils.generic import get_reference_genome
from ukbb_qc.resources import *


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("call_sex")
logger.setLevel(logging.INFO)

def impute_sex_ploidy(
        mt: hl.MatrixTable,
        included_intervals: Optional[hl.Table],
        normalization_contig: str = 'chr20',
        chr_x: Optional[str] = None,
        chr_y: Optional[str] = None,
) -> hl.Table: # TODO: For exomes, calling intervals need to be added
    """
    Imputes sex ploidy from a sparse Matrix Table by normalizing the coverage of chromosomes X and Y using
    the coverage of an autosomal chromosome (by default chr20).
    Coverage is computed using the median block coverage (summed over the block size) and the non-ref coverage at non-ref genotypes.

    :param MatrixTable mt: Input sparse Matrix Table
    :param Table excluded_intervals: An optional table of intervals to exclude from the computation.
    :param str normalization_contig: Which chromosome to normalize by
    :param str chr_x: Optional X Chromosome contig name (by default uses the X contig in the reference)
    :param str chr_y: Optional Y Chromosome contig name (by default uses the Y contig in the reference)
    :return: Table with mean coverage over chromosomes 20, X and Y and sex chromosomes ploidy based on normalized coverage.
    :rtype: Table
    """

    ref = get_reference_genome(mt.locus, add_sequence=True)
    if chr_x is None:
        if len(ref.x_contigs) != 1:
            raise NotImplementedError(
                "Found {0} X chromosome contigs ({1}) in Genome reference. sparse_impute_sex_ploidy currently only supports a single X chromosome contig. Please use the `chr_x` argument to  specify which X chromosome contig to use ".format(
                    len(ref.x_contigs),
                    ",".join(ref.x_contigs)
                )
            )
        chr_x = ref.x_contigs[0]
    if chr_y is None:
        if len(ref.y_contigs) != 1:
            raise NotImplementedError(
                "Found {0} Y chromosome contigs ({1}) in Genome reference. sparse_impute_sex_ploidy currently only supports a single Y chromosome contig. Please use the `chr_y` argument to  specify which Y chromosome contig to use ".format(
                    len(ref.y_contigs),
                    ",".join(ref.y_contigs)
                )
            )
        chr_y = ref.y_contigs[0]

    def get_contig_size(contig: str) -> int:
        contig = hl.utils.range_table(ref.contig_length(contig), n_partitions=int(ref.contig_length(contig) / 500_000))
        contig.show()
        contig = contig.annotate(
            locus=hl.locus(contig=contig, pos=contig.idx + 1, reference_genome=ref)
        )
        contig = contig.filter(contig.locus.sequence_context().lower() != 'n')

        print(contig.count())
        if contig in ref.y_contigs:
            contig = contig.filter(contig.locus.in_y_nonpar())
        print(contig.count())

        contig = contig.key_by('locus')
        if included_intervals is not None:
            contig = contig.filter(hl.is_defined(included_intervals[contig.key]))
        contig_size = contig.count()
        logger.info(f"Contig {contig} has {contig_size} bases for coverage.")
        return contig_size

    def get_chr_dp_ann(chr: str) -> hl.Table:
        contig_size = get_contig_size(chr)
        chr_mt = hl.filter_intervals(mt, [hl.parse_locus_interval(chr)])
        print(chr_mt.count())
        chr_mt = chr_mt.filter(hl.is_defined(included_intervals[chr_mt.key]))
        print(chr_mt.count())
        return chr_mt.select_cols(**{
            f'{chr}_mean_dp': hl.agg.sum(hl.cond(chr_mt.LGT.is_hom_ref(), chr_mt.DP * (chr_mt.END - chr_mt.locus.position), chr_mt.DP)) / contig_size
        }).cols()

    normalization_chrom_dp = get_chr_dp_ann(normalization_contig)
    chrX_dp = get_chr_dp_ann(chr_x)
    chrY_dp = get_chr_dp_ann(chr_y)

    ht = normalization_chrom_dp.annotate(
        **chrX_dp[normalization_chrom_dp.key],
        **chrY_dp[normalization_chrom_dp.key],
    )

    return ht.annotate(
        **{
            f'{chr_x}_ploidy': ht[f'{chr_x}_mean_dp'] / (ht[f'{normalization_contig}_mean_dp'] / 2),
            f'{chr_y}_ploidy': ht[f'{chr_y}_mean_dp'] / (ht[f'{normalization_contig}_mean_dp'] / 2)
        }
    )

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
    mt = mt.filter_rows((hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99) & (mt.variant_qc.AF[1] > 0.05))

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
    # coverage_ht = get_chrY_cov(mt)
    # mt = mt.annotate_cols(**coverage_ht[mt.col_key])

    ploidy_ht = impute_sex_ploidy(mt, capture_ht_path())
    ploidy_ht.show()
    ploidy_ht = ploidy_ht.checkpoint(sex_ht_path(data_source, freeze), overwrite=True)
    # TODO: inspect ploidy_ht output before fixing line 142 to correctly mark sex aneuploidies
    # mt = mt.annotate_cols(**ploidy_ht[mt.col_key])

    # sex_ht = impute_sex(mt)
    # sex_ht = sex_ht.checkpoint(sex_ht_path(data_source, freeze), overwrite=True)
    # sex_ht.export(sex_tsv_path(data_source, freeze))
