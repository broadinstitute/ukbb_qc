import logging
import hail as hl
from ukbb_qc.resources import *

logging.basicConfig(format='%(asctime)s (%(name)s %(lineno)s): %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('variant_check')
logger.setLevel(logging.INFO)


def remove_unique_samples(b_mt: hl.MatrixTable, r_mt: hl.MatrixTable) -> hl.MatrixTable:
    '''
    Removes regeneron unique samples and variants unique to those samples from input MatrixTable

    :param MatrixTable b_mt: broad mt
    :param MatrixTable r_mt: regeneron mt
    :return: regeneron mt with unique samples/variants removed
    :rtype: MatrixTable
    '''

    # update other function to do this, but general
    logger.info('Removing regeneron unique samples and variants unique to those samples')
    b_s = b_mt.aggregate_cols(hl.agg.collect_as_set(b_mt.s))
    r_s = r_mt.aggregate_cols(hl.agg.collect_as_set(r_mt.s))
    r_only_s = hl.literal(r_s.difference(b_s))
    logger.info(f'# of broad samples: {len(b_s)}')
    logger.info(f'# of regeneron samples: {len(r_s)}')
    logger.info(f'# regeneron unique samples: {hl.eval(hl.len(r_only_s))}')
    
    r_mt = r_mt.filter_cols(~r_only_s.contains(r_mt.s))
    r_mt = r_mt.annotate_rows(non_refs=hl.agg.count_where(r_mt.GT.is_non_ref()))
    r_mt = r_mt.filter_rows(r_mt.non_refs > 0)
    logger.info(f'regeneron mt count after removing unique samples: {r_mt.count()}')
    return r_mt.drop('non_refs')


def process_hardcalls(temp_prefix: str, freeze: int, ow: str) -> None:
    '''
    Reads in broad and regeneron MatrixTables, removes Regneron unique samples, filters to adj/removes AC0 variants, and converts to Tables

    :param str temp_prefix: Path to temp bucket (output directory)
    :param str freeze: Freeze number
    :param str ow: Whether to overwrite
    :return: None
    :rtype: None
    '''
    logger.info('Reading in hardcall mts')
    r_mt = hl.read_matrix_table(hardcalls_mt_path('regeneron', freeze, False))
    b_mt = hl.read_matrix_table(hardcalls_mt_path('broad', freeze, False))
    logger.info(f'regeneron mt count {r_mt.count()}')
    logger.info(f'broad mt count {b_mt.count()}')

    logger.info('Removing Regeneron unique samples')
    r_mt = remove_unique_samples(b_mt, r_mt)

    logger.info('Filtering to adj')
    r_mt = filter_to_adj(r_mt)
    b_mt = filter_to_adj(b_mt)

    logger.info('Re-adding variant qc')
    r_mt = r_mt.drop('variant_qc')
    r_mt = hl.variant_qc(r_mt)
    b_mt = hl.variant_qc(b_mt)
    
    logger.info('Splitting multi hts')
    r_mt = hl.split_multi_hts(r_mt)
    b_mt = hl.split_multi_hts(b_mt)

    logger.info('Repartitioning mts')
    r_mt = r_mt.naive_coalesce(10000)
    b_mt = b_mt.naive_coalesce(10000)
   
    logger.info('Filtering out variants with AC == 0')
    r_mt = r_mt.filter_rows(r_mt.variant_qc.AC[r_mt.a_index] > 0)
    b_mt = b_mt.filter_rows(b_mt.variant_qc.AC[b_mt.a_index] > 0) 

    logger.info('Writing out mts')
    r_mt.write(f'{temp_prefix}/regeneron_hardcalls_adj.mt', overwrite=ow)
    b_mt.write(f'{temp_prefix}/broad_hardcalls_adj.mt', overwrite=ow)

    logger.info('Converting mts to hts')
    #r_ht = r_mt.rows().drop('rsid', 'qual', 'filters', 'info')
    #r_ht = r_ht.transmute(r_index=r_ht.a_index, r_was_split=r_ht.was_split)
    #b_ht = b_mt.rows().drop('rsid', 'qual', 'filters', 'variant_qc')
    #b_ht = b_ht.transmute(b_index=b_ht.a_index, b_was_split=b_ht.was_split)
    r_ht.write(f'{temp_prefix}/regeneron_hardcalls_adj.ht', overwrite=ow)
    b_ht.write(f'{temp_prefix}/broad_hardcalls_adj.ht', overwrite=ow)


def import_capture_regions(temp_prefix: str) -> None:
    '''
    Imports capture region text file into ht

    :param str temp_prefix: Path to temp bucket (output directory)
    :return: None
    :rtype: None
    '''

    logger.info('Importing capture table')
    capture_ht = hl.import_table(
                    ukbb_calling_intervals_summary, impute=True,
                    no_header=False)

    # seqnames	start	end	width	strand	target_type	region_type	target_id
    #chr1	11719	12377	659	*	processed_transcript|transcribed_unprocessed_pseudogene|ice_target_1	gnomad	target_1
    capture_ht = capture_ht.transmute(locus=hl.parse_locus_interval(
                                    hl.format('[%s:%s-%s]', 
                                               capture_ht.seqnames, 
                                               capture_ht.start, 
                                               capture_ht.end),
                                    reference_genome='GRCh38')).drop(
                                    'width', 'strand').key_by('locus')

    logger.info('Writing capture ht')
    capture_ht.write(capture_ht_path('regeneron'), overwrite=True)


def get_variant_breakdown(ht: hl.Table, temp_prefix: str) -> None:
    '''
    Gets variant breakdown of unique variants (# multiallelic/biallelic, # snps, capture regions)

    :param Table ht: Input ht
    :param str data_source: Dataset (broad or regeneron)
    :param str temp_prefix: Path to temp bucket (where capture ht is stored)
    :return: None
    :rtype: None
    '''
    logger.info('Getting capture region breakdown')
    capture_ht = hl.read_table(f'{temp_prefix}/capture.ht')
    ht = ht.annotate(**capture_ht[ht.locus])
    logger.info('Getting variant count breakdowns')
    counts = ht.aggregate(hl.struct(
                            total=hl.agg.count(),
                            multi=hl.agg.count_where(ht.was_split),
                            bi=hl.agg.count_where(~ht.was_split),
                            m_snp=hl.agg.count_where(ht.was_split & hl.is_snp(ht.alleles[0], ht.alleles[1])),
                            b_snp=hl.agg.count_where(~ht.was_split & hl.is_snp(ht.alleles[0], ht.alleles[1])),
                            broad_total=hl.agg.count_where(ht.region_type == 'gnomad'),
                            ukbb_annex=hl.agg.count_where(ht.region_type == 'ukbb_annex'),
                            ukbb_islands=hl.agg.count_where(ht.region_type == 'ukbb_islands')
                            ))
    logger.info(f'Variant count breakdown: {counts}')
    

def subset_calling_regions(t: Union[hl.MatrixTable, hl.Table]) -> Union[hl.MatrixTable, hl.Table]:
    '''
    Subsets input data to calling regions provided by Laura

    :param MatrixTable/Table t: Input data to be subset
    :return: Subset data
    :rtype: MatrixTable/Table
    '''
    calling_ht = hl.read_table('gs://broad-ukbb/resources/xgen_plus_spikein.Homo_sapiens_assembly38.targets.pad50.intervals.ht')
    if isinstance(t, hl.Table):
        t = t.annotate(**calling_ht[t.locus])
        return t.filter(hl.is_defined(t.strand)) 
    else:
        t = t.annotate_rows(**calling_ht[t.locus])
        return t.filter_rows(hl.is_defined(t.strand))


def get_variants(temp_prefix: str, unique: bool, ow: bool) -> None:
    '''
    Gets variants either unique to or shared between broad and regeneron datasets and writes resulting Table(s)
    If getting unique variants, also get variant type breakdown

    :param str temp_prefix: Prefix to temp bucket (output directory)
    :param bool unique: Whether to get unique or shared variants
    :param bool ow: Whether to overwrite output
    :return: None
    :rtype: None
    '''

    logger.info('Reading in hts')
    b_ht = hl.read_table(f'{temp_prefix}/broad_hardcalls_adj.ht')
    r_ht = hl.read_table(f'{temp_prefix}/regeneron_hardcalls_adj.ht')

    if unique:
        logger.info('Getting variants unique to both hts')
        r_uniq_ht = r_ht.anti_join(b_ht)
        b_uniq_ht = b_ht.anti_join(r_ht)

        logger.info('Writing out hts')
        r_uniq_ht = r_uniq_ht.checkpoint(f'{temp_prefix}/regeneron_unique_variants.ht', overwrite=ow)
        logger.info(f'# lines in regeneron unique ht: {r_uniq_ht.count()}')
        b_uniq_ht = b_uniq_ht.checkpoint(f'{temp_prefix}/broad_unique_variants.ht', overwrite=ow)
        logger.info(f'# lines in broad unique ht: {b_uniq_ht.count()}')
    
        logger.info('Subsetting to calling regions')
        r_uniq_ht = subset_calling_regions(r_uniq_ht)
        logger.info(f'# lines in regeneron unique ht in calling regions: {r_uniq_ht.count()}')
        b_uniq_ht = subset_calling_regions(b_uniq_ht)
        logger.info(f'# lines in broad unique ht in calling regions: {b_uniq_ht.count()}')

        logger.info('Getting variant type breakdown')
        get_variant_breakdown(b_uniq_ht, 'broad', temp_prefix)
        get_variant_breakdown(r_uniq_ht, 'regeneron', temp_prefix)

    else:
        logger.info('Getting variants shared between both datasets')
        shared_ht = b_ht.semi_join(r_ht)
        shared_ht = shared_ht.checkpoint(f'{temp_prefix}/broad_reg_shared_variants.ht', overwrite=ow)
        logger.info(f'# shared variants: {shared_ht.count()}')
       
 
def filter_to_type(mt: hl.MatrixTable, variant_type: str, multi) -> hl.MatrixTable:
    '''
    Filters mt to specified type

    :param MatrixTable mt: Input mt
    :param str variant_type: Variant type
    :param bool multi: Whether to filter to multiallelics
    :return: MatrixTable filtered to variant type
    :rtype: MatrixTable
    '''

    logger.info(f'mt count pre-filtration: {mt.count()}')
    if multi:
        mt = mt.filter_rows(mt.was_split)
    else:
        mt = mt.filter_rows(~mt.was_split)

    if variant_type == 'SNP':
        logger.info('Filtering hardcalls mts to snps')
        mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1])) 
    else:
        logger.info('Filtering to indels')
        mt =  mt.filter_rows(hl.is_indel(mt.alleles[0], mt.alleles[1]))

    logger.info(f'Finished filtering to {variant_type}')
    logger.info(f'mt count post-filtration: {mt.count()}')
    return mt


def check_concordance(temp_prefix: str, variant_type: str, multi: bool, ow: str) -> None:
    '''
    Checks concordance between broad and regeneron

    :param str temp_prefix: Path to temp bucket
    :param str variant_type: Variant type (SNP/indel)
    :param bool multi: Whether to filter to multiallelics (T) or biallelics (F)
    :param str ow: Whether to overwrite output files
    :return: None
    :rtype: None
    '''

    logger.info(f'Working on {variant_type} multi={multi} concordance')
    r_mt = hl.read_matrix_table(f'{temp_prefix}/regeneron_hardcalls_adj.mt')
    b_mt = hl.read_matrix_table(f'{temp_prefix}/broad_hardcalls_adj.mt')

    b_mt = filter_to_type(b_mt, variant_type, multi)
    r_mt = filter_to_type(r_mt, variant_type, multi)
   
    logger.info('Repartitioning mts') 
    b_mt = b_mt.repartition(2500, shuffle=False)
    r_mt = r_mt.repartition(2500, shuffle=False)

    b_mt = b_mt.checkpoint(f'{temp_prefix}/broad_concordance.mt', overwrite=True)
    r_mt = r_mt.checkpoint(f'{temp_prefix}/regeneron_concordance.mt', overwrite=True)

    logger.info('Checking concordance between broad and regeneron data...')
    summary, samples, variants = hl.concordance(b_mt, r_mt)
    #samples.write(f'{temp_prefix}/{variant_type}_b_r_sample_{variant_type}_{multi}_concordance_temp.ht', overwrite=ow)
    variants.write(f'{temp_prefix}/{variant_type}_b_r_variant_{variant_type}_{multi}_concordance_temp.ht', overwrite=ow)

    # concordance: 
    # index 0 is no data (missing variant), index 1 is no call (missing GT call), 
    # index 2 is homref, index 3 is het, index 4 is homvar

 
def main(args):

    #hl.init(default_reference='GRCh38', log='/variant_check.log', tmp_dir='hdfs:///variant_check.tmp/')
    hl.init(default_reference='GRCh38', log='/variant_check.log')    
    #hl._set_flags(newaggs=None)

    temp_prefix = 'gs://broad-ukbb/broad.freeze_4/temp'

    if args.process_hardcalls:
        process_hardcalls(temp_prefix, args.freeze, args.overwrite)

    if args.import_regions:
        import_capture_regions(temp_prefix)

    if args.get_unique_var:
        get_variants(temp_prefix, True, args.overwrite)

    if args.get_shared_var:
        get_variants(temp_prefix, False, args.overwrite)

    if args.check_concordance:
        check_concordance(temp_prefix, 'SNP', True, args.overwrite)
        check_concordance(temp_prefix, 'SNP', False, args.overwrite)
        check_concordance(temp_prefix, 'indel', True, args.overwrite)
        check_concordance(temp_prefix, 'indel', False, args.overwrite)
 

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--freeze', help='Current freeze #', default=CURRENT_FREEZE, type=int)
    parser.add_argument('--process_hardcalls', help='Read in and process hardcall mts', action='store_true')
    parser.add_argument('--import_regions', help='Import capture regions and write out capture ht', action='store_true')
    parser.add_argument('--get_unique_var', help='Get variants unique to each dataset', action='store_true')
    parser.add_argument('--get_shared_var', help='Get variants shared between broad/regeneron datasets', action='store_true')
    parser.add_argument('--check_concordance', help='Check concordance between broad/regeneron datasets', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
