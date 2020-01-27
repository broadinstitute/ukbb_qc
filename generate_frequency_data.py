from gnomad_hail.utils.annotations import * 
from ukbb_qc.resources import * 


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("generate_frequency_data")
logger.setLevel(logging.INFO)


def generate_frequency_data(mt: hl.MatrixTable, POPS_TO_REMOVE_FOR_POPMAX: List[str],
                            calculate_by_platform: bool = False,
                            ) -> hl.Table:
    '''
    Generates frequency struct annotation containing AC, AF, AN, and homozygote count for dataset stratified by population. Optional to stratify by tranche

    :param MatrixTable mt: Input MatrixTable
    :param list POPS_TO_REMOVE_FOR_POPMAX: List of populations to exclude from popmax calculations
    :param bool calculate_by_platform: Calculate frequencies per tranche
    :return: Table with frequency annotations in struct named `freq` and metadata in globals named `freq_meta`
    :rtype: Table
    '''
    if calculate_by_platform:
        additional_strata_expr = {'tranche': mt.tranche}
    else:
        additional_strata_expr = None

    logger.info('Generating frequency data...')
    mt = annotate_freq(
        mt,
        sex_expr=mt.meta.sex,
        pop_expr=mt.meta.hybrid_pop,
        additional_strata_expr=additional_strata_expr
    )

    # Select freq, FAF and popmax
    faf, faf_meta = faf_expr(mt.freq, mt.freq_meta, mt.locus, POPS_TO_REMOVE_FOR_POPMAX)
    mt = mt.select_rows(
        'freq',
        faf=faf,
        popmax=pop_max_expr(mt.freq, mt.freq_meta, POPS_TO_REMOVE_FOR_POPMAX)
    )
    mt = mt.annotate_globals(faf_meta=faf_meta)

    # Annotate quality metrics histograms, as these also require densifying
    # NOTE skipping for tranche 2, as these were already generated
    mt = mt.annotate_rows(
        **qual_hist_expr(mt.GT, mt.GQ, mt.DP, mt.AD)
    )

    return mt.rows()


def join_gnomad(ht: hl.Table, data_type: str) -> hl.Table:
    """
    Joins UKBB ht to gnomAD ht and adds gnomAD freq, popmax, and faf as annotation

    :param Table ht: Input UKBB ht
    :param str data_type: One of exomes or genomes
    :return: UKBB ht with gnomAD frequency information added as annotation
    :rtype: Table
    """
    if data_type == 'exomes':
        gnomad_ht = hl.read_table(get_gnomad_liftover_data_path(f'{data_type}', '2.1.1')).select(
            'freq', 'popmax', 'faf').select_globals(
            'freq_meta', 'popmax_index_dict', 'faf_index_dict')
        ht = ht.join(gnomad_ht, how='left')
        ht = ht.rename({'freq_1': f'gnomad_{data_type}_freq', 'popmax_1': f'gnomad_{data_type}_popmax',
                   'faf_1': f'gnomad_{data_type}_faf', 'freq_meta_1': f'gnomad_{data_type}_freq_meta',
                   'popmax_index_dict': f'gnomad_{data_type}_popmax_index_dict',
                   'faf_index_dict': f'gnomad_{data_type}_faf_index_dict'})

    else:
        gnomad_ht = hl.read_table('gs://gnomad-public/release/3.0/ht/genomes/gnomad.genomes.r3.0.sites.ht').select(
            'freq', 'popmax', 'faf').select_globals(
            'freq_meta', 'faf_index_dict')

        ht = ht.join(gnomad_ht, how='left')
        ht = ht.rename({'freq_1': f'gnomad_{data_type}_freq', 'popmax_1': f'gnomad_{data_type}_popmax',
                   'faf_1': f'gnomad_{data_type}_faf', 'freq_meta_1': f'gnomad_{data_type}_freq_meta',
                   'faf_index_dict': f'gnomad_{data_type}_faf_index_dict'})
    return ht


def main(args):
    hl.init(log='/frequency_generation.log')

    data_source = args.data_source
    freeze = args.freeze
    POPS_TO_REMOVE_FOR_POPMAX = set(args.pops.split(','))
    logger.info(f'Excluding {POPS_TO_REMOVE_FOR_POPMAX} from popmax and faf calculations')

    mt = get_ukbb_data(data_source, freeze, meta_root='meta')
    logger.info(f'mt count before filtering out low quality and non-releasable samples: {mt.count()}')
    mt = mt.filter_cols(mt.meta.releasable & mt.meta.high_quality)

    # NOTE: adapted from gnomAD v3, not sure this works yet on WES data
    # will need to finish adapting code after completing tranche 2 for sparse data
    if args.densify:
        logger.info('Reading in capture ht')
        capture_ht = hl.read_table(capture_ht_path('broad', freeze))

        # sparse mt only keyed by locus
        mt = mt.key_by('locus', 'alleles')
        mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

        logger.info('Computing adj and sex adjusted genotypes')
        mt = mt.annotate_entries(
            GT=adjusted_sex_ploidy_expr(mt.locus, mt.GT, mt.meta.sex_karyotype),
            adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD)
        )

        mt = hl.experimental.densify(mt)
        mt = mt.filter_rows(hl.len(mt.alleles) > 1)
        mt = mt.filter_rows(hl.is_defined(capture_ht[mt.locus]))

    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    logger.info(f'mt count after filtering out low quality samples and their variants: {mt.count()}')

    if args.filter_related:
        mt = mt.filter_cols(~mt.meta.related_filter)
        mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
        logger.info(f'mt count after filtering out related samples and their variants: {mt.count()}')

    if args.by_platform:
        logger.info('Annotating sample with tranche information (from array sample map ht)')
        sample_map_ht = hl.read_table(array_sample_map_ht(data_source, freeze)).select('batch.c')
        mt = mt.annotate_cols(**sample_map_ht[mt.s])
        mt = mt.transmute_cols(tranche=mt['batch.c'])

    if args.calculate_frequencies:
        logger.info('Calculating frequencies')
        ht = generate_frequency_data(mt, POPS_TO_REMOVE_FOR_POPMAX, args.by_platform)
        write_temp_gcs(ht, var_annotations_ht_path(data_source, freeze, f'ukb_freq_hybrid{"_unrelated" if args.filter_related else ""}'), args.overwrite)

    if args.join_gnomad:
        ht = hl.read_table(var_annotations_ht_path(data_source, freeze, f'ukb_freq{"_unrelated" if args.filter_related else ""}'))

        logger.info('Joining UKBB ht to gnomAD exomes and genomes liftover hts')
        ht = join_gnomad(ht, 'exomes')
        ht = join_gnomad(ht, 'genomes')
        write_temp_gcs(ht, var_annotations_ht_path(data_source, freeze, f'join_freq_hybrid{"_unrelated" if args.filter_related else ""}'), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)
    parser.add_argument('--pops', help='Pops to exclude for popmax', default='asj,fin,oth')
    parser.add_argument('-d', '--densify', help='Densify data (sparse data only)', action='store_true')
    parser.add_argument('--calculate_frequencies', help='Calculate most frequency data', action='store_true')
    parser.add_argument('--by_platform', help='Also calculate frequencies by platform', action='store_true')
    parser.add_argument('--filter_related', help='Filter related individuals before calculating allele frequencies', action='store_true')
    parser.add_argument('-j', '--join_gnomad', help='Join table with gnomAD tables to get gnomAD frequencies', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
