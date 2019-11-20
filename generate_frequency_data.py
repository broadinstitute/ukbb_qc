from gnomad_hail import *
from gnomad_qc.annotations.generate_frequency_data import add_faf_expr,generate_downsamplings_cumulative
from ukbb_qc.resources import * 
from collections import Counter


DOWNSAMPLINGS = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000,
                 55000, 60000, 65000, 70000, 75000, 80000, 85000, 90000, 95000, 100000, 110000, 120000]
POPS_TO_REMOVE_FOR_POPMAX = ['asj', 'fin', 'oth']
F_CUTOFF = 0.05


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("generate_frequency_data")
logger.setLevel(logging.INFO)


def generate_frequency_data(mt: hl.MatrixTable, calculate_downsampling: bool = False,
                            calculate_by_platform: bool = False) -> Tuple[hl.Table, hl.Table]:
    """
    :param MatrixTable mt: Input MatrixTable
    :param bool calculate_downsampling: Calculate frequencies for downsampled data
    :param bool calculate_by_platform: Calculate frequencies for PCR-free data
    """
    if calculate_downsampling:
        mt, downsamplings = generate_downsamplings_cumulative(mt)
        logger.info(f'Got {len(downsamplings)} downsamplings: {downsamplings}')
    cut_dict = {'pop': hl.agg.filter(hl.is_defined(mt.meta.gnomad_pc_project_pop), hl.agg.counter(mt.meta.gnomad_pc_project_pop)),
                'sex': hl.agg.filter(hl.is_defined(mt.meta.sex), hl.agg.collect_as_set(mt.meta.sex)),
                'subpop': hl.agg.filter(hl.is_defined(mt.meta.hybrid_pop) & hl.is_defined(mt.meta.gnomad_pc_project_pop),
                                        hl.agg.collect_as_set(hl.struct(subpop=mt.meta.hybrid_pop, pop=mt.meta.gnomad_pc_project_pop)))
                }

    if calculate_by_platform:
        cut_dict['platform'] = hl.agg.filter(hl.is_defined(mt.meta.qc_platform),
                                             hl.agg.collect_as_set(mt.meta.qc_platform))
    cut_data = mt.aggregate_cols(hl.struct(**cut_dict))

    sample_group_filters = [({}, True)]
    sample_group_filters.extend([
        ({'pop': pop}, mt.meta.gnomad_pc_project_pop == pop) for pop in cut_data.pop
    ] + [
        ({'sex': sex}, mt.meta.sex == sex) for sex in cut_data.sex
    ] + [
        ({'pop': pop, 'sex': sex}, (mt.meta.sex == sex) & (mt.meta.gnomad_pc_project_pop == pop))
        for sex in cut_data.sex for pop in cut_data.pop
    ] + [
        ({'subpop': subpop.subpop, 'pop': subpop.pop},
         mt.meta.hybrid_pop == subpop.subpop)
        for subpop in cut_data.subpop
    ])

    if calculate_by_platform:
        sample_group_filters.extend([
            ({'platform': str(platform)}, mt.meta.qc_platform == platform)
            for platform in cut_data.platform
        ])

    if calculate_downsampling:
        sample_group_filters.extend([
            ({'downsampling': str(ds), 'pop': 'global'},
             mt.downsampling.global_idx < ds) for ds in downsamplings
        ])
        sample_group_filters.extend([
            ({'downsampling': str(ds), 'pop': pop},
             (mt.downsampling[f'{pop_field}_idx'] < ds) & (mt.meta.gnomad_pc_project_pop == pop))
            for ds in downsamplings for pop, pop_count in cut_data.gnomad_pc_project_pop.items() if ds <= pop_count
        ])
    mt = mt.select_cols(group_membership=[x[1] for x in sample_group_filters], age=mt.meta.age)
    mt = mt.select_rows()

    def get_meta_expressions(sample_group_filters):
        meta_expressions = []
        for i in range(len(sample_group_filters)):
            subgroup_dict = sample_group_filters[i][0]
            subgroup_dict['group'] = 'adj'
            meta_expressions.append(subgroup_dict)
        meta_expressions.insert(1, {'group': 'raw'})
        return meta_expressions

    def get_freq_expressions(mt, n_groups):

        adj_freq_expressions = hl.agg.array_agg(
            lambda i: hl.agg.filter(mt.group_membership[i] & mt.adj, hl.agg.call_stats(mt.GT, mt.alleles)),
            hl.range(n_groups)
        )

        # Insert raw as the second element of the array
        return adj_freq_expressions[:1].extend([
            hl.agg.call_stats(mt.GT, mt.alleles)
        ]).extend(
            adj_freq_expressions[1:]
        ).map(
            lambda cs: cs.annotate(
                AC=cs.AC[1],
                AF=cs.AF[1],
                homozygote_count=cs.homozygote_count[1]
            )
        )

    frequency_expression = get_freq_expressions(mt, len(sample_group_filters))
    print(f'Calculating {len(sample_group_filters) + 1} aggregators...')
    global_expression = {
        'freq_meta': get_meta_expressions(sample_group_filters)
    }
    mt = mt.annotate_rows(freq=frequency_expression,
                          age_hist_het=hl.agg.filter(mt.adj & mt.GT.is_het(), hl.agg.hist(mt.age, 30, 80, 10)),
                          age_hist_hom=hl.agg.filter(mt.adj & mt.GT.is_hom_var(), hl.agg.hist(mt.age, 30, 80, 10)))
    if calculate_downsampling: global_expression['downsamplings'] = downsamplings
    mt = mt.annotate_globals(**global_expression)
    sample_data = mt.cols()

    pops = set(cut_data.pop.keys())
    [pops.discard(x) for x in POPS_TO_REMOVE_FOR_POPMAX]

    mt = mt.annotate_rows(popmax=add_popmax_expr(mt.freq, mt.freq_meta, populations=pops),
                          faf=add_faf_expr(mt.freq, mt.freq_meta, mt.locus, populations=pops))

    return mt.rows(), sample_data


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
            'freq_meta', 'freq_index_dict', 'popmax_index_dict', 'faf_index_dict')
        ht = ht.join(gnomad_ht, how='left')
        ht = ht.rename({'freq_1': f'gnomad_{data_type}.freq', 'popmax_1': f'gnomad_{data_type}.popmax',
                   'faf_1': f'gnomad_{data_type}.faf', 'freq_meta_1': f'gnomad_{data_type}.freq_meta',
                   'freq_index_dict': f'gnomad_{data_type}.freq_index_dict',
                   'popmax_index_dict': f'gnomad_{data_type}.popmax_index_dict',
                   'faf_index_dict': f'gnomad_{data_type}.faf_index_dict'})

    else:
        gnomad_ht = hl.read_table('gs://gnomad-public/release/3.0/ht/genomes/gnomad.genomes.r3.0.sites.ht').select(
            'freq', 'popmax', 'faf').select_globals(
            'freq_meta', 'freq_index_dict', 'faf_index_dict')

        ht = ht.join(gnomad_ht, how='left')
        ht = ht.rename({'freq_1': f'gnomad_{data_type}.freq', 'popmax_1': f'gnomad_{data_type}.popmax',
                   'faf_1': f'gnomad_{data_type}.faf', 'freq_meta_1': f'gnomad_{data_type}.freq_meta',
                   'freq_index_dict': f'gnomad_{data_type}.freq_index_dict',
                   'faf_index_dict': f'gnomad_{data_type}.faf_index_dict'})
    return ht


def main(args):
    hl.init(log='/frequency_annotations.log')

    data_source = args.data_source
    freeze = args.freeze
    mt = get_ukbb_data(data_source, freeze, meta_root='meta')
    logger.info(f'mt count before filtering out low quality samples: {mt.count()}')
    mt = mt.filter_cols(mt.meta.high_quality)
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    logger.info(f'mt count after filtering out low quality samples and their variants: {mt.count()}')

    if args.calculate_frequencies:
        logger.info('Calculating frequencies')
        ht, sample_table = generate_frequency_data(mt, args.downsampling, args.by_platform)

        write_temp_gcs(ht, var_annotations_ht_path(data_source, freeze, 'ukb_freq'), args.overwrite)
        if args.downsampling:
            sample_table.write(sample_annotations_table_path(data_type, 'downsampling'), args.overwrite)

    if args.join_gnomad:
        ht = hl.read_table(var_annotations_ht_path(data_source, freeze, 'ukb_freq'))

        logger.info('Joining UKBB ht to gnomAD exomes and genomes liftover hts')
        ht = join_gnomad(ht, 'exomes')
        ht = join_gnomad(ht, 'genomes')
        write_temp_gcs(ht, var_annotations_ht_path(data_source, freeze, 'join_freq'), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)
    parser.add_argument('--downsampling', help='Also calculate downsampling frequency data', action='store_true')
    parser.add_argument('--calculate_frequencies', help='Calculate most frequency data', action='store_true')
    parser.add_argument('--by_platform', help='Also calculate frequencies by platform', action='store_true')
    parser.add_argument('-j', '--join_gnomad', help='Join table with gnomAD tables to get gnomAD frequencies', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
