from gnomad_hail import *
from gnomad_qc.annotations.generate_frequency_data import generate_downsamplings_cumulative
from ukbb_qc.resources import * 
from collections import Counter


DOWNSAMPLINGS = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000,
                 55000, 60000, 65000, 70000, 75000, 80000, 85000, 90000, 95000, 100000, 110000, 120000]
POPS_TO_REMOVE_FOR_POPMAX = ['asj', 'fin', 'oth']
F_CUTOFF = 0.05


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("generate_frequency_data")
logger.setLevel(logging.INFO)


def add_faf_expr(freq: hl.expr.ArrayExpression, freq_meta: hl.expr.ArrayExpression, locus: hl.expr.LocusExpression, populations: Set[str]) -> hl.expr.ArrayExpression:
    """
    Calculates popmax (add an additional entry into freq with popmax: pop)

    :param ArrayExpression freq: ArrayExpression of Structs with ['ac', 'an', 'hom']
    :param ArrayExpression freq_meta: ArrayExpression of meta dictionaries corresponding to freq
    :param LocusExpression locus: LocusExpression
    :param set of str populations: Set of populations over which to calculate popmax
    :return: Frequency data with annotated popmax
    :rtype: ArrayExpression
    """
    pops_to_use = hl.literal(populations)
    freq = hl.map(lambda x: x[0].annotate(meta=x[1]), hl.zip(freq, freq_meta))
    freqs_to_use = hl.filter(lambda f:
                             ((f.meta.size() == 1) & (f.meta.get('group') == 'adj')) |
                             ((f.meta.size() == 2) & (f.meta.get('group') == 'adj') & pops_to_use.contains(f.meta.get('pop'))) |
                             (~locus.in_autosome_or_par() & (
                                     ((f.meta.size() == 2) & (f.meta.get('group') == 'adj') & f.meta.contains('sex')) |
                                     ((f.meta.size() == 3) & (f.meta.get('group') == 'adj') & pops_to_use.contains(f.meta.get('pop')) & f.meta.contains('sex')))),
                             freq)
    return freqs_to_use.map(lambda f: hl.struct(
        meta=f.meta,
        faf95=hl.experimental.filtering_allele_frequency(f.AC, f.AN, 0.95),
        faf99=hl.experimental.filtering_allele_frequency(f.AC, f.AN, 0.99)
    ))


def generate_frequency_data(mt: hl.MatrixTable, calculate_downsampling: bool = False,
                            calculate_by_platform: bool = False,
                            pop_field: str) -> Tuple[hl.Table, hl.Table]:
    """
    :param MatrixTable mt: Input MatrixTable
    :param bool calculate_downsampling: Calculate frequencies for downsampled data
    :param bool calculate_by_platform: Calculate frequencies for PCR-free data
    :param str pop_field: Name of population field ('hybrid_pop' or 'gnomad_pc_project_pop')
    """
    if calculate_downsampling:
        mt, downsamplings = generate_downsamplings_cumulative(mt)
        logger.info(f'Got {len(downsamplings)} downsamplings: {downsamplings}')
    cut_dict = {'pop': hl.agg.filter(hl.is_defined(mt.meta[f'{pop_field}']), hl.agg.counter(mt.meta[f'{pop_field}'])),
                'sex': hl.agg.filter(hl.is_defined(mt.meta.sex), hl.agg.collect_as_set(mt.meta.sex)),
                'subpop': hl.agg.filter(hl.is_defined(mt.meta.subpop) & hl.is_defined(mt.meta[f'{pop_field}']),
                                        hl.agg.collect_as_set(hl.struct(subpop=mt.meta.subpop, pop=mt.meta[f'{pop_field}'])))
                }
    logger.info(cut_dict)

    if calculate_by_platform:
        cut_dict['platform'] = hl.agg.filter(hl.is_defined(mt.meta.qc_platform),
                                             hl.agg.collect_as_set(mt.meta.qc_platform))
    cut_data = mt.aggregate_cols(hl.struct(**cut_dict))

    sample_group_filters = [({}, True)]
    sample_group_filters.extend([
        ({'pop': pop}, mt.meta[f'{pop_field}'] == pop) for pop in cut_data[f'{pop_field}']
    ] + [
        ({'sex': sex}, mt.meta.sex == sex) for sex in cut_data.sex
    ] + [
        ({'pop': pop, 'sex': sex}, (mt.meta.sex == sex) & (mt.meta[f'{pop_field}'] == pop))
        for sex in cut_data.sex for pop in cut_data[f'{pop_field}']
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
             (mt.downsampling[f'{pop_field}']_idx < ds) & (mt.meta[f'{pop_field}'] == pop))
            for ds in downsamplings for pop, pop_count in cut_data[f'{pop_field}'].items() if ds <= pop_count
        ])
    mt = mt.select_cols(group_membership=tuple(x[1] for x in sample_group_filters), project_id=mt.meta.project_id, age=mt.meta.age)
    mt = mt.select_rows()

    frequency_expression = []
    meta_expressions = []
    for i in range(len(sample_group_filters)):
        subgroup_dict = sample_group_filters[i][0]
        subgroup_dict['group'] = 'adj'
        call_stats = hl.agg.filter(mt.group_membership[i] & mt.adj, hl.agg.call_stats(mt.GT, mt.alleles))
        call_stats_bind = hl.bind(lambda cs: cs.annotate(
            AC=cs.AC[1], AF=cs.AF[1], homozygote_count=cs.homozygote_count[1]
        ), call_stats)
        frequency_expression.append(call_stats_bind)
        meta_expressions.append(subgroup_dict)

    raw_stats = hl.agg.call_stats(mt.GT, mt.alleles)
    raw_stats_bind = hl.bind(lambda cs: cs.annotate(
        AC=cs.AC[1], AF=cs.AF[1], homozygote_count=cs.homozygote_count[1]
    ), raw_stats)
    frequency_expression.insert(1, raw_stats_bind)
    meta_expressions.insert(1, {'group': 'raw'})

    logger.info(f'Calculating {len(frequency_expression)} aggregators...')
    global_expression = {
        'freq_meta': meta_expressions
    }
    mt = mt.annotate_rows(freq=frequency_expression,
                          age_hist_het=hl.agg.filter(mt.adj & mt.GT.is_het(), hl.agg.hist(mt.age, 30, 80, 10)),
                          age_hist_hom=hl.agg.filter(mt.adj & mt.GT.is_hom_var(), hl.agg.hist(mt.age, 30, 80, 10)))
    if calculate_downsampling: global_expression['downsamplings'] = downsamplings
    mt = mt.annotate_globals(**global_expression)
    sample_data = mt.cols()

    pops = set(cut_data[f'{pop_field}'].keys())
    [pops.discard(x) for x in POPS_TO_REMOVE_FOR_POPMAX]

    mt = mt.annotate_rows(popmax=add_popmax_expr(mt.freq, mt.freq_meta, populations=pops),
                          faf=add_faf_expr(mt.freq, mt.freq_meta, mt.locus, populations=pops))
    mt = get_projectmax(mt, mt.project_id)

    return mt.rows(), sample_data


def join_gnomad(ht: hl.Table, data_type: str) -> hl.Table:
    """
    Joins UKBB ht to gnomAD ht and adds gnomAD freq, popmax, and faf as annotation

    :param Table ht: Input UKBB ht
    :param str data_type: One of exomes or genomes
    :return: UKBB ht with gnomAD frequency information added as annotation
    :rtype: Table
    """
    gnomad_ht = hl.read_table(get_gnomad_liftover_data_path(f'{data_type}', '2.1.1')).select('freq', 'popmax', 'faf')
    return ht.annotate(hl.format('gnomad_%s', data_type)=gnomad_ht[ht.locus, ht.alleles])


def main(args):
    hl.init(log='/frequency_annotations.log')

    data_source = args.data_source
    freeze = args.freeze
    mt = get_ukbb_data(data_source, freeze, meta_root='meta')
    mt.describe()

    # get desired  population field from metadata
    if args.gnomad:
        pop_field = 'gnomad_pc_project_pop'
    else:
        pop_field = 'hybrid_pop'

    if args.calculate_frequencies:
        ht, sample_table = generate_frequency_data(mt, args.downsampling, args.by_platform, pop_field)

        write_temp_gcs(ht, annotations_ht_path(data_source, freeze, 'ukb_freq'), args.overwrite)
        if args.downsampling:
            sample_table.write(sample_annotations_table_path(data_type, 'downsampling'), args.overwrite)

    if args.join_gnomad:
        exomes_ht = hl.read_table(get_gnomad_liftover_data_path('exomes', '2.1.1')).select('freq', 'popmax', 'faf')
        genomes_ht = hl.read_table(get_gnomad_liftover_data_path('genomes', '2.1.1')).select('freq', 'popmax', 'faf')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)
    parser.add_argument('-g', '--gnomad', help='Calculate frequencies using gnomAD pc project pop', action='store_true')
    parser.add_argument('--downsampling', help='Also calculate downsampling frequency data', action='store_true')
    parser.add_argument('--calculate_frequencies', help='Calculate most frequency data', action='store_true')
    parser.add_argument('--by_platform', help='Also calculate frequencies by platform', action='store_true')
    parser.add_argument('-j', '--join_gnomad', help='Join table with gnomAD tables to get gnomAD frequencies', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified')

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
