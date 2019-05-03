import argparse
from gnomad_hail.utils.sample_qc import add_filters_expr
from ukbb_qc.call_sex import * 


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("apply_hard_filters")
logger.setLevel(logging.INFO)


def apply_hard_filters_expr(ht: hl.Table, min_callrate: float = 0.99, min_depth: float = 20.0) -> hl.Table:
    """
    Creates hard filters expression and annotates ht with expression (creates hard_filters column)
    :param Table ht: Table to be annotated 
    :param float min_callrate: Callrate threshold to be used to filter samples; default is 0.99
    :param float min_depth: Depth threshold to be used to filter samples; default is 20.0
    :return: Table with hard_filters column
    :rtype: Table
    """
    
    # the default coverage/depth cutoffs were set visually using plots:
    # p = hl.plot.histogram(mt.sample_qc.dp_stats.mean, range=(10,120), legend='Mean Sample DP')
    # p = hl.plot.histogram(mt.sample_qc.call_rate, range=(0.991, 0.997), legend='Mean Sample Callrate')

    hard_filters = {
        # we don't have contamination/chimera for regeneron vcf
        #'contamination': ht.freemix > 0.05,
        #'chimera': ht.pct_chimeras > 0.05,
        'low_callrate': ht.sample_qc.call_rate < min_callrate,
        'ambiguous_sex': ht.sex == 'ambiguous_sex',
        'sex_aneuploidy': ht.sex == 'sex_aneuploidy',
        'low_coverage': ht.sample_qc.dp_stats.mean < min_depth
    }

    ht = ht.annotate(hard_filters=add_filters_expr(hard_filters, None))
    return ht


def main(args):

    data_source = args.data_source
    freeze = args.freeze

    logger.info('Reading in raw MT...')
    mt = hl.read_matrix_table(raw_mt_path(data_source, freeze))

    logger.info('Computing raw sample QC metrics...')
    mt = hl.sample_qc(mt)
    mt = mt.transmute_cols(raw_sample_qc=mt.sample_qc)
    ht = mt.cols()
 
    logger.info('Annotating hard filters...')
    ht = apply_hard_filters_expr(ht)
    
    logger.info('Writing out hard filters HT...')
    ht = ht.checkpoint(hard_filters_ht_path(data_source, freeze), overwrite=args.overwrite)
    
    logger.info('Checking number of samples flagged with hard filters...')
    ht = ht.explode(ht.hard_filters)
    filters = ht.aggregate(hl.agg.counter(ht.hard_filters))
    for filt in filters:
        logger.info(f'Samples flagged due to {filt}: {filters[filt]}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This script applies hard filters to UKBB data')
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE)
    parser.add_argument('-c', '--callrate', help='Minimum callrate', default=0.99)
    parser.add_argument('-d', '--depth', help='Minimum depth', default=20)
    parser.add_argument('-o', '--overwrite', help='Overwrite pre-existing data', action='store_true', default=True)
    args = parser.parse_args()

    main(args)
