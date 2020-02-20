import argparse
import hail as hl
import logging
from gnomad_hail.utils.sample_qc import add_filters_expr
from ukbb_qc.resources.basics import get_ukbb_data
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.utils.utils import annotate_interval_qc_filter 


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("apply_hard_filters")
logger.setLevel(logging.INFO)


def apply_hard_filters_expr(ht: hl.Table, min_callrate: float, min_depth: float) -> hl.Table:
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

    logger.info('Callrate cutoff for hard filters: {}'.format(min_callrate))
    logger.info('Depth cutoff for hard filters: {}'.format(min_depth))

    hard_filters = {
        # we don't have contamination/chimera for regeneron vcf
        #'contamination': ht.freemix > 0.05,
        #'chimera': ht.pct_chimeras > 0.05,
        'low_callrate': ht.raw_sample_qc.call_rate < min_callrate,
        'ambiguous_sex': ht.sex == 'ambiguous_sex',
        'sex_aneuploidy': ht.sex == 'sex_aneuploidy',
        'low_coverage': ht.raw_sample_qc.dp_stats.mean < min_depth
    }

    ht = ht.annotate(hard_filters=add_filters_expr(hard_filters, None))
    return ht


def main(args):
    hl.init(log='/apply_hard_filters.log', default_reference='GRCh38')

    data_source = args.data_source
    freeze = args.freeze
    temp_ht_path = f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/interval_qc_sample_qc.ht'

    if args.intervals:
        logger.info('Adding sex imputation annotations...')
        logger.info('Filtering to well covered intervals...')
        mt = get_ukbb_data(data_source, freeze, split=False, raw=True)
        mt = annotate_interval_qc_filter(data_source, freeze, mt, autosomes_only=True)
        mt = mt.filter_rows(mt.interval_qc_pass)

        logger.info('Adding sex imputation annotations...')
        ht = hl.read_table(sex_ht_path(data_source, freeze))
        mt = mt.annotate_cols(**ht[mt.col_key])

        logger.info('Computing raw sample QC metrics...')
        mt = mt.annotate_entries(GT=hl.experimental.lgt_to_gt(mt.LGT, mt.LA))
        mt = mt.key_rows_by('locus', 'alleles')
        mt = hl.sample_qc(mt)
        ht = mt.transmute_cols(raw_sample_qc=mt.sample_qc.select('call_rate', 'gq_stats', 'dp_stats')).cols()
        ht.write(temp_ht_path, overwrite=args.overwrite)

    ht = hl.read_table(temp_ht_path)
    ht = apply_hard_filters_expr(ht, args.callrate, args.depth)
    
    logger.info('Writing out hard filters HT...')
    ht = ht.naive_coalesce(args.n_partitions)
    ht = ht.checkpoint(hard_filters_ht_path(data_source, freeze), overwrite=args.overwrite)
    
    logger.info('Checking number of samples flagged with hard filters...')
    ht = ht.explode(ht.hard_filters)
    filters = ht.aggregate(hl.agg.counter(ht.hard_filters))
    for filt in filters:
        logger.info(f'Samples flagged due to {filt}: {filters[filt]}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This script applies hard filters to UKBB data')
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)
    parser.add_argument('-c', '--callrate', help='Minimum callrate', default=0.99, type=float)
    parser.add_argument('-d', '--depth', help='Minimum depth', default=20.0, type=float)
    parser.add_argument('--n_partitions', help='Desired number of partitions', default=15000, type=int)
    parser.add_argument('-o', '--overwrite', help='Overwrite pre-existing data', action='store_true', default=True)
    parser.add_argument('-i', '--intervals', help='Filter to intervals to check distributions', action='store_true')
    args = parser.parse_args()

    main(args)
