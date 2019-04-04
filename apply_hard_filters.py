import argparse
from call_sex import *
from gnomad_hail.utils.generic import *
from gnomad_hail.utils.sample_qc import add_filters_expr
import logging
from resources import * 


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("apply_hard_filters")
logger.setLevel(logging.INFO)


def apply_hard_filters_expr(mt: hl.MatrixTable, min_callrate: float = 0.99, min_depth: float = 20.0) -> hl.MatrixTable:
    """
    Creates hard filters expression and annotates mt with expression (creates hard_filters column)
    :param MatrixTable mt: MatrixTable to be annotated 
    :param float min_callrate: Callrate threshold to be used to filter samples; default is 0.99
    :param float min_depth: Depth threshold to be used to filter samples; default is 20.0
    :return: MatrixTable with updated filter field
    :rtype: MatrixTable
    """
    
    # the default coverage/depth cutoffs were set visually using plots:
    # p = hl.plot.histogram(mt.sample_qc.dp_stats.mean, range=(10,120), legend='Mean Sample DP')
    # p = hl.plot.histogram(mt.sample_qc.call_rate, range=(0.991, 0.997), legend='Mean Sample Callrate')

    hard_filters = {
        # we don't have contamination/chimera for regeneron vcf
        #'contamination': ht.freemix > 0.05,
        #'chimera': ht.pct_chimeras > 0.05,
        'low_callrate': mt.sample_qc.call_rate < min_callrate,
        'ambiguous_sex': mt.sex == 'ambiguous_sex',
        'sex_aneuploidy': mt.sex == 'sex_aneuploidy',
        'low_coverage': mt.sample_qc.dp_stats.mean < min_depth
    }

    mt = mt.annotate_cols(hard_filters = add_filters_expr(hard_filters, None))

    return mt


def annotate_sex(mt: hl.MatrixTable, sex_ht: hl.Table) -> hl.MatrixTable:
    """
    Annotates mt with imputed sex calculated using call_sex.py
    :param MatrixTable mt: MatrixTable containing samples to be annotated and filtered
    :param Table sex_ht: Table containing sex annotations from Kristen's code call_sex.py
    :return: MatrixTable with ambiguous sex/sex aneuploidies flagged
    :rtype: MatrixTable
    """

     # s	is_female	f_stat	n_called	expected_homs	observed_homs	sex	y_cov	twenty_cov	normalized_y_coverage
    sex_colnames = ['f_stat', 'is_female', 'sex']
    sex_ht = sex_ht.select(*sex_colnames)
    mt = mt.annotate_cols(**sex_ht[mt.col_key])
    return mt


def main(args):

    # hl.init() isn't necessary? according to Laurent
    #hl.init()

    datasource = args.datasource
    if args.freeze:
        freeze = args.freeze
    else:
        freeze = CURRENT_FREEZE

    logger.info('Reading in mt with adj, variant_qc, and sample_qc annotations')
    mt = hl.read_matrix_table(adj_mt_path(datasource, freeze))

    logger.info('Reading in qc mt')
    qc_mt = hl.read_matrix_table(qc_mt_path(datasource, freeze, False))

    logger.info('Adding variant_qc to qc mt')
    # add hail's variant qc to qc mt - this is to filter on frequency for imputing sex
    qc_mt = hl.variant_qc(mt)

    logger.info('Getting build of qc mt')
    reference = get_reference_genome(qc_mt.locus).name

    logger.info('Imputing sex (using call_sex.py) on qc mt')
    sex_ht = impute_sex(qc_mt, reference, f'{sample_qc_prefix}/{datasource}.freeze_{freeze}/sex_check')

    logger.info('Annotate mt with sex information')
    mt = annotate_sex(mt, sex_ht)

    logger.info('Adding hard filters to mt')
    mt = apply_hard_filters_expr(mt)
    
    logger.info('Writing raw mt with annotations')
    mt = mt.checkpoint(hard_filters_mt_path(datasource, freeze), overwrite = args.overwrite)

    logger.info('Checking number of samples flagged with hard filters')
    ht = mt.cols()
    ht = ht.explode(ht.hard_filters)
    filters = ht.aggregate(hl.agg.counter(ht.hard_filters))
    for filt in filters:
        logger.info(f'Samples flagged due to {filt}: {filters[filt]}')
    logger.info('Complete')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This script applies hard filters to UKBB data')
    parser.add_argument('-s', '--datasource', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Current freeze #', type=int)
    parser.add_argument('-c', '--callrate', help='Minimum callrate', default=0.99)
    parser.add_argument('-d', '--depth', help='Minimum depth', default=20)
    parser.add_argument('-o', '--overwrite', help='Overwrite pre-existing data', action='store_true', default=True)
    args = parser.parse_args()

    main(args)
