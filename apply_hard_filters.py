import argparse
from call_sex import *
from gnomad_hail import *
import logging
from resources import * 


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("apply_hard_filters")
logger.setLevel(logging.INFO)


def flag_ambiguous(mt: hl.MatrixTable, sex_ht: hl.Table) -> hl.MatrixTable:
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


def flag_low_qual(mt: hl.MatrixTable, min_callrate: float = 0.99, min_depth: float = 20.0) -> hl.MatrixTable:
    """
    Annotates samples with low callrate and depth
    :param MatrixTable mt: MatrixTable containing samples to be filtered
    :param float min_callrate: Callrate threshold to be used to filter samples; default is 0.99
    :return: MatrixTable with low quality samples flagged
    :rtype: MatrixTable
    """

    # the default cutoffs were set visually using plots:
    # p = hl.plot.histogram(mt.sample_qc.dp_stats.mean, range=(10,120), legend='Mean Sample DP')
    # p = hl.plot.histogram(mt.sample_qc.call_rate, range=(0.991, 0.997), legend='Mean Sample Callrate')
    mt = mt.annotate_cols(low_call_rate = hl.cond(mt.sample_qc.call_rate < min_callrate, True, False))
    mt = mt.annotate_cols(low_dp = hl.cond(mt.sample_qc.dp_stats.mean < min_depth, True, False))
    return mt
    

def main(args):

    hl.init()

    datasource = args.datasource
    if args.freeze:
        freeze = args.freeze
    else:
        freeze = CURRENT_FREEZE

    logger.info('Reading in mt with adj, variant_qc, and sample_qc annotations')
    mt = hl.read_matrix_table(adj_mt_path(datasource, freeze))

    logger.info('Reading in qc mt')
    qc_mt = hl.read_matrix_table(qc_mt_path(datasource, freeze, False))

    # TODO - make common function for returning reference genome + call that instead of this hail command
    logger.info('Getting build of qc mt')
    reference = qc_mt.locus.dtype.reference_genome.name

    logger.info('Imputing sex (using call_sex.py) on qc mt')
    sex_ht = impute_sex(qc_mt, reference, f'{sample_qc_prefix}/{datasource}.freeze_{freeze}/sex_check')
    
    logger.info('Flag samples with ambiguous sex or sex aneuploidy in raw mt')
    mt = flag_ambiguous(mt, sex_ht)

    logger.info('Flag low quality samples (samples with low callrate and depth) in raw mt')
    mt = flag_low_qual(mt, args.callrate, args.depth)

    logger.info('Writing raw mt with annotations')
    mt.write(hard_filters_mt_path(datasource, freeze), overwrite=args.overwrite)
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
