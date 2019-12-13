import hail as hl
from gnomad_hail import *
from gnomad_hail.utils.generic import *
from gnomad_hail.utils.sparse_mt import *
from ukbb_qc.resources import *

logging.basicConfig(format='%(asctime)s (%(name)s %(lineno)s): %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('densify_sparse_mt')
logger.setLevel(logging.INFO)


def rep_on_read(path, n_desired):
     mt = hl.read_matrix_table(path)
     intervals = mt._calculate_new_partitions(n_desired)
     return hl.read_matrix_table(path, _intervals=intervals)


def main(args):

    hl.init(default_reference='GRCh38', log='/densify_sparse_mt.log')    
    freeze = args.freeze
    n = args.n_partitions
    out_path = f'gs://broad-ukbb/broad.freeze_{freeze}/temp/' 

    logger.info('Reading in input mt (raw sparse mt)')
    mt = hl.read_matrix_table(args.input)
    logger.info(f'Sparse mt count: {mt.count()}')

    # subset of 10 samples
    samples = ['UKB_1000026_233144425', 'UKB_1000031_236163960', 'UKB_1000083_232997991', 'UKB_1000122_232973530', 'UKB_1000139_232977596', 'UKB_1000218_232941722', 'UKB_1000276_236163521', 'UKB_1000447_232921119', 'UKB_1000486_233000970', 'UKB_1000494_236259219']

    logger.info('Filtering to chr20 and sample subset')
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval('chr20')])
    mt = mt.filter_cols(hl.literal(samples).contains(mt.s))
    logger.info(f'Count after filtration: {mt.count()}')

    if args.before_densify:
        out_path = out_path + 'dense_subset.mt'
        mt = mt.repartition(n)
    else:
        out_path = out_path + 'temp_dense_subset.mt'
        
    logger.info('Densifying mt')
    mt = hl.experimental.densify(mt)

    logger.info('Filtering out lines that are only reference')
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)
    logger.info(f'Count after filtration: {mt.count()}')

    logger.info('Writing out mt')
    mt.write(out_path, args.overwrite)

    if args.on_read:
        logger.info('Reading in densified sparse mt and repartitioning')
        mt = rep_on_read(out_path, n)

        logger.info('Filtering to calling intervals')
        capture_ht = hl.read_table(capture_ht_path('broad', freeze))
        mt = mt.filter_rows(hl.is_defined(capture_ht[mt.locus])) 
        logger.info(f'Count after filtration: {mt.count()}')

        logger.info('Writing out mt')
        out_path = out_path.replace('temp_', '')
        mt.write(out_path, args.overwrite)
   

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Path to input mt', required=True)
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)
    parser.add_argument('-n', '--n_partitions', help='Number of desired partitions', type=int)
    parser.add_argument('--before_densify', help='Repartition before writing', action='store_true')
    parser.add_argument('--on_read', help='Repartition on read', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
