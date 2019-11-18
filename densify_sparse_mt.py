import hail as hl
from gnomad_hail import *
from gnomad_hail.utils.generic import *
from gnomad_hail.utils.sparse_mt import *
from ukbb_qc.resources import *

logging.basicConfig(format='%(asctime)s (%(name)s %(lineno)s): %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('densify_sparse_mt')
logger.setLevel(logging.INFO)


def main(args):

    hl.init(default_reference='GRCh38', log='/densify_sparse_mt.log')    
    freeze = args.freeze

    logger.info('Reading in input mt (raw sparse mt)')
    #mt = get_ukbb_data('broad', freeze, adj=False, split=False, raw=True)
    mt = hl.read_matrix_table(args.input)
    logger.info(f'Sparse mt count: {mt.count()}')

    logger.info('Densifying mt')
    mt = hl.experimental.densify(mt)

    logger.info('Filtering out lines that are only reference')
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)
    logger.info(f'Count after filtration: {mt.count()}')

    logger.info('Writing out mt')
    mt.write(raw_mt_path('broad', freeze), args.overwrite)
   

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Path to input mt', required=True)
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)
    parser.add_argument('--test', help='Filter to chr20 to test', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
