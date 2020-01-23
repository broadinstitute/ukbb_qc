import hail as hl
from gnomad_hail import *
from gnomad_hail.utils.generic import rep_on_read
from ukbb_qc.resources import *


logging.basicConfig(format='%(asctime)s (%(name)s %(lineno)s): %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('densify_sparse_mt')
logger.setLevel(logging.INFO)


def main(args):

    hl.init(default_reference='GRCh38', log='/densify_sparse_mt.log')    
    freeze = args.freeze

    # NOTE: for tranche 2, autoscaled to 1000 preemptibles with 400GB disk on workers and master 
    logger.info('Reading in capture ht')
    capture_ht = hl.read_table(capture_ht_path('broad', freeze))

    # NOTE: set to 25000 partitions for tranche 2
    logger.info('Reading in input mt (raw sparse mt)')
    mt = rep_on_read(args.input, args.n_partitions)
    
    logger.info('Densifying mt')
    mt = hl.experimental.densify(mt)

    # NOTE: Do NOT do any counts between densifying and writing
    # NOTE: Used gs://broad-ukbb/resources/xgen_plus_spikein.Homo_sapiens_assembly38.targets.pad50.intervals.ht on tranche 2
    logger.info('Filtering out lines that are only reference or not covered in capture intervals')
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)
    mt = mt.filter_rows(hl.is_defined(capture_ht[mt.locus]))

    # NOTE: this job took approximately 24 hours for tranche 2
    logger.info('Writing out mt')
    mt.write(raw_mt_path('broad', freeze), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Path to input mt', required=True)
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)
    parser.add_argument('-n', '--n_partitions', help='Number of desired partitions', type=int, default=25000)
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
