from gnomad_hail import *
from ukbb_qc.resources import *


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("import_vcfs")
logger.setLevel(logging.INFO)


def main(args):
    hl.init(log='/import_vcfs.log')

    data_source = args.data_source
    freeze = args.freeze

    logger.info('Importing VCFs...')
    mt = hl.import_vcf(args.import_regex, force_bgz=True, reference_genome='GRCh38').naive_coalesce(args.num_partitions)
    mt = mt.checkpoint(raw_mt_path(data_source, freeze, is_temp=True), overwrite=args.overwrite)
    mt.summarize(entries=False)

    logger.info('Dropping missing variant annotations...')
    # TODO: identify and drop highly missing annotations automatically
    drop_annots = [x for x in mt['info'] if x.startswith('AS_')]
    drop_annots.extend(['SB_TABLE', 'RAW_MQandDP', 'QUALapprox', 'InbreedingCoeff', 'ExcessHet', 'END'])
    logger.info(f'Dropping the following empty variant annotations: {drop_annots}')
    mt = mt.annotate_rows(info=mt.info.drop(*drop_annots))
    mt = mt.checkpoint(raw_mt_path(data_source, freeze), overwrite=args.overwrite)
    mt.summarize(entries=False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--import_regex', help='Regex string containing VCF(s) for import', type=str)
    parser.add_argument('-n', '--num_partitions', help='Number of partitions to be used for raw MT', type=int)
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE)

    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite pre-existing data', action='store_true')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
