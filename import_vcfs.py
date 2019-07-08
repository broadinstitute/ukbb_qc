from gnomad_hail import *
from ukbb_qc.resources import *
from ukbb_qc.sanity_checks import *


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("import_vcfs")
logger.setLevel(logging.INFO)


def main(args):
    hl.init(log='/import_vcfs.log')

    data_source = args.data_source
    freeze = args.freeze

    if args.import_vcfs:
        logger.info('Importing genotyped VCFs for inspection...')
        mt = hl.import_vcf(args.import_regex, force_bgz=True, reference_genome='GRCh38').naive_coalesce(args.num_partitions)
        mt = mt.checkpoint(raw_mt_path(data_source, freeze, is_temp=True), overwrite=args.overwrite)
        mt.summarize(entries=False)

        logger.info('Running sanity checks on raw callset...')
        sample_qc_path = sample_qc_prefix(datasource, freeze)
        var_summary = summarize_mt(mt)
        logger.info('\nVariant summary: {}'.format(var_summary))

        ht = hl.import_table(sample_list_path('regeneron', freeze), no_header=True)
        if sample_check(mt, ht, sample_qc_path):
            logger.warning('\nSample mismatch between dataset and list of expected samples')
        else:
            logger.info('\nAll samples in dataset match provided list of samples')
        # NOTE: stop here to inspect summaries, missing variant annotations, and sample IDs -- there were problems with Broad callset

    if args.import_vqsr:
        logger.info('Importing VQSR annotations...')
        mt = hl.import_vcf(args.import_vqsr_regex, force_bgz=True, reference_genome='GRCh38').naive_coalesce(args.num_vqsr_partitions)
        ht = mt.rows()
        ht = hl.split_multi_hts(ht).checkpoint(var_annotations_ht_path(data_source, freeze, 'vqsr'), overwrite=args.overwrite)
        row_count = ht.count()
        logger.info(f'Found {row_count} variants with VQSR annotations')

    if args.write_raw_mt:
        mt = hl.read_matrix_table(raw_mt_path(data_source, freeze, is_temp=True))
        mt = mt.key_cols_by(s=mt.s[1:])
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
    parser.add_argument('--import_vcfs', help='Import and write temporary raw MT', action='store_true')
    parser.add_argument('--import_vqsr', help='Import and write VQSR annotations to HT', action='store_true')
    parser.add_argument('--write_raw_mt', help='Trim temporary raw MT into final raw MT', action='store_true')

    parser.add_argument('-i', '--import_regex', help='Regex string containing VCF(s) for import', type=str)
    parser.add_argument('-v', '--import_vqsr_regex', help='Regex string containing VQSR VCF(s) for import', type=str)
    parser.add_argument('-n', '--num_partitions', help='Number of partitions to be used for raw MT', type=int)
    parser.add_argument('-m', '--num_vqsr_partitions', help='Number of partitions to be used for VQSR HT', type=int)
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE)

    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite pre-existing data', action='store_true')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
