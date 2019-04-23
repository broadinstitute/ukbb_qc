import argparse
import logging
from resources import * 


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("create_meta_ht")
logger.setLevel(logging.INFO)


def main(args):

    data_source = args.data_source
    if args.freeze:
        freeze = args.freeze
    else:
        freeze = CURRENT_FREEZE

    logger.info('Reading in mt with hard filters and sex annotations')
    mt = hl.read_matrix_table(hard_filters_mt_path(data_source, freeze))
    mt = mt.annotate_cols(raw_sample_qc = mt.sample_qc)
    mt = mt.drop(mt.sample_qc)

    logger.info('Converting hard filters mt to ht')
    left_ht = mt.cols()

    logger.info('Reading in qc ht')
    right_ht = hl.read_table(qc_ht_path(data_source, freeze))

    logger.info('Joining hard filters ht and qc ht to start making meta ht')
    left_ht = table_join(left_ht, 's', right_ht, 's', 'outer')

    logger.info('Reading in platform PCA ht')
    right_ht = hl.read_table(platform_pca_results_ht_path(data_source, freeze))
    right_ht = right_ht.annotate(platform_pca_scores = right_ht.scores)
    right_ht = right_ht.drop('scores')

    logger.info('Joining platform PCA ht into current join (meta ht)')
    left_ht = table_join(left_ht, 's', right_ht, 's', 'outer')

    logger.info('Reading in population PC ht')
    right_ht = hl.read_table(ancestry_pc_project_scores_ht_path(data_source, freeze, None))

    logger.info('Joining population PC ht into current join')
    left_ht = table_join(left_ht, 's', right_ht, 's', 'outer')

    logger.info('Creating checkpoint')
    left_ht = left_ht.checkpoint(get_ht_checkpoint_path(data_source, freeze, 'intermediate_ht_join'), overwrite = True)

    logger.info('Reading in array sample concordance ht')
    right_ht = hl.read_table(array_sample_concordance_path(data_source, freeze))

    logger.info('Joining array sample concordance ht to current join')
    left_ht = left_ht.annotate(array_id = left_ht.s.split("_")[1])
    left_ht = table_join(left_ht, 'array_id', right_ht, 's', 'left')

    logger.info('Reading in dups ht')
    right_ht = hl.read_table(duplicates_ht_path(data_source, freeze, args.dup_sets))
    right_ht = right_ht.drop('qc_mt_params')

    logger.info('Joining dups ht to current join')
    left_ht = table_join(left_ht, 's', right_ht, 's', 'left')

    logger.info('Reading in related samples ht')
    related_samples_to_drop_ht = hl.read_table(related_drop_path(data_source, freeze))

    logger.info('Annotating meta ht with related samples')
    left_ht = left_ht.annotate(related_filter = hl.is_defined(related_samples_to_drop_ht[left_ht.s]))
   
    logger.info('Reading in outlier ht')
    right_ht = hl.read_table(platform_pop_outlier_ht_path(data_source, freeze)) 

    logger.info('Joining outlier ht to current join')
    left_ht = table_join(left_ht, 's', right_ht, 's', 'left')

    logger.info('Writing out meta ht')
    left_ht.write(meta_ht_path(data_source, freeze), overwrite = True)
    logger.info('Complete') 


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This script applies hard filters to UKBB data')
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Current freeze #', type=int)
    parser.add_argument('-d', '--dup_sets', help='Dup sets for duplicate ht table', action='store_true', default=False)
    args = parser.parse_args()

    main(args)
