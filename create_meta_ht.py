import argparse
import logging
from resources import * 


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("create_meta_ht")
logger.setLevel(logging.INFO)


def sample_comparison(ht1, ht2) -> bool:
    """
    Checks if two sample counts are the same

    :param Table ht1: First Table to be checked
    :param Table ht2: Second Table to be checked
    :return: Whether the sample counts are the same
    :rtype: bool
    """

    s_count1 = ht1.count()
    s_count2 = ht2.count()
    logger.info(f'{s_count1} samples in left table; {s_count2} samples in right table')

    return s_count1 == s_count2


def table_join(
    left_ht: hl.Table, left_key: str, right_ht: hl.Table, right_key: str, join_type: str
    ) -> hl.Table:
    """
    Joins two tables and returns joined table. Also prints warning if sample counts are not the same.

    :param Table left_ht: left Table to be joined
    :param str left_key: key for left Table
    :param Table right_ht: right Table to be joined
    :param str right_key: key for right Table
    :param str join_type: how to join the tables (left, right, inner, outer)
    :return: joined Table
    :rtype: Table
    """


    sample_count_match = sample_comparison(left_ht, right_ht)

    if not sample_count_match:
        logger.warning('Sample counts in left and right tables do not match')

    return left_ht.key_by(left_key).join(right_ht.key_by(right_key), how = join_type)


def main(args):

    data_source = args.data_source
    if args.freeze:
        freeze = args.freeze
    else:
        freeze = CURRENT_FREEZE

    logger.info('Reading hard filters ht')
    left_ht = hl.read_table(hard_filters_ht_path(data_source, freeze))

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
