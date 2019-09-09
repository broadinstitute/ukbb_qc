from gnomad_hail import *
from gnomad_hail.utils.sample_qc import *
from ukbb_qc.resources import *
import hail as hl
import argparse


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("relatedness")
logger.setLevel(logging.INFO)


def get_ped(relatedness_ht: hl.Table, dups_ht: hl.Table, sex_ht: hl.Table) -> hl.Pedigree:
    relatedness_ht = relatedness_ht.key_by(  # Needed for familial inference at this point -- should be generalized
        i=relatedness_ht.i.s,
        j=relatedness_ht.j.s
    )
    dups_to_remove = dups_ht.aggregate(hl.agg.filter(dups_ht.duplicate == True, hl.agg.collect_as_set(dups_ht.s)))
    logger.info(f"Removing {len(dups_to_remove)} duplicates from family creation.")
    sex = {row.s: row.is_female for row in sex_ht.to_pandas().itertuples()}

    # check for empty set
    if len(dups_to_remove) == 0:
        dups_to_remove = hl.empty_set(hl.tstr)
    ped = infer_families(relatedness_ht, sex, dups_to_remove)
    logger.info(f"Found {len(ped.complete_trios())} complete trios.")
    return ped


# Note: Will need to change this to work with gnomAD, we will want to default to gnomAD ranking, unless tied pair is gnomAD, UKBB pair
def rank_related_samples(
        relatedness_ht: hl.Table,
        qc_ht: hl.Table
) -> Tuple[hl.Table, Callable[[hl.expr.Expression, hl.expr.Expression], hl.expr.NumericExpression]]:

    def annotate_related_pairs(related_pairs: hl.Table, index_col: str) -> hl.Table:
        related_pairs = related_pairs.key_by(**related_pairs[index_col])
        return related_pairs.annotate(
            **{
                index_col: related_pairs[index_col].annotate(
                    dp_mean=hl.or_else(qc_ht[related_pairs.key].sample_qc.dp_stats.mean, -1.0)
                )
            }
        ).key_by()

    relatedness_ht = annotate_related_pairs(relatedness_ht, "i")
    relatedness_ht = annotate_related_pairs(relatedness_ht, "j")

    def tie_breaker(l, r):
        return (l.dp_mean - r.dp_mean)

    return relatedness_ht, tie_breaker


def main(args):
    hl.init(log='/relatedness.log', tmp_dir='hdfs:///pc_relate.tmp/')

    data_source = args.data_source
    freeze = args.freeze

    if not args.skip_pc_relate:
        logger.info('Running PCA for PC-Relate...')
        pruned_qc_mt = remove_hard_filter_samples(data_source, freeze, 
                                            hl.read_matrix_table(qc_mt_path(data_source, freeze, ld_pruned=True)).unfilter_entries()
        eig, scores, _ = hl.hwe_normalized_pca(pruned_qc_mt.GT, k=10, compute_loadings=False)
        scores.write(relatedness_pca_scores_ht_path(data_source, freeze), args.overwrite)

        logger.info('Running PC-Relate...')
        # NOTE: This needs SSDs on your workers (for the temp files) and no preemptible workers while the BlockMatrix writes
        relatedness_ht = hl.pc_relate(pruned_qc_mt.GT, min_individual_maf=0.05, scores_expr=scores[pruned_qc_mt.col_key].scores,
                                      block_size=4096, min_kinship=args.min_emission_kinship, statistics='all')
        relatedness_ht.write(relatedness_ht_path(data_source, freeze), args.overwrite)

    if not args.skip_filter_dups:
        logger.info("Filtering duplicate samples...")
        sample_qc_ht = remove_hard_filter_samples(data_source, freeze, hl.read_table(qc_ht_path(data_source, freeze)))
        samples_rankings_ht = sample_qc_ht.select(rank=-1 * sample_qc_ht.sample_qc.dp_stats.mean)
        relatedness_ht = hl.read_table(relatedness_ht_path(data_source, freeze))
        if len(get_duplicated_samples(relatedness_ht)) > 0:
            dups_ht = filter_duplicate_samples(relatedness_ht, samples_rankings_ht)
            dups_ht.write(duplicates_ht_path(data_source, freeze, dup_sets=True), overwrite=args.overwrite)
            dups_ht = flatten_duplicate_samples_ht(dups_ht)
            dups_ht = sample_qc_ht.select(
                duplicate=hl.is_defined(dups_ht[sample_qc_ht.key].dup_filtered) & dups_ht[sample_qc_ht.key].dup_filtered,
            )
        else:
            dups_ht = sample_qc_ht.select(duplicate=False)
        dups_ht.write(duplicates_ht_path(data_source, freeze), overwrite=args.overwrite)

    if not args.skip_infer_families:
        logger.info("Inferring families")
        ped = get_ped(
            hl.read_table(relatedness_ht_path(data_source, freeze)),
            hl.read_table(duplicates_ht_path(data_source, freeze)),
            hl.read_table(sex_ht_path(data_source, freeze))
        )
        ped.write(inferred_ped_path(data_source, freeze))

    if not args.skip_filter_related_samples:
        logger.info("Filtering related samples")
        relatedness_ht = hl.read_table(relatedness_ht_path(data_source, freeze))
        relatedness_ht = relatedness_ht.filter(relatedness_ht.kin > args.min_filtering_kinship)
        related_pairs_ht, related_pairs_tie_breaker = rank_related_samples(
            relatedness_ht,
            hl.read_table(qc_ht_path(data_source, freeze))
        )

        related_samples_to_drop_ht = hl.maximal_independent_set(related_pairs_ht.i, related_pairs_ht.j,
                                                                keep=False, tie_breaker=related_pairs_tie_breaker)
        related_samples_to_drop_ht = related_samples_to_drop_ht.key_by()
        related_samples_to_drop_ht = related_samples_to_drop_ht.select(**related_samples_to_drop_ht.node)
        related_samples_to_drop_ht = related_samples_to_drop_ht.key_by('s')
        related_samples_to_drop_ht.write(related_drop_path(data_source, freeze), overwrite=args.overwrite)
        logger.info(f'{related_samples_to_drop_ht.count()} samples flagged in callset using maximal independent set')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)

    parser.add_argument('--skip_pc_relate', help='Skip running PC-relate on all samples. NOTE: This needs SSDs on your workers (for the temp files) and no pre-emptibles while the BlockMatrix writes',
                        action='store_true')
    parser.add_argument('--min_emission_kinship', help='Minimum kinship threshold for emitting a pair of samples in PC relate and filtering related individuals.',
                        default=0.05,
                        type=float)
    parser.add_argument('--min_filtering_kinship',
                        help='Minimum kinship threshold for filtering a pair of samples in PC relate and filtering related individuals. (Default = 0.08838835; 2nd degree relatives)',
                        default=0.08838835, type=float)
    parser.add_argument('--skip_filter_dups', help='Skip filtering duplicated samples', action='store_true')
    parser.add_argument('--skip_infer_families',
                        help='Skip extracting duplicate samples and infers families samples based on PC-relate results',
                        action='store_true')
    parser.add_argument('--skip_filter_related_samples',
                        help='Skip Filter related samples (based on the pairs present from the --run_pc_relate and using the --min_filtering_kinship value for that run)',
                        action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
