from resources import *
from gnomad_hail import *
from platform_pca import *
import hail as hl
import argparse


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("relatedness")
logger.setLevel(logging.INFO)


# Note: The following functions are slightly modified from Laurent's myoseq sample qc should probably move to a common location
def compute_qc_mt(mt: hl.MatrixTable, min_af: float = 0.001, min_callrate: float = 0.99) -> hl.MatrixTable:
    """
    Returns MatrixTable for sample QC purposes
    Default criteria: callrate > 0.99, AF > 0.001, SNPs only, bi-allelics only
    :param MatrixTable mt: Raw MatrixTable to be filtered
    :param float min_af: Minimum allele frequency for variant filtering
    :param float min_callrate: Minimum callrate for variant filtering
    :return: MatrixTable for sample QC purposes
    :rtype: MatrixTable
    """
    qc_mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]) &
                        (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > min_af) &
                        (hl.agg.fraction(hl.is_defined(mt.GT)) > min_callrate))
    qc_mt = qc_mt.annotate_globals(
        qc_mt_params=hl.struct(
            min_af=min_af,
            min_callrate=min_callrate
        )
    )
    qc_mt.annotate_cols(sample_callrate=hl.agg.fraction(hl.is_defined(qc_mt.GT)))
    return qc_mt


def ld_prune_qc_mt(qc_mt: hl.MatrixTable, ld_r2: float = 0.1) -> hl.MatrixTable:
    """
    Returns MatrixTable for sample QC purposes
    :param MatrixTable mt: Raw MatrixTable to be filtered
    :param float ld_r2: r2 to LD prune the QC matrix to
    :return: LD pruned MatrixTable for sample QC purposes
    :rtype: MatrixTable
    """
    pruned_ht = hl.ld_prune(qc_mt.GT, r2=ld_r2)
    qc_mt = qc_mt.filter_rows(hl.is_defined(pruned_ht[qc_mt.row_key]))
    qc_mt = qc_mt.annotate_globals(qc_mt_params=qc_mt.qc_mt_params.annotate(ld_r2=ld_r2))
    qc_mt.annotate_cols(sample_callrate=hl.agg.fraction(hl.is_defined(qc_mt.GT)))
    return qc_mt


# Note: everything below was copied exactly from Laurent's myoseq qc

def filter_dups(
        relatedness_ht: hl.Table,
        dup_ranking_func: Callable[..., Tuple[hl.Table, hl.expr.Expression]],
        dup_ranking_func_args: Optional[List]
):
    """
    Creates a HT with duplicated samples sets.
    Each row is indexed by the sample that is kept and also contains the set of duplicate samples that should be filtered.
    dup_ranking_func is a function to decide which duplicate to keep.
    It should take a table keyed by samples key and any number of additional arguments given through `dup_ranking_func_args`
    It should return a tuple containing:
        - the input table with any modifications (e.g. annotations) needed for ranking
        - A hail expression of a type that can be sorted giving the corresponding rank (where smaller is better)
    :param relatedness_ht: Input relatedness HT
    :param dup_ranking_func: Ranking function to pick amongst duplicates.
    :param dup_ranking_func_args: Optional additional arguments for `dup_ranking_func`
    :return:
    """
    logger.info("Getting duplicate samples")
    dups = get_duplicated_samples(relatedness_ht)
    logger.info(f"Found {len(dups)} duplicate sets.")
    dups_ht = hl.Table.parallelize([hl.struct(dup_set=i, dups=dups[i]) for i in range(0, len(dups))])
    dups_ht = dups_ht.explode(dups_ht.dups, name='_dup')
    dups_ht = dups_ht.key_by(**dups_ht._dup)
    if dup_ranking_func_args is None:
        dup_ranking_func_args = []
    dups_ht, rank_expr = dup_ranking_func(dups_ht, *dup_ranking_func_args)
    dups_cols = hl.bind(
        lambda x: hl.struct(
            kept=x[0],
            filtered=x[1:]
        ),
        hl.sorted(hl.agg.collect(hl.tuple([dups_ht._dup, rank_expr])), key=lambda x: x[1]).map(lambda x: x[0])
    )
    dups_ht = dups_ht.group_by(dups_ht.dup_set).aggregate(
        **dups_cols
    )

    dups_ht = dups_ht.key_by(**{f'{x}_kept': expr for x, expr in dups_ht.kept.items()}).drop('kept')
    return dups_ht


def get_ped(relatedness_ht: hl.Table, dups_ht: hl.Table, sex_ht: hl.Table) -> hl.Pedigree:
    relatedness_ht = relatedness_ht.key_by(  # Needed for familial inference at this point -- should be generalized
        i=relatedness_ht.i.s,
        j=relatedness_ht.j.s
    )
    dups_to_remove = dups_ht.aggregate(hl.agg.explode(lambda x: hl.agg.collect_as_set(x.s), dups_ht.filtered))
    logger.info(f"Removing {len(dups_to_remove)} duplicates from family creation.")
    sex = {row.s: row.is_female for row in sex_ht.to_pandas().itertuples()}
    ped = infer_families(relatedness_ht, sex, dups_to_remove)
    logger.info(f"Found {len(ped.complete_trios())} complete trios.")
    return ped


def rank_dup_samples(dups_ht: hl.Table, qc_ht: hl.Table) -> Tuple[hl.Table, hl.expr.Expression]:
    dups_ht = dups_ht.annotate(rank=-1 * qc_ht[dups_ht.key].sample_qc.dp_stats.mean)
    return dups_ht, dups_ht.rank


def get_related_samples_to_drop(
        relatedness_ht: hl.Table,
        min_filtering_kinship: float,
        rank_func: Callable[..., Tuple[hl.Table, Callable[[hl.expr.Expression, hl.expr.Expression], hl.expr.NumericExpression]]],
        rank_func_args: Optional[List]
) -> hl.Table:
    """
    Use the maximal independence function in Hail to intelligently prune clusters of related individuals, removing
    less desirable samples while maximizing the number of unrelated individuals kept in the sample set.
    :param relatedness_ht:
    :param min_filtering_kinship:
    :param rank_func:
    :param rank_func_args:
    :return:
    """
    # Define maximal independent set, using rank list
    related_pairs = relatedness_ht.filter(relatedness_ht.kin > min_filtering_kinship).key_by().select('i', 'j')
    n_related_samples = related_pairs.annotate(ij=[related_pairs.i, related_pairs.j]).explode('ij').key_by('ij').distinct().count()
    logger.info('{} samples with at least 2nd-degree relatedness found in callset'.format(n_related_samples))

    if rank_func_args is None:
        rank_func_args = []
    related_pairs, tie_breaker = rank_func(relatedness_ht, *rank_func_args)

    related_samples_to_drop_ranked = hl.maximal_independent_set(related_pairs.i, related_pairs.j,
                                                                keep=False, tie_breaker=tie_breaker)
    related_samples_to_drop_ranked = related_samples_to_drop_ranked.key_by()
    related_samples_to_drop_ranked = related_samples_to_drop_ranked.select(**related_samples_to_drop_ranked.node)
    return related_samples_to_drop_ranked.key_by(*relatedness_ht.i)


# Note: Will need to change this to work with gnomAD, we will want to default to gnomAD ranking, unless tied pair is gnomAD, UKBB pair
def rank_related_samples(
        related_pairs: hl.Table,
        qc_ht: hl.Table
) -> Tuple[hl.Table, Callable[[hl.expr.Expression, hl.expr.Expression], hl.expr.NumericExpression]]:

    def annotate_related_pairs(related_pairs: hl.Table, index_col: str) -> hl.Table:
        related_pairs = related_pairs.key_by(**related_pairs[index_col])
        related_pairs = related_pairs.filter(hl.is_missing(case_parents[related_pairs.key]))
        return related_pairs.annotate(
            **{
                index_col: related_pairs[index_col].annotate(
                    dp_mean=hl.or_else(qc_ht[related_pairs.key].sample_qc.dp_stats.mean, -1.0)
                )
            }
        ).key_by()

    related_pairs = annotate_related_pairs(related_pairs, "i")
    related_pairs = annotate_related_pairs(related_pairs, "j")

    def tie_breaker(l, r):
        return ( l.dp_mean - r.dp_mean)

    return related_pairs, tie_breaker


def main(args):
    hl.init(log='/relatedness.log', tmp_dir='hdfs:///pc_relate.tmp/')

    data_source = args.data_source
    freeze = args.freeze

    # Note: should we move this part and the one below it into the script that KC has for sample/variant/adj filtering?
    if not args.skip_compute_qc_mt:
        logger.info("Filtering to bi-allelic, high-callrate, common SNPs for sample QC...")
        qc_mt = compute_qc_mt(get_ukbb_data(data_source, freeze, raw=True, split=False))
        qc_mt = qc_mt.naive_coalesce(5000)
        qc_mt.write(qc_mt_path(data_source, freeze), overwrite=args.overwrite)

        qc_mt = hl.read_matrix_table(qc_mt_path(data_source, freeze))
        qc_mt = filter_to_autosomes(qc_mt, reference_genome='GRCh38')
        qc_ht = hl.sample_qc(qc_mt).cols().select('sample_qc')
        qc_ht.write(qc_ht_path(data_source, freeze), overwrite=args.overwrite)

    if not args.skip_ld_prune_qc_mt:
        logger.info("LD prune qc MatrixTable sample QC...")
        qc_mt = hl.read_matrix_table(qc_mt_path(data_source, freeze))
        pruned_qc_mt = ld_prune_qc_mt(qc_mt)
        variants, samples = pruned_qc_mt.count()
        logger.info(f'{variants} variants found in LD-pruned MT')
        pruned_qc_mt = filter_to_adj(pruned_qc_mt)
        pruned_qc_mt.write(qc_mt_path(data_source, freeze, ld_pruned=True), overwrite=args.overwrite)

    if not args.skip_pc_relate:
        logger.info('Running PCA for PC-Relate...')
        pruned_qc_mt = hl.read_matrix_table(qc_mt_path(data_source, freeze, ld_pruned=True)).unfilter_entries()
        eig, scores, _ = hl.hwe_normalized_pca(pruned_qc_mt.GT, k=10, compute_loadings=False)
        scores.write(relatedness_pca_scores_ht_path(data_source, freeze), args.overwrite)

        logger.info('Running PC-Relate...')
        # NOTE: This needs SSDs on your workers (for the temp files) and no preemptible workers while the BlockMatrix writes
        relatedness_ht = hl.pc_relate(pruned_qc_mt.GT, min_individual_maf=0.05, scores_expr=scores[pruned_qc_mt.col_key].scores,
                                      block_size=4096, min_kinship=args.min_emission_kinship, statistics='all')
        relatedness_ht.write(relatedness_ht_path(data_source, freeze), args.overwrite)

    if not args.skip_filter_dups:
        logger.info("Filtering duplicate samples...")
        dups_ht = filter_dups(
            hl.read_table(relatedness_ht_path(data_source, freeze)),
            rank_dup_samples,
            [hl.read_table(qc_ht_path(data_source, freeze))]
        )
        dups_ht.write(duplicates_ht_path(data_source, freeze), overwrite=args.overwrite)

    if not args.skip_infer_families:
        logger.info("Inferring families")
        ped = get_ped(
            hl.read_table(relatedness_ht_path(data_source, freeze)),
            hl.read_table(duplicates_ht_path(data_source, freeze)),
            hl.read_table(sex_ht_path(data_source, freeze))
        )
        ped.write(inferred_ped_path(data_source, freeze))

    if args.skip_filter_related_samples:
        logger.info("Filtering related samples")
        related_samples_to_drop = get_related_samples_to_drop(
            hl.read_table(relatedness_ht_path(data_source, freeze)),
            args.min_filtering_kinship,
            rank_related_samples,
            [hl.read_table(qc_ht_path(data_source, freeze))]
        )
        related_samples_to_drop.write(related_drop_path(data_source, freeze), overwrite=args.overwrite)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--data_source', help='Source of the data, either broad or regeneron')
    parser.add_argument('--freeze', help='Data freeze to use', default=CURRENT_FREEZE)

    parser.add_argument('--skip_compute_qc_mt', help='Skop Compute matrix to be used in sample qc', action='store_true')
    parser.add_argument('--skip_ld_prune_qc_mt', help='Skip LD prunes the qc matrix', action='store_true')
    parser.add_argument('--skip_pc_relate', help='Skip runing PC-relate on all samples. NOTE: This needs SSDs on your workers (for the temp files) and no pre-emptibles while the BlockMatrix writes',
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
