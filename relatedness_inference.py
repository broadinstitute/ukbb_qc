from resources import *
from gnomad_hail import *
from gnomad_hail.utils.sample_qc import *
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
    #hl.init(log='/relatedness.log', tmp_dir='hdfs:///pc_relate.tmp/')

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
        sample_qc_ht = hl.read_table(qc_ht_path(data_source, freeze))
        samples_rankings_ht = sample_qc_ht.select(rank=-1 * sample_qc_ht.sample_qc.dp_stats.mean)
        dups_ht = filter_duplicate_samples(
            hl.read_table(relatedness_ht_path(data_source, freeze)),
            samples_rankings_ht
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
        logger.info(f'{related_samples_to_drop_ht.count()} samples dropped from callset using maximal independent set')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE)

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
