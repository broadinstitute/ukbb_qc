import sys
sys.path.append('/home/hail')
print(sys.path)

from gnomad_hail import *
from ukbb_qc.resources import *
import hail as hl
import argparse



logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("relatedness")
logger.setLevel(logging.INFO)

mt = get_gnomad_data(data_type, raw=True,
                     split=False)  # NOTE: using full calls since hardcalls doesn't exist at this stage


mt.annotate_cols(callrate=hl.agg.fraction(hl.is_defined(mt.GT))).naive_coalesce(5000).write(qc_mt_path(data_type),
                                                                                            overwrite=args.overwrite)
qc_mt = hl.read_matrix_table(qc_mt_path(data_type))

logger.info("Importing metadata...")
meta_ht = hl.import_table(qc_meta_path(data_type), impute=True, types={'age': hl.tfloat64}).key_by('s')
qc_mt = qc_mt.annotate_cols(**meta_ht[qc_mt.s])




def get_qc_mt(mt: hl.MatrixTable, min_af: float = 0.001, min_callrate: float = 0.99, ld_r2: float = 0.1) -> hl.MatrixTable:
    logger.info("Filtering to bi-allelic, high-callrate, common SNPs for sample QC...")
    mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]) &
                        (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > min_af) &
                        (hl.agg.fraction(hl.is_defined(mt.GT)) > min_callrate))

    pruned_ht = hl.ld_prune(qc_mt.GT, r2=ld_r2)
    qc_mt = qc_mt.filter_rows(hl.is_defined(pruned_ht[qc_mt.row_key]))
    qc_mt = qc_mt.annotate_globals(
        qc_mt_params=hl.struct(
            min_af=min_af,
            min_callrate=min_callrate,
            ld_r2=ld_r2
        )
    )
    qc_mt.annotate_cols(sample_callrate=hl.agg.fraction(hl.is_defined(qc_mt.GT)))
    qc_mt = filter_to_adj(qc_mt)
    return qc_mt



def main(args):
    hl.init(log='/relatedness.log')

    if args.compute_qc_mt:
        logger.info("Creating QC MatrixTable")
        qc_mt = get_qc_mt(hl.read_matrix_table(args.input_mt))
        qc_mt = qc_mt.repartition(n_partitions=200)
        qc_mt.write(f'{output_prefix}.qc.mt', overwrite=args.overwrite)

    if args.run_pc_relate:
        logger.info('Running PCA for PC-Relate')
        qc_mt = hl.read_matrix_table(f'{output_prefix}.qc.mt')._unfilter_entries()
        if args.test_run:
            qc_mt = hl.filter_intervals(qc_mt, [hl.parse_locus_interval('chr20', reference_genome='GRCh38')])
        eig, scores, _ = hl.hwe_normalized_pca(qc_mt.GT, k=10, compute_loadings=False)
        scores.write(f'{output_prefix}.pruned.pca_scores.ht', args.overwrite)

        logger.info('Running PC-Relate')
        scores = hl.read_table(f'{output_prefix}.pruned.pca_scores.ht')
        # NOTE: This needs SSDs on your workers (for the temp files) and no preemptible workers while the BlockMatrix writes
        relatedness_ht = hl.pc_relate(qc_mt.GT, min_individual_maf=0.05, scores_expr=scores[qc_mt.col_key].scores,
                                      block_size=4096, min_kinship=args.min_emission_kinship, statistics='all')
        relatedness_ht.write(relatedness_ht_path(), args.overwrite)

    if args.filter_dups:
        logger.info("Filtering duplicate samples")
        dups_ht = filter_dups(
            hl.read_table(f'{output_prefix}.relatedness.ht'),
            rank_dup_samples,
            [hl.read_table(f'{output_prefix}.sample_qc.ht')]
        )
        dups_ht.write(f'{output_prefix}.dups.ht', overwrite=args.overwrite)

    if args.infer_families:
        logger.info("Inferring families")
        ped = get_ped(
            hl.read_table(f'{output_prefix}.relatedness.ht'),
            hl.read_table(f'{output_prefix}.dups.ht'),
            hl.read_table(f'{output_prefix}.sex.ht')
        )
        ped.write(f"{output_prefix}.ped")

    if args.filter_related_samples:
        logger.info("Filtering related samples")
        related_samples_to_drop = get_related_samples_to_drop(
            hl.read_table(f'{output_prefix}.relatedness.ht'),
            args.min_filtering_kinship,
            rank_related_samples,
            [
                hl.read_table(args.meta),
                hl.read_table(f'{output_prefix}.sample_qc.ht'),
                hl.import_fam(f"{output_prefix}.ped", delimiter="\t")
            ]
        )
        related_samples_to_drop.write(f'{output_prefix}.related_samples_to_drop.ht', overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)',
                        action='store_true')
    parser.add_argument('--test_run', 'Test with chr20',action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--data_source', help='Source of the data, either broad or regeneron')
    parser.add_argument('--freeze', help='Data freeze to use', default=CURRENT_FREEZE)
    parser.add_argument('--output_dir', help='Directory to output PCA files to.')

    parser.add_argument('--run_pc_relate', help='Runs PC-relate on all samples. NOTE: This needs SSDs on your workers (for the temp files) and no pre-emptibles while the BlockMatrix writes', action='store_true')
    parser.add_argument('--min_emission_kinship', help='Minimum kinship threshold for emitting a pair of samples in PC relate and filtering related individuals.', default=0.05, type=float)
    parser.add_argument('--min_filtering_kinship', help='Minimum kinship threshold for filtering a pair of samples in PC relate and filtering related individuals. (Default = 0.08838835; 2nd degree relatives)', default=0.08838835, type=float)
    parser.add_argument('--filter_dups', help='Filter duplicated samples', action='store_true')
    parser.add_argument('--infer_families', help='Extracts duplicate samples and infers families samples based on PC-relate results', action='store_true')
    parser.add_argument('--filter_related_samples', help='Filter related samples (based on the pairs present from the --run_pc_relate and using the --min_filtering_kinship value for that run)', action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)