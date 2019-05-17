from gnomad_hail import *
from gnomad_hail.utils.sample_qc import *
import gnomad_hail.resources.sample_qc as gres
from ukbb_qc.resources import *
import hail as hl
import argparse
import pickle
import numpy as np
import hdbscan


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("population_pca")
logger.setLevel(logging.INFO)


# The following are from Mike with some changes made
def prep_meta(ht: hl.Table) -> hl.Table:
    """
    Preps gnomAD exome and genome metadata for population PC
    :param Table ht: Hail table from gnomAD's Column annotations
    :return: Hail Table ready for joining
    :rtype: Table
    """
    ht = ht.key_by('s')
    ht = ht.rename({'qc_pop':'pop_for_rf'})
    ht = ht.select('pop_for_rf', **{f'PC{i}': ht[f'PC{i}'] for i in range(1, 11)})
    ht = ht.filter(hl.is_defined(ht.PC1))
    return ht


def liftover_ht(ht: hl.Table, ref_from: str = 'GRCh37', ref_to: str = 'GRCh38') -> hl.Table:
    """
    Liftover input ht from one build to another
    :param Table ht: Table to be lifted over
    :param str ref_from: Current reference of ht
    :param str ref_to: Reference to liftover to
    :return: lifted over Table
    :rtype: Table
    """
    logger.info(f'Lifting a {ref_from} genome over to {ref_to}')
    rgfrom = hl.get_reference(ref_from)
    rgto = hl.get_reference(ref_to)
    rgfrom.add_liftover(f'gs://hail-common/references/{ref_from.lower()}_to_{ref_to.lower()}.over.chain.gz', rgto)
    ht = ht.annotate(new_locus=hl.liftover(ht.locus, ref_to, include_strand=True),
                     old_locus=ht.locus)
    ht = ht.filter(hl.is_defined(ht.new_locus) & ~ht.new_locus.is_negative_strand)
    ht = ht.key_by(locus=ht.new_locus.result, alleles=ht.alleles)

    return ht


def project_on_gnomad_pop_pcs(mt: hl.MatrixTable) -> hl.Table:
    """
    Performs pc_project on a mt using gnomAD population pca loadings and known pops
    :param MatrixTable mt: raw matrix table to perform pc_project and pop assignment on
    :return: pc_project scores Table and
    :rtype: tuple(Table, Table)
    """
    mt = mt.select_entries('GT')

    # Load gnomAD metadata and population pc loadings
    gnomad_meta_exomes_ht = prep_meta(hl.read_table(metadata_exomes_ht_path(version=CURRENT_EXOME_META)))
    gnomad_meta_genomes_ht = prep_meta(hl.read_table(metadata_genomes_ht_path(version=CURRENT_GENOME_META)))
    gnomad_meta_ht = gnomad_meta_exomes_ht.union(gnomad_meta_genomes_ht)
    gnomad_loadings_ht = hl.read_table(gres.ancestry_pca_loadings_ht_path())
    gnomad_loadings_ht = liftover_ht(gnomad_loadings_ht)

    scores_ht = pc_project(mt, gnomad_loadings_ht)
    scores_ht = scores_ht.annotate(pop_for_rf=hl.null(hl.tstr))
    scores_ht = scores_ht.select('pop_for_rf', **{f'PC{i + 1}': scores_ht.scores[i] for i in range(10)})

    joint_scores_ht = gnomad_meta_ht.union(scores_ht)

    return joint_scores_ht


def assign_cluster_pcs(pca_ht: hl.Table, pc_scores_ann: str = 'scores', hdbscan_min_cluster_size: int = 50, hdbscan_min_samples: int = None) -> hl.Table:
    """

    :param MatrixTable pca_ht:
    :param str pc_scores_ann:
    :param int hdbscan_min_cluster_size:
    :param int hdbscan_min_samples:
    :return: cluster assignment Table
    :rtype: Table
    """
    pca_pd = pca_ht.to_pandas()
    pc_scores = np.matrix(pca_pd[pc_scores_ann].tolist())
    logger.info('Assigning clusters to {} exome samples in MT...'.format(len(pc_scores)))

    clusterer = hdbscan.HDBSCAN(min_cluster_size=hdbscan_min_cluster_size, min_samples=hdbscan_min_samples)
    cluster_labels = clusterer.fit_predict(pc_scores)
    n_clusters = len(set(cluster_labels)) - (-1 in cluster_labels)  # NOTE: -1 is the label for noisy (un-classifiable) data points
    logger.info('Found {} unique clusters...'.format(n_clusters))

    pca_pd['ancestry_cluster'] = cluster_labels
    ht = hl.Table.from_pandas(pca_pd, key=[*pca_ht.key])
    ht = ht.annotate(ancestry_cluster=hl.int32(ht.ancestry_cluster))
    return ht


def main(args):
    hl.init(log='/population_pca.log')

    data_source = args.data_source
    freeze = args.freeze

    if args.run_pca:
        logger.info("Running population PCA")
        pca_evals, pop_pca_scores_ht, pop_pca_loadings_ht = run_pca_with_relateds(
            hl.read_matrix_table(qc_mt_path(data_source, freeze, ld_pruned=True)),
            hl.read_table(related_drop_path(data_source, freeze)),
            args.n_pcs
        )
        pop_pca_scores_ht = pop_pca_scores_ht.annotate_globals(n_pcs=args.n_pcs)
        pop_pca_loadings_ht.write(ancestry_pca_loadings_ht_path(data_source, freeze), args.overwrite)
        pop_pca_scores_ht.write(ancestry_pca_scores_ht_path(data_source, freeze), args.overwrite)

    if args.assign_clusters:
        logger.info("Assigning PCA clustering")
        scores_ht = hl.read_table(ancestry_pca_scores_ht_path(data_source, freeze))
        pops_ht = assign_cluster_pcs(scores_ht,
                                     hdbscan_min_cluster_size=args.hdbscan_min_cluster_size,
                                     hdbscan_min_samples=args.hdbscan_min_samples)
        pops_ht.annotate_globals(hdbscan_min_cluster_size=args.hdbscan_min_cluster_size)
        pops_ht.write(ancestry_cluster_ht_path(data_source, freeze), args.overwrite)

    if args.run_pc_project:
        # Note: I used all workers for this as it kept failing with preemptibles
        mt = get_ukbb_data(data_source, freeze, split=False, adj=True)
        joint_scores_ht = project_on_gnomad_pop_pcs(mt)
        joint_scores_ht.write(ancestry_pc_project_scores_ht_path(data_source, freeze, "joint"), overwrite=args.overwrite)

    if args.run_rf:
        joint_scores_ht = hl.read_table(ancestry_pc_project_scores_ht_path(data_source, freeze, "joint"))
        joint_scores_pd = joint_scores_ht.to_pandas()
        joint_pops_pd, joint_pops_rf_model = assign_population_pcs(
            joint_scores_pd,
            pc_cols=[f'PC{i + 1}' for i in range(args.n_pcs)],
            known_col='pop_for_rf',
            min_prob=args.min_pop_prob
        )

        joint_pops_ht = hl.Table.from_pandas(joint_pops_pd, key=list(joint_scores_ht.key))
        scores_ht = joint_scores_ht.filter(hl.is_missing(joint_scores_ht.pop_for_rf))
        scores_ht = scores_ht.drop('pop_for_rf')
        joint_pops_ht = joint_pops_ht.drop('pop_for_rf')
        scores_ht = scores_ht.annotate(pop=joint_pops_ht[scores_ht.key])
        scores_ht = scores_ht.annotate_globals(n_pcs=args.n_pcs,
                                               min_prob=args.min_pop_prob)
        scores_ht = scores_ht.checkpoint(ancestry_pc_project_scores_ht_path(data_source, freeze), overwrite=args.overwrite)

        logger.info(
            "Found the following sample count after population assignment (reduced to only input samples): {}".format(
                scores_ht.aggregate(hl.agg.counter(scores_ht.pop['pop']))))

        with hl.hadoop_open(qc_temp_data_prefix(data_source, freeze) + 'project_gnomad_pop_rf_model.pkl', 'wb') as out:
            pickle.dump(joint_pops_rf_model, out)

    if args.assign_hybrid_ancestry:
        pc_project_ht = hl.read_table(ancestry_pc_project_scores_ht_path(data_source, freeze))
        pc_cluster_ht = hl.read_table(ancestry_cluster_ht_path(data_source, freeze))
        pc_scores_ht = hl.read_table(ancestry_pca_scores_ht_path(data_source, freeze))

        pop_ht = pc_project_ht.select(gnomad_pc_project_pop=pc_project_ht.pop.pop,
                                      **{f'gnomad_pc_project_PC{i+1}': pc_project_ht[f'PC{i+1}']
                                         for i in range(hl.eval(pc_project_ht.n_pcs))
                                         }
                                      )

        pop_ht = pop_ht.annotate(HDBSCAN_pop_cluster=pc_cluster_ht[pop_ht.key].ancestry_cluster)
        pop_ht = pop_ht.annotate(hybrid_pop=hl.case()
                                 .when(pop_ht.HDBSCAN_pop_cluster == -1, pop_ht.gnomad_pc_project_pop)
                                 .default(hl.str(pop_ht.HDBSCAN_pop_cluster)),
                                 **{f'pop_pca_PC{i+1}': pc_scores_ht[pop_ht.key].scores[i]
                                    for i in range(hl.eval(pc_scores_ht.n_pcs))
                                    }
                                 )
        pop_ht.write(ancestry_hybrid_ht_path(data_source, freeze), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE)

    parser.add_argument('--run_pc_project', help='Runs pc project and population assignment using gnomAD data', action='store_true')
    parser.add_argument('--n_pcs', help='Number of PCs to compute (default: 10)', default=10, type=int)
    parser.add_argument('--run_rf', 
                        help='Create random forest model to assign population labels based on PCA results', action='store_true')
    parser.add_argument('--min_pop_prob',
                     help='Minimum probability of belonging to a given population for assignment (if below, the sample is labeled as "oth" (default: 0.9)',
                     default=0.5,
                     type=float)

    parser.add_argument('--run_pca', help='Runs pop PCA on pruned qc MT', action='store_true')
    parser.add_argument('--assign_clusters', help='Assigns clusters based on PCA results using HDBSCAN', action='store_true')
    parser.add_argument('--hdbscan_min_samples', help='Minimum samples parameter for HDBSCAN. If not specified, --hdbscan_min_cluster_size is used.', type=int, required=False)
    parser.add_argument('--hdbscan_min_cluster_size', help='Minimum cluster size parameter for HDBSCAN.', type=int, default=50)
    parser.add_argument('--assign_hybrid_ancestry', help='Assigns samples to HDBSCAN clusters where available otherwise uses pc_project assignments.',
                        action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
