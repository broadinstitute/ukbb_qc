from gnomad_hail import *
from gnomad_hail.utils.sample_qc import *
import gnomad_hail.resources.sample_qc as gres
from ukbb_qc.resources import *
import hail as hl
import argparse
import pickle


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("platform_pca")
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


def project_on_gnomad_pop_pcs(qc_mt: hl.MatrixTable) -> hl.Table:
    """
    Performs pc_project on a qc_mt using gnomAD population pca loadings and known pops
    :param MatrixTable qc_mt: raw matrix table to perform pc_project and pop assignment on
    :return: pc_project scores Table and
    :rtype: tuple(Table, Table)
    """
    qc_mt = qc_mt.select_entries('GT')

    # Load gnomAD metadata and population pc loadings
    gnomad_meta_exomes_ht = prep_meta(hl.read_table(metadata_exomes_ht_path(version=CURRENT_EXOME_META)))
    gnomad_meta_genomes_ht = prep_meta(hl.read_table(metadata_genomes_ht_path(version=CURRENT_GENOME_META)))
    gnomad_meta_ht = gnomad_meta_exomes_ht.union(gnomad_meta_genomes_ht)
    gnomad_loadings_ht = hl.read_table(gres.ancestry_pca_loadings_ht_path())
    gnomad_loadings_ht = liftover_ht(gnomad_loadings_ht)

    scores_ht = pc_project(qc_mt, gnomad_loadings_ht)
    scores_ht = scores_ht.annotate(pop_for_rf=hl.null(hl.tstr))
    scores_ht = scores_ht.select('pop_for_rf', **{f'PC{i + 1}': scores_ht.scores[i] for i in range(10)})

    joint_scores_ht = gnomad_meta_ht.union(scores_ht)

    return joint_scores_ht


def main(args):
    hl.init(log='/population_pca.log')

    if args.run_pc_project:
        qc_mt = get_ukbb_data(args.data_source, args.freeze, split=False, raw=True)
        joint_scores_ht = project_on_gnomad_pop_pcs(qc_mt)
        joint_scores_ht.checkpoint(ancestry_pc_project_scores_ht_path(args.data_source, args.freeze, "joint"), overwrite=args.overwrite)
        joint_scores_pd = joint_scores_ht.to_pandas()

        joint_pops_pd, joint_pops_rf_model = assign_population_pcs(
            joint_scores_pd,
            pc_cols=[f'PC{i + 1}' for i in range(args.n_pcs)],
            known_col='pop_for_rf',
            min_prob=args.min_pop_prob
        )

        joint_pops_ht = hl.Table.from_pandas(joint_pops_pd, key=list(joint_scores_ht.key))
        scores_ht = joint_scores_ht.filter(hl.is_missing(joint_scores_ht.pop_for_rf))
        scores_ht = scores_ht.annotate(pop=joint_pops_ht[scores_ht.key])
        scores_ht.checkpoint(ancestry_pc_project_scores_ht_path(args.data_source, args.freeze), overwrite=args.overwrite)

        logger.info(
            "Found the following sample count after population assignment (reduced to only input samples): {}".format(
                scores_ht.aggregate(hl.agg.counter(scores_ht.pop['pop']))))

        with hl.hadoop_open(qc_temp_data_prefix(args.data_source, args.freeze) + 'project_gnomad_pop_rf_model.pkl', 'wb') as out:
            pickle.dump(joint_pops_rf_model, out)


    # TODO: Need to add in the code for doing PCA on the full dataset when we have gnomAD + UKBB callset


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE)

    parser.add_argument('--run_pc_project', help='Runs pc project and population assignment using gnomAD data', action='store_true')
    parser.add_argument('--n_pcs', help='Number of PCs to compute (default: 10)', default=10, type=int)
    parser.add_argument('--min_pop_prob',
                     help='Minimum probability of belonging to a given population for assignment (if below, the sample is labeled as "oth" (default: 0.9)',
                     default=0.9,
                     type=float)

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)