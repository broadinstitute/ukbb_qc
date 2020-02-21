import argparse
import hail as hl
import hdbscan
import logging
import numpy as np
import pickle
from gnomad_hail.utils.liftover import annotate_snp_mismatch, get_liftover_genome, lift_data
from gnomad_hail.utils.sample_qc import run_pca_with_relateds
from gnomad_hail.utils.slack import try_slack
from gnomad_qc.v2.resources.sample_qc import (
                                    ancestry_cluster_ht_path, ancestry_cluster_array_ht_path,
                                    ancestry_cluster_joint_scratch_array_ht_path,
                                    ancestry_pc_project_scores_ht_path,
                                    ancestry_pca_loadings_ht_path, ancestry_pca_scores_ht_path,
                                    ancestry_hybrid_ht_path,
                                    gnomad_ancestry_loadings_liftover_path,
                                    qc_mt_path, qc_temp_data_prefix
                                    )
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import related_drop_path
from ukbb_qc.utils.utils import remove_hard_filter_samples


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("population_pca")
logger.setLevel(logging.INFO)


# TODO: The two following functions are from Mike with some changes made - will change to use Mike"s code when it is committed
def prep_meta(ht: hl.Table) -> hl.Table:
    """
    Preps gnomAD exome and genome metadata for population PC
    :param Table ht: Hail table from gnomAD"s Column annotations
    :return: Hail Table ready for joining
    :rtype: Table
    """
    ht = ht.key_by("s")
    ht = ht.rename({"qc_pop":"pop_for_rf"})
    ht = ht.select("pop_for_rf", **{f"PC{i}": ht[f"PC{i}"] for i in range(1, 11)})
    ht = ht.filter(hl.is_defined(ht.PC1))
    return ht


def project_on_gnomad_pop_pcs(mt: hl.MatrixTable) -> hl.Table:
    """
    Performs pc_project on a mt using gnomAD population pca loadings and known pops
    :param MatrixTable mt: raw matrix table to perform pc_project and pop assignment on
    :return: pc_project scores Table and
    :rtype: tuple(Table, Table)
    """

    # Load gnomAD metadata and population pc loadings
    gnomad_meta_exomes_ht = prep_meta(hl.read_table(metadata_exomes_ht_path(version=CURRENT_EXOME_META)))
    gnomad_meta_genomes_ht = prep_meta(hl.read_table(metadata_genomes_ht_path(version=CURRENT_GENOME_META)))
    gnomad_meta_ht = gnomad_meta_exomes_ht.union(gnomad_meta_genomes_ht)
    gnomad_loadings_ht = hl.read_table(gnomad_ancestry_loadings_liftover_path())
    gnomad_loadings_ht = gnomad_loadings_ht.filter(
        ~gnomad_loadings_ht.reference_mismatch & ~gnomad_loadings_ht.new_locus.is_negative_strand)

    scores_ht = pc_project(mt, gnomad_loadings_ht)
    scores_ht = scores_ht.annotate(pop_for_rf=hl.null(hl.tstr))
    scores_ht = scores_ht.select("pop_for_rf", **{f"PC{i + 1}": scores_ht.scores[i] for i in range(10)})

    joint_scores_ht = gnomad_meta_ht.union(scores_ht)

    return joint_scores_ht


def assign_cluster_pcs(pca_ht: hl.Table, pc_scores_ann: str = "scores", hdbscan_min_cluster_size: int = 50,
                       hdbscan_min_samples: int = None) -> hl.Table:
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
    logger.info("Assigning clusters to {} exome samples in MT...".format(len(pc_scores)))

    clusterer = hdbscan.HDBSCAN(min_cluster_size=hdbscan_min_cluster_size, min_samples=hdbscan_min_samples)
    cluster_labels = clusterer.fit_predict(pc_scores)
    n_clusters = len(set(cluster_labels)) - (
                -1 in cluster_labels)  # NOTE: -1 is the label for noisy (un-classifiable) data points
    logger.info("Found {} unique clusters...".format(n_clusters))

    pca_pd["ancestry_cluster"] = cluster_labels
    ht = hl.Table.from_pandas(pca_pd, key=[*pca_ht.key])
    ht = ht.annotate(ancestry_cluster=hl.int32(ht.ancestry_cluster))
    return ht


def load_ukbb_array_pcs(data_source, freeze, num_pcs=16):
    sample_map_ht = hl.read_table(array_sample_map_ht(data_source, freeze))
    array_pcs = hl.import_table("gs://broad-ukbb/resources/array/ukb_sqc_v2_ukb26041.txt", impute=True, delimiter="\t",
                                no_header=True)

    sample_map_ht = sample_map_ht.key_by("ukbb_app_26041_id")
    array_pcs = array_pcs.key_by(s=sample_map_ht[hl.str(array_pcs.f1)].s)
    array_pcs = array_pcs.filter(hl.is_defined(array_pcs.s))
    array_pcs = array_pcs.annotate(scores=hl.array([array_pcs[f"f{i+26}"] for i in range(0, num_pcs)]))

    return array_pcs


def main(args):
    hl.init(log="/population_pca.log", default_reference="GRCh38")

    data_source = args.data_source
    freeze = args.freeze
    n_pcs = args.n_pcs

    if args.liftover_gnomad_ancestry_loadings:
        gnomad_loadings_ht = hl.read_table(ancestry_pca_loadings_ht_path())

        logger.info("Preparing reference genomes for liftover")
        source, target = get_liftover_genome(gnomad_loadings_ht)

        logger.info(f"Lifting data to {target.name}")
        gnomad_loadings_ht = lift_data(
            gnomad_loadings_ht, gnomad=False, data_type=None,
            path=gnomad_ancestry_loadings_liftover_path(checkpoint=True), rg=target,
            overwrite=args.overwrite
        )

        logger.info("Checking SNPs for reference mismatches")
        gnomad_loadings_ht = annotate_snp_mismatch(gnomad_loadings_ht, data_type=None, rg=target)
        gnomad_loadings_ht.write(gnomad_ancestry_loadings_liftover_path())

    if args.run_pca:
        logger.info("Running population PCA")
        qc_mt = remove_hard_filter_samples(
            data_source, freeze,
            hl.read_matrix_table(qc_mt_path(data_source, freeze, ld_pruned=True))
        )
        related_ht = related_drop_path(data_source, freeze)
        related_ht = related_ht.filter(related_ht.relationship == "2nd degree relatives")
        pca_evals, pop_pca_scores_ht, pop_pca_loadings_ht = run_pca_with_relateds(
            qc_mt,
            related_ht,
            n_pcs
        )
        pop_pca_scores_ht = pop_pca_scores_ht.annotate_globals(n_pcs=n_pcs)
        pop_pca_loadings_ht.write(ancestry_pca_loadings_ht_path(data_source, freeze), args.overwrite)
        pop_pca_scores_ht.write(ancestry_pca_scores_ht_path(data_source, freeze), args.overwrite)

    if args.assign_clusters:
        logger.info("Assigning PCA clustering")
        scores_ht = hl.read_table(ancestry_pca_scores_ht_path(data_source, freeze))
        pops_ht = assign_cluster_pcs(
            scores_ht,
            hdbscan_min_cluster_size=args.hdbscan_min_cluster_size,
            hdbscan_min_samples=args.hdbscan_min_samples
        )
        pops_ht.annotate_globals(hdbscan_min_cluster_size=args.hdbscan_min_cluster_size)
        pops_ht.write(ancestry_cluster_ht_path(data_source, freeze), args.overwrite)

    if args.assign_clusters_array_pcs:
        logger.info("Load UKBB array PC data")
        array_pc_ht = load_ukbb_array_pcs(datasource, freeze)
        logger.info("Assigning PCA clustering")
        pops_ht = assign_cluster_pcs(
            array_pc_ht,
            hdbscan_min_cluster_size=args.hdbscan_min_cluster_size,
            hdbscan_min_samples=args.hdbscan_min_samples
        )
        pops_ht.annotate_globals(hdbscan_min_cluster_size=args.hdbscan_min_cluster_size)
        pops_ht.write(ancestry_cluster_array_ht_path(data_source, freeze), args.overwrite)

    if args.assign_clusters_joint_scratch_array_pcs:
        logger.info("Load exome scratch PC data")
        scores_ht = hl.read_table(ancestry_pca_scores_ht_path(data_source, freeze))

        logger.info("Load UKBB array PC data")
        array_pc_ht = load_ukbb_array_pcs(data_source, freeze)
        array_pc_ht = array_pc_ht.annotate(exome_scores=scores_ht[array_pc_ht.key].scores)
        array_pc_ht = array_pc_ht.filter(hl.is_defined(array_pc_ht.exome_scores))
        array_pc_ht = array_pc_ht.select(
            **{f"PC{i + 1}": array_pc_ht.exome_scores[i] for i in range(n_pcs)},
            **{f"PC{i + n_pcs + 1}": array_pc_ht.scores[i] for i in range(args.n_array_pcs)}
        )
        array_pc_ht = array_pc_ht.annotate(
            scores=hl.array([array_pc_ht[f"PC{i + 1}"] for i in range(n_pcs + args.n_array_pcs)])
        )

        logger.info("Assigning PCA clustering")
        pops_ht = assign_cluster_pcs(
            array_pc_ht,
            hdbscan_min_cluster_size=args.hdbscan_min_cluster_size,
            hdbscan_min_samples=args.hdbscan_min_samples
        )
        pops_ht.annotate_globals(hdbscan_min_cluster_size=args.hdbscan_min_cluster_size)
        pops_ht.write(ancestry_cluster_joint_scratch_array_ht_path(data_source, freeze), args.overwrite)

    if args.run_pc_project:
        # NOTE: I used all workers for this as it kept failing with preemptibles
        # Using the split mt led to better clustering and fewer others in Broad freeze 4
        mt = get_ukbb_data(data_source, freeze, split=True, adj=True)
        mt = remove_hard_filter_samples(data_source, freeze, mt)
        joint_scores_ht = project_on_gnomad_pop_pcs(mt)
        joint_scores_ht.write(
            ancestry_pc_project_scores_ht_path(data_source, freeze, "joint"), 
            overwrite=args.overwrite
        )

    if args.run_rf:
        joint_scores_ht = hl.read_table(ancestry_pc_project_scores_ht_path(data_source, freeze, "joint"))
        joint_scores_pd = joint_scores_ht.to_pandas()
        joint_pops_pd, joint_pops_rf_model = assign_population_pcs(
            joint_scores_pd,
            pc_cols=[f"PC{i + 1}" for i in range(n_pcs)],
            known_col="pop_for_rf",
            min_prob=args.min_pop_prob
        )

        joint_pops_ht = hl.Table.from_pandas(joint_pops_pd, key=list(joint_scores_ht.key))
        scores_ht = joint_scores_ht.filter(hl.is_missing(joint_scores_ht.pop_for_rf))
        scores_ht = scores_ht.drop("pop_for_rf")
        joint_pops_ht = joint_pops_ht.drop("pop_for_rf")
        scores_ht = scores_ht.annotate(pop=joint_pops_ht[scores_ht.key])
        scores_ht = scores_ht.annotate_globals(
            n_pcs=n_pcs,
            min_prob=args.min_pop_prob
        )
        scores_ht = scores_ht.checkpoint(
            ancestry_pc_project_scores_ht_path(data_source, freeze), 
            overwrite=args.overwrite
        )

        logger.info(
            "Found the following sample count after population assignment (reduced to only input samples): {}".format(
                scores_ht.aggregate(hl.agg.counter(scores_ht.pop["pop"])))
        )

        with hl.hadoop_open(
            qc_temp_data_prefix(data_source, freeze) + "project_gnomad_pop_rf_model.pkl", "wb"
        ) as out:
            pickle.dump(joint_pops_rf_model, out)

    if args.assign_hybrid_ancestry:
        pc_project_ht = hl.read_table(ancestry_pc_project_scores_ht_path(data_source, freeze))
        pc_cluster_ht = hl.read_table(ancestry_cluster_ht_path(data_source, freeze))
        pc_scores_ht = hl.read_table(ancestry_pca_scores_ht_path(data_source, freeze))

        pop_ht = pc_project_ht.select(
            gnomad_pc_project_pop=pc_project_ht.pop.pop,
            **{
                f"gnomad_pc_project_PC{i+1}": pc_project_ht[f"PC{i+1}"]
                for i in range(hl.eval(pc_project_ht.n_pcs))
            }
        )

        pop_ht = pop_ht.annotate(HDBSCAN_pop_cluster=pc_cluster_ht[pop_ht.key].ancestry_cluster)
        pop_ht = pop_ht.annotate(
            hybrid_pop=hl.case()
                .when(pop_ht.HDBSCAN_pop_cluster == -1, pop_ht.gnomad_pc_project_pop)
                .default(hl.str(pop_ht.HDBSCAN_pop_cluster)),
                **{
                    f"pop_pca_PC{i+1}": pc_scores_ht[pop_ht.key].scores[i]
                    for i in range(hl.eval(pc_scores_ht.n_pcs))
                }
        )
        pop_ht.write(ancestry_hybrid_ht_path(data_source, freeze), args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-s", "--data_source", help="Data source", choices=["regeneron", "broad"], default="broad"
    )
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )

    parser.add_argument(
        "--liftover_gnomad_ancestry_loadings", 
        help="Performs liftover on the gnomad joint ancestry loadings", action="store_true"
    )
    parser.add_argument(
        "--run_pca", help="Runs pop PCA on pruned qc MT", action="store_true"
    )
    parser.add_argument(
        "--assign_clusters", help="Assigns clusters based on PCA results using HDBSCAN", action="store_true"
    )
    parser.add_argument(
        "--assign_clusters_array_pcs", 
        help="Assigns clusters based on PCA results using HDBSCAN", action="store_true"
    )
    parser.add_argument(
        "--n_array_pcs", help="Number of array PCs to use in clustering (default: 16)", 
        default=16, type=int
    )
    parser.add_argument(
        "--assign_clusters_joint_scratch_array_pcs", 
        help="Assigns clusters based on PCA results using HDBSCAN", action="store_true"
    )
    parser.add_argument(
        "--hdbscan_min_samples", 
        help="Minimum samples parameter for HDBSCAN. If not specified, --hdbscan_min_cluster_size is used.", 
        type=int, required=False
    )
    parser.add_argument(
        "--hdbscan_min_cluster_size", help="Minimum cluster size parameter for HDBSCAN.", 
        type=int, default=50
    )

    parser.add_argument(
        "--run_pc_project", 
        help="Runs pc project and population assignment using gnomAD data", action="store_true"
    )
    parser.add_argument(
        "--n_pcs",
        help="Number of PCs to compute or use in population inference (default: 10)",
        default=10,
        type=int
    )
    parser.add_argument(
        "--run_rf", 
        help="Create random forest model to assign population labels based on PCA results", 
        action="store_true"
    )
    parser.add_argument(
        "--min_pop_prob",
        help="Minimum probability of belonging to a given population for assignment (if below, the sample is labeled as "oth" (default: 0.5)",
        default=0.5,
        type=float
    )

    parser.add_argument(
        "--assign_hybrid_ancestry", 
        help="Assigns samples to HDBSCAN clusters where available otherwise uses pc_project assignments.",
        action="store_true"
    )

    parser.add_argument("-o", "--overwrite", help="Overwrite all data from this subset (default: False)", action="store_true")
    parser.add_argument("--slack_channel", help="Slack channel to post results and notifications to.")
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
