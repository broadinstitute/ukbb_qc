import argparse
import logging
import pickle

import hail as hl
import hdbscan
import numpy as np

from gnomad.sample_qc.ancestry import (
    assign_population_pcs,
    pc_project,
    run_pca_with_relateds,
)
from gnomad.utils.liftover import (
    annotate_snp_mismatch,
    get_liftover_genome,
    lift_data,
)
from gnomad.utils.slack import try_slack
from gnomad.utils.sparse_mt import densify_sites
from gnomad_qc.v2.resources.basics import get_gnomad_meta
from gnomad_qc.v2.resources.sample_qc import (
    ancestry_pca_loadings_ht_path as gnomad_ancestry_pca_loadings_ht_path,
)
from ukbb_qc.resources.basics import (
    array_sample_map_ht_path,
    get_ukbb_data,
    last_END_positions_ht_path,
    logging_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import (
    ancestry_cluster_ht_path,
    ancestry_hybrid_ht_path,
    ancestry_pc_project_scores_ht_path,
    ancestry_pca_loadings_ht_path,
    ancestry_pca_scores_ht_path,
    get_ukbb_array_pcs_ht_path,
    get_ukbb_self_reported_ancestry_path,
    gnomad_ancestry_loadings_liftover_path,
    qc_mt_path,
    qc_temp_data_prefix,
    related_drop_path,
)
from ukbb_qc.utils.utils import remove_hard_filter_samples


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("population_pca")
logger.setLevel(logging.INFO)


def project_on_gnomad_pop_pcs(
    mt: hl.MatrixTable, last_END_ht: hl.Table, n_pcs: int = 10
) -> hl.Table:
    """
    Performs pc_project on a mt using gnomAD population pca loadings and known pops

    :param MatrixTable mt: Raw MatrixTable to perform `pc_project` and pop assignment on
    :param Table last_END_ht: Table with most upstream informative position (last END) annotated for each locus.
        Used when densifying data.
    :return: pc_project scores Table
    :param int n_pcs: Number of PCs to keep from gnomAD for `pc_project`. Default is 10.
    :rtype: Table
    """
    # Load gnomAD metadata and population pc loadings
    gnomad_meta_ht = get_gnomad_meta("joint", full_meta=True)
    gnomad_meta_ht = gnomad_meta_ht.key_by("s")
    # NOTE: Filtering to defined qc_pop to only keep gnomAD samples that have inferred ancestries
    # This is to avoid samples that were hard filtered prior to population inference
    gnomad_meta_ht = gnomad_meta_ht.filter(hl.is_defined(gnomad_meta_ht.qc_pop))
    gnomad_meta_ht = gnomad_meta_ht.rename({"qc_pop": "pop_for_rf"})
    gnomad_meta_ht = gnomad_meta_ht.select(
        "pop_for_rf", scores=[gnomad_meta_ht[f"PC{i}"] for i in range(1, n_pcs + 1)]
    )

    gnomad_loadings_ht = hl.read_table(gnomad_ancestry_loadings_liftover_path())
    gnomad_loadings_ht = gnomad_loadings_ht.filter(
        ~gnomad_loadings_ht.reference_mismatch
        & ~gnomad_loadings_ht.new_locus.is_negative_strand
    )

    # Need to densify before running pc project
    mt = densify_sites(mt, gnomad_loadings_ht, last_END_ht)

    scores_ht = pc_project(mt, gnomad_loadings_ht)
    scores_ht = scores_ht.annotate(pop_for_rf=hl.null(hl.tstr))
    scores_ht = scores_ht.select("pop_for_rf", scores=scores_ht.scores[:n_pcs])

    joint_scores_ht = gnomad_meta_ht.union(scores_ht)

    return joint_scores_ht


# TODO: Think about modifying this function and assign_platform_from_pcs in gnomad.utils.sample_qc to be more general
def assign_cluster_from_pcs(
    pca_ht: hl.Table,
    pc_scores_ann: str = "scores",
    hdbscan_min_cluster_size: int = 50,
    hdbscan_min_samples: int = None,
) -> hl.Table:
    """
    Assigns ancestry cluster using HBDSCAN on the `pc_scores_ann` from PCA.

    :param MatrixTable pca_ht: Input Table with the PCA score for each sample
    :param str pc_scores_ann: Field in `pca_ht` containing the scores
    :param int hdbscan_min_cluster_size: HDBSCAN `min_cluster_size` parameter. If not specified the smallest of 500 and 0.1*n_samples will be used.
    :param int hdbscan_min_samples: HDBSCAN `min_samples` parameter
    :return: A Table with an `ancestry_cluster` annotation containing the platform based on HDBSCAN clustering
    :rtype: Table
    """
    logger.info("Assigning ancestry based on PCA clustering")
    # Read and format data for clustering
    pca_pd = pca_ht.to_pandas()
    pc_scores = np.matrix(pca_pd[pc_scores_ann].tolist())
    logger.info(f"Assigning clusters to {len(pc_scores)} exome samples in MT...")

    # Cluster data using HDBSCAN
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=hdbscan_min_cluster_size, min_samples=hdbscan_min_samples
    )
    cluster_labels = clusterer.fit_predict(pc_scores)
    n_clusters = len(set(cluster_labels)) - (
        -1 in cluster_labels
    )  # NOTE: -1 is the label for noisy (un-classifiable) data points
    logger.info(f"Found {n_clusters} unique clusters...")

    # Add ancestry cluster to PCA DataFrame and convert back to Table
    pca_pd["ancestry_cluster"] = cluster_labels
    ht = hl.Table.from_pandas(pca_pd, key=[*pca_ht.key])
    ht = ht.annotate(ancestry_cluster=hl.int32(ht.ancestry_cluster))

    return ht


def get_array_pcs_mapped_to_exome_ids(freeze: int, n_pcs: int = 20) -> hl.Table:
    """
    Map UKBB genotype array PC file to exome IDs

    :param int freeze: One of data freezes
    :param int n_pcs: The number of PCs to keep in the Table. Default is 20.
    :return: array ancestry PC HT keyed by exome IDs
    :rtype: Table
    """
    array_pc_ht = hl.read_table(get_ukbb_array_pcs_ht_path())
    sample_map_ht = hl.read_table(array_sample_map_ht_path(freeze))
    sample_map_ht = sample_map_ht.key_by("ukbb_app_26041_id")
    array_pc_ht = array_pc_ht.key_by(s=sample_map_ht[array_pc_ht.key].s)
    array_pc_ht = array_pc_ht.filter(hl.is_defined(array_pc_ht.s))
    array_pc_ht = array_pc_ht.annotate(scores=array_pc_ht.scores[:n_pcs])

    return array_pc_ht


def main(args):
    hl.init(log="/population_pca.log", default_reference="GRCh38")

    data_source = "broad"
    freeze = args.freeze
    n_exome_pcs = args.n_exome_pcs
    n_array_pcs = args.n_array_pcs
    n_project_pcs = args.n_project_pcs
    if not args.hdbscan_min_samples:
        hdbscan_min_samples = args.hdbscan_min_cluster_size

    try:
        # Note: This code only needed to be run once and has already been run
        if args.liftover_gnomad_ancestry_loadings:
            gnomad_loadings_ht = hl.read_table(gnomad_ancestry_pca_loadings_ht_path())

            logger.info("Preparing reference genomes for liftover")
            source, target = get_liftover_genome(gnomad_loadings_ht)

            logger.info(f"Lifting data to {target.name}")
            gnomad_loadings_ht = lift_data(
                gnomad_loadings_ht,
                gnomad=False,
                data_type=None,
                path=gnomad_ancestry_loadings_liftover_path(checkpoint=True),
                rg=target,
                overwrite=args.overwrite,
            )

            logger.info("Checking SNPs for reference mismatches")
            gnomad_loadings_ht = annotate_snp_mismatch(
                gnomad_loadings_ht, data_type=None, rg=target
            )
            gnomad_loadings_ht.write(
                gnomad_ancestry_loadings_liftover_path(), overwrite=args.overwrite
            )

        if args.run_pca:
            logger.info("Running population PCA...")
            qc_mt = remove_hard_filter_samples(
                data_source,
                freeze,
                hl.read_matrix_table(qc_mt_path(data_source, freeze, ld_pruned=True)),
                gt_field="GT",
            )
            related_ht = hl.read_table(related_drop_path(data_source, freeze))
            pca_evals, pop_pca_scores_ht, pop_pca_loadings_ht = run_pca_with_relateds(
                qc_mt, related_ht, n_exome_pcs
            )
            pop_pca_scores_ht = pop_pca_scores_ht.annotate_globals(
                n_exome_pcs=n_exome_pcs
            )
            pop_pca_loadings_ht = pop_pca_loadings_ht.naive_coalesce(args.n_partitions)
            pop_pca_loadings_ht.write(
                ancestry_pca_loadings_ht_path(data_source, freeze), args.overwrite
            )
            pop_pca_scores_ht = pop_pca_scores_ht.repartition(args.n_partitions)
            pop_pca_scores_ht.write(
                ancestry_pca_scores_ht_path(data_source, freeze), args.overwrite
            )

        if args.assign_clusters:
            logger.info("Assigning PCA clustering...")
            scores_ht = hl.read_table(ancestry_pca_scores_ht_path(data_source, freeze))
            n_exome_pcs = hl.eval(scores_ht.n_exome_pcs)
            pops_ht = assign_cluster_from_pcs(
                scores_ht,
                hdbscan_min_cluster_size=args.hdbscan_min_cluster_size,
                hdbscan_min_samples=hdbscan_min_samples,
            )
            pops_ht.annotate_globals(
                n_exome_pcs=n_exome_pcs,
                hdbscan_min_cluster_size=args.hdbscan_min_cluster_size,
                hdbscan_min_samples=hdbscan_min_samples,
            )
            pops_ht = pops_ht.repartition(args.n_partitions)
            pops_ht.write(ancestry_cluster_ht_path(data_source, freeze), args.overwrite)

        if args.assign_clusters_array_pcs:
            logger.info("Loading UKBB array PC data...")
            array_pc_ht = get_array_pcs_mapped_to_exome_ids(freeze)

            logger.info("Assigning PCA clustering...")
            pops_ht = assign_cluster_from_pcs(
                array_pc_ht,
                hdbscan_min_cluster_size=args.hdbscan_min_cluster_size,
                hdbscan_min_samples=hdbscan_min_samples,
            )
            pops_ht.annotate_globals(
                n_array_pcs=n_array_pcs,
                hdbscan_min_cluster_size=args.hdbscan_min_cluster_size,
                hdbscan_min_samples=hdbscan_min_samples,
            )
            pops_ht = pops_ht.repartition(args.n_partitions)
            pops_ht.write(
                ancestry_cluster_ht_path(data_source, freeze, "array"), args.overwrite
            )

        if args.assign_clusters_joint_scratch_array_pcs:
            logger.info("Load exome scratch PC data...")
            scores_ht = hl.read_table(ancestry_pca_scores_ht_path(data_source, freeze))

            logger.info("Load UKBB array PC data...")
            array_pc_ht = get_array_pcs_mapped_to_exome_ids(freeze)
            array_pc_ht = array_pc_ht.annotate(
                scores=scores_ht[array_pc_ht.key].scores.extend(
                    array_pc_ht.scores[:n_array_pcs]
                )
            )

            logger.info("Assigning PCA clustering...")
            pops_ht = assign_cluster_from_pcs(
                array_pc_ht,
                hdbscan_min_cluster_size=args.hdbscan_min_cluster_size,
                hdbscan_min_samples=hdbscan_min_samples,
            )
            pops_ht.annotate_globals(
                hdbscan_min_cluster_size=args.hdbscan_min_cluster_size,
                hdbscan_min_samples=hdbscan_min_samples,
                n_exome_pcs=scores_ht.n_exome_pcs,
                n_array_pcs=n_array_pcs,
            )
            pops_ht = pops_ht.repartition(args.n_partitions)
            pops_ht.write(
                ancestry_cluster_ht_path(data_source, freeze, "joint"), args.overwrite,
            )

        if args.run_pc_project:
            # NOTE: I used all workers for this as it kept failing with preemptibles
            # Using the split mt led to better clustering and fewer others in Broad freeze 4
            logger.info("Running PC project...")
            mt = get_ukbb_data(data_source, freeze, split=True, adj=True)
            last_END_ht = hl.read_table(last_END_positions_ht_path(freeze))
            joint_scores_ht = project_on_gnomad_pop_pcs(mt, last_END_ht, n_project_pcs)
            joint_scores_ht = joint_scores_ht.repartition(args.n_partitions)
            joint_scores_ht.write(
                ancestry_pc_project_scores_ht_path(data_source, freeze, "joint"),
                overwrite=args.overwrite,
            )

        if args.run_rf:
            # NOTE: I needed to switch to n1-standard-16s for 300k
            logger.info("Running random forest after projection on gnomAD PCs...")
            joint_scores_ht = hl.read_table(
                ancestry_pc_project_scores_ht_path(data_source, freeze, "joint")
            )
            joint_scores_ht = joint_scores_ht.annotate(
                scores=joint_scores_ht.scores[:n_project_pcs]
            )
            joint_pops_ht, joint_pops_rf_model = assign_population_pcs(
                joint_scores_ht,
                pc_cols=joint_scores_ht.scores,
                known_col="pop_for_rf",
                min_prob=args.min_pop_prob,
            )

            scores_ht = joint_scores_ht.filter(
                hl.is_missing(joint_scores_ht.pop_for_rf)
            )
            scores_ht = scores_ht.select("scores")
            joint_pops_ht = joint_pops_ht.drop("pop_for_rf")
            scores_ht = scores_ht.annotate(pop=joint_pops_ht[scores_ht.key])
            scores_ht = scores_ht.annotate_globals(
                n_project_pcs=n_project_pcs, min_prob=args.min_pop_prob
            )

            # NOTE: Removing hard filtered samples here to avoid sample filtration pre-densify
            scores_ht = remove_hard_filter_samples(data_source, freeze, scores_ht)
            scores_ht = scores_ht.repartition(args.n_partitions)
            scores_ht = scores_ht.checkpoint(
                ancestry_pc_project_scores_ht_path(data_source, freeze),
                overwrite=args.overwrite,
            )

            logger.info(
                "Found the following sample count after population assignment (reduced to only input samples): {}".format(
                    scores_ht.aggregate(hl.agg.counter(scores_ht.pop["pop"]))
                )
            )

            logger.info("Writing out random forest model...")
            with hl.hadoop_open(
                qc_temp_data_prefix(data_source, freeze)
                + "project_gnomad_pop_rf_model.pkl",
                "wb",
            ) as out:
                pickle.dump(joint_pops_rf_model, out)

        if args.assign_hybrid_ancestry:
            logger.info(
                "Assigning hybrid population based on PCA clustering and pop assignments after projection on gnomAD PCs..."
            )
            pc_project_ht = hl.read_table(
                ancestry_pc_project_scores_ht_path(data_source, freeze)
            )
            pc_cluster_ht = hl.read_table(ancestry_cluster_ht_path(data_source, freeze))
            pc_scores_ht = hl.read_table(
                ancestry_pca_scores_ht_path(data_source, freeze)
            )

            pop_ht = pc_project_ht.select(
                gnomad_pc_project_pop=pc_project_ht.pop.pop,
                gnomad_pc_project_scores=pc_project_ht.scores,
            )

            pop_ht = pop_ht.annotate(
                HDBSCAN_pop_cluster=pc_cluster_ht[pop_ht.key].ancestry_cluster
            )
            pop_ht = pop_ht.annotate(
                hybrid_pop=hl.case()
                .when(pop_ht.HDBSCAN_pop_cluster == -1, pop_ht.gnomad_pc_project_pop)
                .default(hl.str(pop_ht.HDBSCAN_pop_cluster)),
                pop_pca_scores=pc_scores_ht[pop_ht.key].scores,
            )

            logger.info("Getting self reported ancestries...")
            ukbb_ancestry_ht = hl.read_table(
                get_ukbb_self_reported_ancestry_path(freeze)
            )
            pop_ht = pop_ht.annotate(
                self_reported_ancestry=ukbb_ancestry_ht[
                    pop_ht.key
                ].self_reported_ancestry
            )

            pop_ht = pop_ht.repartition(args.n_partitions)
            pop_ht.write(
                ancestry_hybrid_ht_path(data_source, freeze), overwrite=args.overwrite,
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(data_source, freeze))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--n_partitions",
        help="Number of desired partitions for output",
        default=5000,
        type=int,
    )
    parser.add_argument(
        "--liftover_gnomad_ancestry_loadings",
        help="Performs liftover on the gnomad joint ancestry loadings.",
        action="store_true",
    )
    parser.add_argument(
        "--hdbscan_min_samples",
        help="Minimum samples parameter for HDBSCAN. If not specified, --hdbscan_min_cluster_size is used.",
        type=int,
        required=False,
    )
    parser.add_argument(
        "--hdbscan_min_cluster_size",
        help="Minimum cluster size parameter for HDBSCAN.",
        type=int,
        default=50,
    )
    hdbscan_ancestry = parser.add_argument_group(
        "hdbscan_ancestry",
        description="Arguments relevant to infering ancestry using PCA from scratch on the exomes followed by HDBSCAN clustering.",
    )
    hdbscan_ancestry.add_argument(
        "--run_pca", help="Runs pop PCA on pruned QC MT.", action="store_true"
    )
    hdbscan_ancestry.add_argument(
        "--assign_clusters",
        help="Assigns clusters based on PCA results using HDBSCAN.",
        action="store_true",
    )
    hdbscan_ancestry.add_argument(
        "--n_exome_pcs",
        help="Number of PCs to compute on the exome data and use for HDBSCAN population inference (default: 20).",
        default=20,
        type=int,
    )
    array_ancestry = parser.add_argument_group(
        "array_ancestry",
        description="Arguments relevant to infering ancestry using array PCs and HDBSCAN clustering.",
    )
    array_ancestry.add_argument(
        "--assign_clusters_array_pcs",
        help="Assigns clusters based on UKBB genotype array PCs using HDBSCAN.",
        action="store_true",
    )
    array_ancestry.add_argument(
        "--assign_clusters_joint_scratch_array_pcs",
        help="Assigns clusters using HDBSCAN based on PCs from both the UKBB genotype array data and the exome PCs",
        action="store_true",
    )
    array_ancestry.add_argument(
        "--n_array_pcs",
        help="Number of array PCs to use in clustering (default: 20).",
        default=20,
        type=int,
    )
    gnomad_project_ancestry = parser.add_argument_group(
        "gnomad_project_ancestry",
        description="Arguments relevant to infering ancestry using projection onto gnomAD PCs followed by random forest with gnomAD known pop labels.",
    )
    gnomad_project_ancestry.add_argument(
        "--run_pc_project",
        help="Runs pc project and population assignment using gnomAD data.",
        action="store_true",
    )
    gnomad_project_ancestry.add_argument(
        "--n_project_pcs",
        help="Number of PCs to use in pc_project and population inference (default: 10).",
        default=10,
        type=int,
    )
    gnomad_project_ancestry.add_argument(
        "--run_rf",
        help="Create random forest model to assign population labels based on PCA results.",
        action="store_true",
    )
    gnomad_project_ancestry.add_argument(
        "--min_pop_prob",
        help="Minimum probability of belonging to a given population for assignment (if below, the sample is labeled as 'oth' (default: 0.5).",
        default=0.5,
        type=float,
    )
    parser.add_argument(
        "--assign_hybrid_ancestry",
        help="Assigns samples to HDBSCAN clusters where available otherwise uses pc_project assignments.",
        action="store_true",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False).",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
