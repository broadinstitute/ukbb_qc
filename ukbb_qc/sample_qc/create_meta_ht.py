import argparse
import logging

import hail as hl

# from gnomad.slack_creds import slack_token
from gnomad.sample_qc.relatedness import (
    DUPLICATE_OR_TWINS,
    PARENT_CHILD,
    SIBLINGS,
)
from gnomad.utils.slack import slack_notifications
from ukbb_qc.resources.basics import array_sample_map_ht_path, get_checkpoint_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import (
    ancestry_hybrid_ht_path,
    array_concordance_results_path,
    hard_filters_ht_path,
    platform_pca_assignments_ht_path,
    platform_pop_outlier_ht_path,
    qc_ht_path,
    related_drop_path,
    relatedness_ht_path,
    sex_ht_path,
    meta_ht_path,
)
from ukbb_qc.resources.variant_qc import TRUTH_SAMPLES
from ukbb_qc.utils.utils import (
    get_age_ht,
    get_array_sex_ht,
    get_relatedness_set_ht,
    join_tables,
)


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_meta_ht")
logger.setLevel(logging.INFO)


def main(args):

    hl.init(log="create_meta.log", default_reference="GRCh38")

    data_source = "broad"
    freeze = args.freeze
    logging_statement = "Reading in {} and joining with meta HT"

    logger.info(
        "Reading in array map/pharma meta HT and annotating with age HT to start creating meta HT"
    )
    left_ht = hl.read_table(array_sample_map_ht_path(freeze))
    left_ht = left_ht.annotate(
        ukbb_meta=hl.struct(**left_ht.row.drop("eid_sample"))
    ).select("ukbb_meta")
    right_ht = get_age_ht(freeze)
    left_ht = join_tables(left_ht, "s", right_ht, "s", "left")
    left_ht = left_ht.annotate(
        ukbb_meta=left_ht.ukbb_meta.annotate(age=right_ht[left_ht.key].age)
    )

    logger.info(logging_statement.format("array sample concordance HT"))
    right_ht = hl.read_table(array_concordance_results_path(data_source, freeze))
    right_ht = right_ht.annotate(
        array_concordance=hl.struct(**right_ht.row.drop("s"))
    ).select("array_concordance")
    right_ht = right_ht.annotate_globals(
        array_concordance_sites_cutoffs=right_ht.globals
    ).select_globals("array_concordance_sites_cutoffs")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "left")

    logger.info(logging_statement.format("sex HT"))
    right_ht = hl.read_table(sex_ht_path(data_source, freeze))
    # Create struct for join
    right_ht = right_ht.transmute(
        #sex_imputation=hl.struct(**right_ht.row.drop("s", "array_sex"))
        sex_imputation=hl.struct(**right_ht.row.drop("s"))
    )
    right_ht = right_ht.annotate_globals(
        sex_imputation_ploidy_cutoffs=right_ht.globals
    ).select_globals("sex_imputation_ploidy_cutoffs")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "right")
    # left_ht = left_ht.transmute(
    #    ukbb_meta=left_ht.ukbb_meta.annotate(array_sex=left_ht.array_sex)
    # )
    array_sex_ht = get_array_sex_ht(freeze)
    left_ht = left_ht.annotate(
        ukbb_meta=left_ht.ukbb_meta.annotate(
            array_sex=array_sex_ht[left_ht.key].array_sex
        )
    )

    logger.info(logging_statement.format("sample QC HT"))
    right_ht = hl.read_table(qc_ht_path(data_source, freeze))
    left_ht = join_tables(left_ht, "s", right_ht, "s", "right")

    logger.info(logging_statement.format("platform PCA HT"))
    right_ht = hl.read_table(platform_pca_assignments_ht_path(data_source, freeze))
    # Put platform info into struct for join
    right_ht = right_ht.transmute(
        platform_inference=hl.struct(
            callrate_pcs=right_ht.scores, qc_platform=right_ht.qc_platform
        )
    )
    right_ht = right_ht.annotate_globals(
        platform_inference_hdbscan_parameters=right_ht.globals
    ).select_globals("platform_inference_hdbscan_parameters")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "outer")

    logger.info(logging_statement.format("population PCA HT"))
    right_ht = hl.read_table(
        "gs://broad-ukbb/broad.freeze_6/sample_qc/population_pca/hybrid_pop_assignments_10pcs.ht"
    )
    pca_ht = hl.read_table(
        "gs://broad-ukbb/broad.freeze_6/sample_qc/population_pca/pca_scores_10pcs.ht"
    )
    right_ht = right_ht.annotate_globals(**pca_ht.index_globals())
    # right_ht = hl.read_table(ancestry_hybrid_ht_path(data_source, freeze))
    right_ht = right_ht.transmute(
        gnomad_PC_project_pop_data=hl.struct(
            gnomad_PCs=right_ht.gnomad_pc_project_scores,
            gnomad_pc_project_pop=right_ht.gnomad_pc_project_pop,
        ),
        hybrid_pop_data=hl.struct(
            pop_PCs=right_ht.pop_pca_scores,
            HDBSCAN_pop_cluster=right_ht.HDBSCAN_pop_cluster,
            hybrid_pop=right_ht.hybrid_pop,
        ),
    )
    right_ht = right_ht.transmute_globals(
        population_inference_pca_metrics=hl.struct(
            n_project_pcs=right_ht.n_project_pcs,
            min_prob=right_ht.min_prob,
            n_exome_pcs=right_ht.n_exome_pcs,
        )
    )
    left_ht = join_tables(left_ht, "s", right_ht, "s", "outer")
    left_ht = left_ht.transmute(
        ukbb_meta=left_ht.ukbb_meta.annotate(
            self_reported_ancestry=left_ht.self_reported_ancestry
        )
    )

    logger.info("Creating checkpoint")
    left_ht = left_ht.checkpoint(
        get_checkpoint_path(data_source, freeze, "intermediate_ht_join", mt=False),
        overwrite=True,
    )

    logger.info(
        "Reading hard filters HT, renaming hard filters struct to sample_filters, and joining with meta HT"
    )
    right_ht = hl.read_table(hard_filters_ht_path(data_source, freeze)).select(
        "hard_filters"
    )
    left_ht = join_tables(left_ht, "s", right_ht, "s", "outer")
    left_ht = left_ht.transmute(sample_filters=left_ht.hard_filters)

    logger.info(
        "Reading in related samples to drop HT and preparing to annotate meta HT's sample_filter struct with relatedness booleans"
    )
    related_samples_to_drop_ht = hl.read_table(related_drop_path(data_source, freeze))
    relatedness_ht = get_relatedness_set_ht(
        hl.read_table(relatedness_ht_path(data_source, freeze))
    )
    related_samples_to_drop_ht = related_samples_to_drop_ht.annotate(
        relationships=relatedness_ht[related_samples_to_drop_ht.s].relationship
    )

    def _get_relationship_filter(relationship: str) -> hl.expr.builders.CaseBuilder:
        """
        Returns case statement to populate relatedness filters in sample_filters struct

        :param str relationship: Relationship to check for. One of DUPLICATE_OR_TWINS, PARENT_CHILD, or SIBLINGS.
        :return: Case statement used to population sample_filters related filter field.
        :rtype: hl.expr.builders.CaseBuilder
        """
        return (
            hl.case()
            .when(left_ht.sample_filters.hard_filtered, hl.null(hl.tbool))
            .when(
                hl.is_defined(related_samples_to_drop_ht[left_ht.key]),
                related_samples_to_drop_ht[left_ht.s].relationships.contains(
                    relationship
                ),
            )
            .default(False)
        )

    left_ht = left_ht.annotate(
        sample_filters=left_ht.sample_filters.annotate(
            related=hl.if_else(
                left_ht.sample_filters.hard_filtered,
                hl.null(hl.tbool),
                hl.is_defined(related_samples_to_drop_ht[left_ht.key]),
            ),
            duplicate=_get_relationship_filter(DUPLICATE_OR_TWINS),
            parent_child=_get_relationship_filter(PARENT_CHILD),
            sibling=_get_relationship_filter(SIBLINGS),
        )
    )
    left_ht = left_ht.annotate(
        relatedness_inference=hl.struct(
            relationships=relatedness_ht[left_ht.s].relationship
        )
    )

    logger.info("Adding relatedness globals (cutoffs)")
    left_ht = left_ht.annotate_globals(
        relatedness_inference_cutoffs=hl.struct(
            **related_samples_to_drop_ht.index_globals()
        )
    )

    logger.info(logging_statement.format("outlier HT"))
    right_ht = hl.read_table(
        platform_pop_outlier_ht_path(
            data_source,
            freeze,
            args.pop_assignment_method,
            args.platform_assignment_method,
        )
    )
    right_ht = right_ht.transmute_globals(
        outlier_detection_metrics=right_ht.qc_metrics_stats
    )
    left_ht = join_tables(left_ht, "s", right_ht, "s", "outer")
    left_ht = left_ht.transmute(
        sample_filters=left_ht.sample_filters.annotate(
            fail_n_snp=left_ht.fail_n_snp,
            fail_r_ti_tv=left_ht.fail_r_ti_tv,
            fail_r_insertion_deletion=left_ht.fail_r_insertion_deletion,
            fail_n_insertion=left_ht.fail_n_insertion,
            fail_n_deletion=left_ht.fail_n_deletion,
            fail_r_het_hom_var=left_ht.fail_r_het_hom_var,
        )
    )

    logger.info("Annotating high_quality field")
    left_ht = left_ht.annotate(
        sample_filters=left_ht.sample_filters.annotate(
            high_quality=(
                (~left_ht.sample_filters.low_callrate)
                & (~left_ht.sample_filters.ambiguous_sex)
                & (~left_ht.sample_filters.sex_aneuploidy)
                & (~left_ht.sample_filters.low_coverage)
                & (hl.len(left_ht.qc_metrics_filters) == 0)
            )
        )
    )
    left_ht = left_ht.drop("qc_metrics_filters")

    logger.info("Annotating releasable field")
    # Control samples: CHMI_CHMI3_Nex1 and Coriell_NA12878_NA12878
    left_ht = left_ht.annotate(
        sample_filters=left_ht.sample_filters.annotate(
            release=hl.if_else(
                (
                    hl.is_defined(left_ht.ukbb_meta.batch)
                    & left_ht.sample_filters.high_quality
                ),
                True,
                False,
            )
        )
    )

    logger.info("Annotating control samples")
    left_ht = left_ht.annotate(
        sample_filters=left_ht.sample_filters.annotate(
            control=(hl.literal(TRUTH_SAMPLES).contains(left_ht.s))
        )
    )
    logger.info(
        "Release and control sample counts:"
        f"{left_ht.aggregate(hl.struct(release=hl.agg.count_where(left_ht.sample_filters.release), control=hl.agg.count_where(left_ht.sample_filters.control)))}"
    )

    logger.info("Removing duplicate samples and writing out meta ht")
    left_ht = left_ht.distinct()
    left_ht = left_ht.repartition(args.n_partitions)
    left_ht = left_ht.checkpoint(
        meta_ht_path(data_source, freeze), overwrite=args.overwrite
    )
    left_ht.summarize()
    logger.info(f"Final count: {left_ht.count()}")
    logger.info("Complete")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="This script creates the sample meta table for UKBB data"
    )
    parser.add_argument(
        "--n_partitions",
        help="Desired number of partitions for output HT",
        default=5000,
        type=int,
    )
    parser.add_argument(
        "-f", "--freeze", help="Current freeze", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--pop_assignment_method",
        help="Population assignment method to use for outlier stratification",
        default="hybrid_pop",
        choices=["gnomad_pc_project_pop", "HDBSCAN_pop_cluster", "hybrid_pop"],
    )
    parser.add_argument(
        "--platform_assignment_method",
        help="Platform assignment method to use for outlier stratification",
        default="batch",
        choices=["batch", "qc_platform"],
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to",
    )
    args = parser.parse_args()

    if args.slack_channel:
        slack_notifications(
            "xoxb-3114021723-1092439063975-G6icxweKHcX1nmfG5FFqbErn", args.slack_channel
        )

    main(args)
