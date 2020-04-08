import argparse
import logging

import hail as hl

from gnomad.sample_qc.relatedness import (
    DUPLICATE_OR_TWINS,
    PARENT_CHILD,
    SECOND_DEGREE_RELATIVES,
    SIBLINGS,
)
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
    sex_ht_path,
    meta_ht_path,
)
from ukbb_qc.resources.variant_qc import TRUTH_SAMPLES
from ukbb_qc.utils.utils import get_age_ht, join_tables


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

    logger.info(
        "Reading in HT mapping ukbb pharma (exome) id to ukbb application 26041 (array) id"
    )
    left_ht = hl.read_table(array_sample_map_ht_path(freeze))
    left_ht = left_ht.annotate(
        pharma_meta=hl.struct(
            ukbb_app_26041_id=left_ht.ukbb_app_26041_id,
            batch=left_ht.batch,
            batch_num=left_ht.batch_num,
        )
    ).select("pharma_meta")

    logger.info("Getting age information from phenotype file")
    right_ht = get_age_ht(data_source, freeze)

    logger.info(
        "Joining array map/pharma meta HT (left) and age HT (right) with a right join to start creating meta HT"
    )
    left_ht = join_tables(left_ht, "s", right_ht, "s", "left")

    logger.info("Reading in array sample concordance HT")
    right_ht = hl.read_table(array_concordance_results_path(data_source, freeze))
    right_ht = right_ht.tramsmute(
        array_concordance=hl.struct(
            concordance=right_ht.concordance,
            n_discordant=right_ht.n_discordant,
            num_gt_con_non_ref=right_ht.num_gt_con_non_ref,
            num_gt_non_ref=right_ht.num_gt_non_ref,
            prop_gt_con_non_ref=right_ht.prop_gt_con_non_ref,
        )
    ).select("array_concordance")

    logger.info("Joining array sample concordance HT with meta HT")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "left")

    logger.info("Reading in sex HT")
    right_ht = hl.read_table(sex_ht_path(data_source, freeze))
    # Create struct for join
    right_ht = right_ht.transmute(
        sex_imputation=hl.struct(
            is_female=right_ht.is_female,
            f_stat=right_ht.f_stat,
            n_called=right_ht.n_called,
            expected_homs=right_ht.expected_homs,
            observed_homs=right_ht.observed_homs,
            chr20_mean_dp=right_ht.chr20_mean_dp,
            chrX_mean_dp=right_ht.chrX_mean_dp,
            chrY_mean_dp=right_ht.chrY_mean_dp,
            chrX_ploidy=right_ht.chrX_ploidy,
            chrY_ploidy=right_ht.chrY_ploidy,
            sex_karyotype=right_ht.sex_karyotype,
        )
    )

    logger.info("Joining sex HT with meta HT")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "right")

    logger.info("Reading in sample QC HT")
    right_ht = hl.read_table(qc_ht_path(data_source, freeze))

    logger.info("Joining sample QC HT with meta HT")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "right")

    logger.info("Reading in platform PCA HT")
    right_ht = hl.read_table(platform_pca_assignments_ht_path(data_source, freeze))
    # Put platform info into struct for join
    right_ht = right_ht.transmute(
        platform_inference=hl.struct(
            callrate_pcs=right_ht.scores, qc_platform=right_ht.qc_platform
        )
    )

    logger.info("Joining platform HT with meta HT")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "outer")

    logger.info("Reading in population PC HT")
    right_ht = hl.read_table(ancestry_hybrid_ht_path(data_source, freeze))
    # Put population info into structs for join
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

    logger.info("Joining population PC HT with meta HT")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "outer")

    logger.info("Creating checkpoint")
    left_ht = left_ht.checkpoint(
        get_checkpoint_path(data_source, freeze, "intermediate_ht_join", mt=False),
        overwrite=True,
    )

    logger.info("Reading hard filters HT")
    right_ht = hl.read_table(hard_filters_ht_path(data_source, freeze))

    logger.info("Joining hard filters HT with meta HT")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "outer")

    logger.info("Renaming hard_filters struct to sample_filter struct")
    left_ht = left_ht.transmute(sample_filters=left_ht.hard_filters)

    logger.info("Reading in related samples HT")
    related_samples_to_drop_ht = hl.read_table(related_drop_path(data_source, freeze))

    logger.info("Annotating sample_filter struct with relatedness booleans")
    left_ht = left_ht.annotate(
        sample_filters=left_ht.sample_filters.annotate(
            related=hl.if_else(
                related_samples_to_drop_ht.relationship[left_ht.s]
                == SECOND_DEGREE_RELATIVES,
                True,
                False,
            ),
            duplicate=hl.if_else(
                related_samples_to_drop_ht.relationship[left_ht.s]
                == DUPLICATE_OR_TWINS,
                True,
                False,
            ),
            parent_child=hl.if_else(
                related_samples_to_drop_ht.relationship[left_ht.s] == PARENT_CHILD,
                True,
                False,
            ),
            sibling=hl.if_else(
                related_samples_to_drop_ht.relationship[left_ht.s] == SIBLINGS,
                True,
                False,
            ),
        )
    )

    logger.info("Adding relatedness globals (cutoffs)")
    left_ht = left_ht.annotate_globals(**related_samples_to_drop_ht.globals)

    logger.info("Reading in outlier HT")
    right_ht = hl.read_table(
        platform_pop_outlier_ht_path(
            data_source,
            freeze,
            args.pop_assignment_method,
            args.platform_assignment_method,
        )
    )

    logger.info("Joining outlier ht to current join")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "left")
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
                (hl.is_defined(left_ht.pharma_meta.batch) & left_ht.high_quality),
                True,
                False,
            )
        )
    )

    logger.info("Annotating control samples")
    left_ht = left_ht.annotate(control=(hl.literal(TRUTH_SAMPLES).contains(left_ht.s)))
    logger.info(
        "Release and control sample counts:"
        f"{left_ht.aggregate(hl.struct(release=hl.agg.count_where(left_ht.release), control=hl.agg.count_where(left_ht.control)))}"
    )

    logger.info("Writing out meta ht")
    left_ht = left_ht.checkpoint(
        meta_ht_path(data_source, freeze), overwrite=args.overwrite
    )
    logger.info(f"Final count: {left_ht.count()}")
    logger.info("Complete")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="This script creates the sample meta table for UKBB data"
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
    args = parser.parse_args()
    main(args)
