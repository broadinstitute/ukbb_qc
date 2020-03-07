import argparse
import hail as hl
import logging
from gnomad.utils.generic import file_exists
from ukbb_qc.load_data.utils import import_phenotype_ht
from ukbb_qc.resources.basics import array_sample_map_ht_path, phenotype_ht_path,  get_checkpoint_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import (
    ancestry_hybrid_ht_path,
    array_sample_concordance_path,
    hard_filters_ht_path,
    platform_pca_results_ht_path,
    platform_pop_outlier_ht_path,
    qc_ht_path,
    related_drop_path,
    sex_ht_path,
    meta_ht_path,
)
from ukbb_qc.utils.utils import join_tables


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_meta_ht")
logger.setLevel(logging.INFO)


def get_age_ht(data_source: str, freeze: int) -> hl.Table:
    """
    Pull age information from UKBB phenotype file

    :param str data_source: One of "regeneron" or "broad"
    :param int freeze: One of the data freezes
    :return: Table with age at recruitment per sample
    :rtype: Table
    """
    # Read in phenotype table
    if not file_exists(f"{phenotype_ht_path()}_SUCCESS"):
        import_phenotype_ht()
    phenotype_ht = hl.read_table(phenotype_ht_path()).select("f.21022.0.0")

    # Re-key phenotype table to UKBB ID using array sample map table and return
    sample_map_ht = hl.read_table(array_sample_map_ht_path(data_source, freeze))
    sample_map_ht = sample_map_ht.key_by("ukbb_app_26041_id")
    phenotype_ht = phenotype_ht.key_by(s=sample_map_ht[phenotype_ht.key].s)
    phenotype_ht = phenotype_ht.rename({"f.21022.0.0": "age"}).select("age")
    return phenotype_ht


def main(args):

    hl.init(log="create_meta.log", default_reference="GRCh38")

    data_source = args.data_source
    if args.freeze:
        freeze = args.freeze
    else:
        freeze = CURRENT_FREEZE

    logger.info(
        "Reading in HT mapping ukbb pharma (exome) id to ukbb application 26041 (array) id"
    )
    left_ht = hl.read_table(array_sample_map_ht_path(data_source, freeze))
    left_ht = left_ht.key_by(eid_sample=left_ht.s)
    left_ht = left_ht.annotate(
        pharma_meta=hl.struct(
            ukbb_app_26041_id=left_ht.ukbb_app_26041_id,
            eid_sample=left_ht.eid_sample,
            batch=left_ht.batch,
            batch_num=left_ht.batch_num,
        )
    ).select("pharma_meta")

    logger.info("Getting age information from phenotype file")
    right_ht = get_age_ht(data_source, freeze)
    right_ht = right_ht.key_by(array_id=right_ht.s).select("age")

    logger.info(
        "Joining array map/pharma meta HT (left) and age HT (right) with a right join to start creating meta HT"
    )
    left_ht = join_tables(left_ht, "eid_sample", right_ht, "array_id", "right")

    logger.info("Reading in array sample concordance ht")
    right_ht = hl.read_table(array_sample_concordance_path(data_source, freeze))
    # Add array ID and struct for join
    right_ht = right_ht.annotate(array_id=right_ht.s.split("_")[1])
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
    left_ht = join_tables(left_ht, "array_id", right_ht, "array_id", "left")
    left_ht = left_ht.drop("array_id", "eid_sample")

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

    logger.info("Reading in qc HT")
    right_ht = hl.read_table(qc_ht_path(data_source, freeze))

    logger.info("Joining sample QC HT with meta HT")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "right")

    logger.info("Reading in platform PCA ht")
    right_ht = hl.read_table(platform_pca_results_ht_path(data_source, freeze))
    # Put platform info into struct for join
    right_ht = right_ht.transmute(
        platform_inference=hl.struct(
            callrate_pcs=right_ht.scores, qc_platform=right_ht.qc_platform
        )
    )

    logger.info("Joining platform HT with meta HT")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "outer")

    logger.info("Reading in population PC ht")
    right_ht = hl.read_table(ancestry_hybrid_ht_path(data_source, freeze))
    # Put population info into structs for join
    right_ht = right_ht.transmute(
        gnomad_PC_project_pop_data=hl.struct(
            gnomad_PCs=hl.array(
                [right_ht[f"gnomad_pc_project_PC{i + 1}"] for i in range(n_pcs)]
            ),
            gnomad_pc_project_pop=right_ht.gnomad_pc_project_pop,
        ),
        hybrid_pop_data=hl.struct(
            pop_PCs=hl.array(
                [right_ht[f"pop_pca_PC{i + 1}"] for i in range(n_pcs)]
            ),
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
    related_samples_to_drop_second_ht = related_samples_to_drop_ht.filter(
        related_samples_to_drop_ht.relationship == "2nd degree relatives"
    )
    related_samples_to_drop_pc_ht = related_samples_to_drop_ht.filter(
        related_samples_to_drop_ht.relationship == "Parent-child"
    )
    related_samples_to_drop_sib_ht = related_samples_to_drop_ht.filter(
        related_samples_to_drop_ht.relationship == "Siblings"
    )
    related_samples_to_drop_dup_ht = related_samples_to_drop_ht.filter(
        related_samples_to_drop_ht.relationship == "Duplicate/twins"
    )

    logger.info("Annotating sample_filter struct with relatedness booleans")
    left_ht = left_ht.annotate(
        sample_filters=left_ht.sample_filters.annotate(
            related=hl.is_defined(related_samples_to_drop_second_ht[left_ht.s]),
            duplicate=hl.is_defined(related_samples_to_drop_dup_ht[left_ht.s]),
            parent_child=hl.is_defined(related_samples_to_drop_pc_ht[left_ht.s]),
            sibling=hl.is_defined(related_samples_to_drop_pc_ht[left_ht.s]),
        )
    )

    logger.info("Reading in outlier HT")
    right_ht = hl.read_table(
        platform_pop_outlier_ht_path(data_source, freeze, args.pop_assignment_method)
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
                (left_ht.s.contains("UKB") & left_ht.high_quality), True, False
            )
        )
    )
    logger.info(
        "Release and control sample counts:"
        f"{left_ht.aggregate(hl.struct(release=hl.agg.count_where(left_ht.release), control=hl.agg.count_where(~left_ht.release)))}"
    )

    logger.info("Writing out meta ht")
    left_ht.write(meta_ht_path(data_source, freeze), overwrite=args.overwrite)
    logger.info(f"Final count: {left_ht.count()}")
    logger.info("Complete")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script applies hard filters to UKBB data"
    )
    parser.add_argument(
        "-s",
        "--data_source",
        help="Data source",
        choices=["regeneron", "broad"],
        default="broad",
    )
    parser.add_argument(
        "-f", "--freeze", help="Current freeze #", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--n_pcs",
        help="Number of PCs used in population inference (default: 10)",
        default=10,
        type=int,
    )
    parser.add_argument(
        "-d",
        "--dup_sets",
        help="Dup sets for duplicate ht table",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--pop_assignment_method",
        help="Population assignment method to use for outlier stratification",
        default="hybrid_pop",
        choices=["gnomad_pc_project_pop", "HDBSCAN_pop_cluster", "hybrid_pop"],
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    args = parser.parse_args()

    main(args)
