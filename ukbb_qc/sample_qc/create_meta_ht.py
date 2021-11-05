import argparse
import logging

import hail as hl

from gnomad.sample_qc.relatedness import (
    DUPLICATE_OR_TWINS,
    PARENT_CHILD,
    SIBLINGS,
)
from gnomad.utils.slack import slack_notifications
from ukbb_qc.resources.basics import (
    array_sample_map_ht_path,
    check_dups_to_remove,
    excluded_samples_path,
    get_checkpoint_path,
    known_dups_ht_path,
    logging_path,
    pan_ancestry_ht_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import (
    ancestry_hybrid_ht_path,
    array_concordance_results_path,
    hard_filters_ht_path,
    platform_pop_outlier_ht_path,
    qc_ht_path,
    related_drop_path,
    relatedness_ht_path,
    sex_ht_path,
)
from ukbb_qc.resources.variant_qc import SYNDIP, NA12878
from ukbb_qc.slack_creds import slack_token
from ukbb_qc.utils.utils import (
    get_age_ht,
    get_relationship_filter_expr,
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

    hl.init(log="/create_meta.log", default_reference="GRCh38")

    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)

    try:
        logging_statement = "Reading in {} and joining with meta HT"
        logger.info(
            "Reading in array map/pharma meta HT and annotating with age HT to start creating meta HT"
        )
        left_ht = hl.read_table(
            array_sample_map_ht_path(freeze), _n_partitions=args.n_partitions
        )
        left_ht = left_ht.annotate(
            ukb_meta=hl.struct(**left_ht.row.drop("s", "eid_sample"))
        ).select("ukbb_meta")
        right_ht = get_age_ht(freeze)
        left_ht = left_ht.annotate(
            ukbb_meta=left_ht.ukbb_meta.annotate(age=right_ht[left_ht.key].age)
        )

        logger.info("Filtering to samples with defined batch...")
        left_ht = left_ht.filter(hl.is_defined(left_ht.ukb_meta.batch))

        logger.info("Removing samples with withdrawn consent...")
        logger.info("Sample count before filtering: %i", left_ht.count())
        withdrawn_ht = hl.import_table(excluded_samples_path(), no_header=True).key_by(
            "f0"
        )

        left_ht = left_ht.filter(hl.is_missing(withdrawn_ht[left_ht.ukbb_app_26041_id]))
        logger.info("Sample count after filtering: %i", left_ht.count())

        logger.info(logging_statement.format("array sample concordance HT"))
        left_ht = left_ht.annotate(ukbb_id=left_ht.ukbb_meta.ukbb_app_26041_id)
        right_ht = hl.read_table(array_concordance_results_path(*tranche_data))
        right_ht = right_ht.annotate(
            array_concordance=hl.struct(**right_ht.row.drop("s"))
        ).select("array_concordance")
        right_ht = right_ht.annotate_globals(
            array_concordance_sites_cutoffs=right_ht.globals
        ).select_globals("array_concordance_sites_cutoffs")
        left_ht = join_tables(left_ht, "ukbb_id", right_ht, "s", "left")
        left_ht = left_ht.key_by("s").drop("ukbb_id")

        logger.info(logging_statement.format("sex HT"))
        right_ht = hl.read_table(sex_ht_path(*tranche_data))
        # Create struct for join
        right_ht = right_ht.transmute(
            sex_imputation=hl.struct(**right_ht.row.drop("s", "array_sex"))
        ).select("sex_imputation", "array_sex")
        right_ht = right_ht.annotate_globals(
            sex_imputation_ploidy_cutoffs=right_ht.globals
        ).select_globals("sex_imputation_ploidy_cutoffs")
        left_ht = join_tables(left_ht, "s", right_ht, "s", "right")
        left_ht = left_ht.transmute(
            ukbb_meta=left_ht.ukbb_meta.annotate(array_sex=left_ht.array_sex)
        )

        logger.info(logging_statement.format("sample QC HT"))
        right_ht = hl.read_table(qc_ht_path(*tranche_data))
        left_ht = join_tables(left_ht, "s", right_ht, "s", "right")

        logger.info(logging_statement.format("population PCA HT"))
        logger.info("Keeping only self reported ancestry information...")
        right_ht = hl.read_table(ancestry_hybrid_ht_path(*tranche_data))
        right_ht = right_ht.select("self_reported_ancestry")
        left_ht = join_tables(left_ht, "s", right_ht, "s", "outer")
        left_ht = left_ht.transmute(
            ukbb_meta=left_ht.ukbb_meta.annotate(
                self_reported_ancestry=left_ht.self_reported_ancestry
            )
        )

        logger.info(logging_statement.format("pan-ancestry HT"))
        right_ht = hl.read_table(pan_ancestry_ht_path())
        left_ht = left_ht.annotate(
            pan_ancestry_pop=right_ht[left_ht.ukbb_meta.ukbb_app_26041_id].pop
        )

        logger.info("Creating checkpoint")
        left_ht = left_ht.checkpoint(
            get_checkpoint_path(*tranche_data, "intermediate_ht_join", mt=False),
            overwrite=True,
        )

        logger.info(
            "Reading hard filters HT, renaming hard filters struct to sample_filters, and joining with meta HT"
        )
        right_ht = hl.read_table(hard_filters_ht_path(*tranche_data)).select(
            "hard_filters"
        )
        left_ht = join_tables(left_ht, "s", right_ht, "s", "outer")
        left_ht = left_ht.transmute(sample_filters=left_ht.hard_filters)

        logger.info(
            "Reading in related samples to drop HT and preparing to annotate meta HT's sample_filter struct with relatedness booleans"
        )
        related_samples_to_drop_ht = hl.read_table(related_drop_path(*tranche_data))
        relatedness_ht = get_relatedness_set_ht(
            hl.read_table(relatedness_ht_path(*tranche_data))
        )
        related_samples_to_drop_ht = related_samples_to_drop_ht.annotate(
            relationships=relatedness_ht[related_samples_to_drop_ht.s].relationships
        )

        # Annotating meta HT with related filter booleans
        # Any sample that is hard filtered will have missing values for these bools
        # Any sample that was filtered for relatedness will have True for sample_filters.related
        # If a filtered related sample had a relationship with a higher degree than second-degree (duplicate, parent-child, sibling),
        # that filter will also be True
        left_ht = left_ht.annotate(
            sample_filters=left_ht.sample_filters.annotate(
                related=hl.if_else(
                    left_ht.sample_filters.hard_filtered,
                    hl.null(hl.tbool),
                    hl.is_defined(related_samples_to_drop_ht[left_ht.key]),
                ),
                duplicate=get_relationship_filter_expr(
                    left_ht.sample_filters.hard_filtered,
                    DUPLICATE_OR_TWINS,
                    related_samples_to_drop_ht[left_ht.key].relationships,
                ),
                parent_child=get_relationship_filter_expr(
                    left_ht.sample_filters.hard_filtered,
                    PARENT_CHILD,
                    related_samples_to_drop_ht[left_ht.key].relationships,
                ),
                sibling=get_relationship_filter_expr(
                    left_ht.sample_filters.hard_filtered,
                    SIBLINGS,
                    related_samples_to_drop_ht[left_ht.key].relationships,
                ),
            )
        )
        left_ht = left_ht.annotate(
            relatedness_inference_relationships=relatedness_ht[left_ht.s].relationships,
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
                *tranche_data,
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

        logger.info("Removing control samples...")
        truth_samples = hl.literal([SYNDIP, NA12878])
        left_ht = left_ht.filter(~truth_samples.contains(left_ht.s))

        logger.info("Annotating releasable field")
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

        logger.info(
            "Release sample counts:"
            f"{left_ht.aggregate(hl.struct(release=hl.agg.count_where(left_ht.sample_filters.release)))}"
        )

        logger.info("Removing duplicate samples and writing out meta ht")
        left_ht = left_ht.distinct()
        known_dups_ht = hl.read_table(known_dups_ht_path(freeze))
        ids_to_remove = known_dups_ht.aggregate(
            hl.agg.collect(known_dups_ht["Sample Name - ID1"]), _localize=False
        )
        # Check number of samples to remove
        num_ids = check_dups_to_remove(ids_to_remove, left_ht)
        logger.info(f"Removing {num_ids} samples...")
        left_ht = left_ht.filter(~ids_to_remove.contains(left_ht.s))

        logger.info("Rekeying table by UKBB ID...")
        left_ht = left_ht.key_by()
        left_ht = left_ht.key_by(s=left_ht.ukbb_meta.ukbb_app_26041_id)

        logger.info(
            "Writing new meta HT to temp (to avoid overwriting meta HT used for QC (written to `meta_ht_path`)..."
        )
        left_ht = left_ht.checkpoint(
            get_checkpoint_path(*tranche_data, "ror_meta", mt=False),
            overwrite=args.overwrite,
        )
        left_ht.summarize()
        logger.info(f"Final count: {left_ht.count()}")
        logger.info("Complete")

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(*tranche_data))


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
        "--platform_assignment_method",
        help="Platform assignment method to use for outlier stratification",
        default="batch",
        choices=["batch", "qc_platform"],
    )
    parser.add_argument(
        "--pop_assignment_method",
        help="Population assignment method to use for outlier stratification",
        default="hybrid_pop_data",
        choices=["gnomad_pc_project_pop_data", "hybrid_pop_data"],
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
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
