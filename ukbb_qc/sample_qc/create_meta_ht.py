import argparse
import hail as hl
import logging
from ukbb_qc.load_data.utils import import_phenotype_ht
from ukbb_qc.resources.basics import array_sample_map_ht_path, phenotype_ht_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import (
    ancestry_hybrid_ht_path,
    array_sample_concordance_path,
    duplicates_ht_path,
    hard_filters_ht_path,
    platform_pca_results_ht_path,
    platform_pop_outlier_ht_path,
    qc_ht_path,
    related_drop_path,
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
    if not file_exists(f"{phenotype_ht_path(data_source, freeze)}_SUCCESS"):
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

    logger.info("Reading hard filters ht")
    left_ht = hl.read_table(hard_filters_ht_path(data_source, freeze))

    logger.info("Reading in qc ht")
    right_ht = hl.read_table(qc_ht_path(data_source, freeze))

    logger.info("Joining hard filters ht and qc ht to start making meta ht")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "outer")

    logger.info("Reading in ht mapping ukbb pharma id to ukbb application 26041 id")
    right_ht = hl.read_table(array_sample_map_ht_path(data_source, freeze))

    logger.info("Joining ukbb application 26041 id map into current join (meta ht)")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "left")

    logger.info("Reading in platform PCA ht")
    right_ht = hl.read_table(platform_pca_results_ht_path(data_source, freeze))
    right_ht = right_ht.transmute(platform_pca_scores=right_ht.scores)

    logger.info("Joining platform PCA ht into current join (meta ht)")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "outer")

    logger.info("Reading in population PC ht")
    right_ht = hl.read_table(ancestry_hybrid_ht_path(data_source, freeze))

    logger.info("Joining population PC ht into current join")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "outer")

    logger.info("Getting age information from phenotype file")
    right_ht = get_age_ht(data_source, freeze)
    left_ht = join_tables(left_ht, "s", right_ht, "s", "left")

    logger.info("Creating checkpoint")
    left_ht = left_ht.checkpoint(
        get_ht_checkpoint_path(data_source, freeze, "intermediate_ht_join"),
        overwrite=True,
    )

    logger.info("Reading in array sample concordance ht")
    right_ht = hl.read_table(array_sample_concordance_path(data_source, freeze))
    right_ht = right_ht.drop("left_col", "right_col")

    logger.info("Joining array sample concordance ht to current join")
    left_ht = left_ht.annotate(array_id=left_ht.s.split("_")[1])
    left_ht = join_tables(left_ht, "array_id", right_ht, "s", "left")

    logger.info("Reading in dups ht")
    right_ht = hl.read_table(duplicates_ht_path(data_source, freeze, args.dup_sets))
    right_ht = right_ht.drop("qc_mt_params")

    logger.info("Joining dups ht to current join")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "left")

    logger.info("Reading in related samples ht")
    related_samples_to_drop_ht = hl.read_table(related_drop_path(data_source, freeze))
    related_samples_to_drop_second_ht = related_samples_to_drop_ht.filter(
        related_samples_to_drop_ht.relationship == "2nd degree relatives"
    )
    related_samples_to_drop_first_ht = related_samples_to_drop_ht.filter(
        related_samples_to_drop_ht.relationship
        == "Parent-child" | related_samples_to_drop_ht.relationship
        == "Siblings"
    )
    related_samples_to_drop_dup_ht = related_samples_to_drop_ht.filter(
        related_samples_to_drop_ht.relationship == "Duplicate/twins"
    )

    logger.info("Annotating meta ht with related samples")
    left_ht = left_ht.annotate(
        related_filter=hl.is_defined(related_samples_to_drop_second_ht[left_ht.s])
    )
    left_ht = left_ht.annotate(
        related_dup_filter=hl.is_defined(related_samples_to_drop_dup_ht[left_ht.s])
    )
    left_ht = left_ht.annotate(
        related_first_filter=hl.is_defined(related_samples_to_drop_first_ht[left_ht.s])
    )
    left_ht = left_ht.annotate_globals(**related_samples_to_drop_ht.globals)

    logger.info("Reading in outlier ht")
    right_ht = hl.read_table(
        platform_pop_outlier_ht_path(data_source, freeze, args.pop_assignment_method)
    )

    logger.info("Joining outlier ht to current join")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "left")

    logger.info("Annotating high_quality field")
    left_ht = left_ht.annotate(
        high_quality=(
            (hl.len(left_ht.hard_filters) == 0)
            & (hl.len(left_ht.qc_metrics_filters) == 0)
        )
    )

    logger.info("Annotating releasable field")
    left_ht = left_ht.annotate(
        release=hl.if_else(
            (left_ht.s.contains("UKB") & left_ht.high_quality), True, False
        )
    )
    logger.info(
        "Releasable and non-releasable counts:"
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
