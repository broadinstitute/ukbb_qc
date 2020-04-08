import argparse
import logging
from typing import Callable, Tuple

import hail as hl

from gnomad.sample_qc.relatedness import (
    DUPLICATE_OR_TWINS,
    explode_duplicate_samples_ht,
    get_duplicated_samples_ht,
    get_duplicated_samples,
    get_relationship_expr,
    infer_families,
    PARENT_CHILD,
    SECOND_DEGREE_RELATIVES,
    SIBLINGS,
)
from gnomad.utils.slack import try_slack
from ukbb_qc.resources.basics import get_checkpoint_path, logging_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import (
    duplicates_ht_path,
    inferred_ped_path,
    related_drop_path,
    relatedness_ht_path,
    relatedness_pca_scores_ht_path,
    sex_ht_path,
    qc_ht_path,
    qc_mt_path,
)
from ukbb_qc.utils.utils import remove_hard_filter_samples


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("relatedness")
logger.setLevel(logging.INFO)


# Note: Will need to change this to work with gnomAD, we will want to default to gnomAD ranking, unless tied pair is gnomAD, UKBB pair
def rank_related_samples(
    relatedness_ht: hl.Table, qc_ht: hl.Table
) -> Tuple[
    hl.Table,
    Callable[[hl.expr.Expression, hl.expr.Expression], hl.expr.NumericExpression],
]:
    """
    Rank related samples based on their mean depths. 

    Ranking is used when determining which sample to filter from a pair of related samples.

    :param Table relatedness_ht: Table of samples with relatedness results from pc_relate.
    :param Table qc_ht: Table of samples and their sample QC metrics, including mean depth.
    :return: Tuple of the relatedness table, annotated with ranking, and a tie breaker 
    :rtype: Tuple[hl.Table, Callable[[hl.expr.Expression, hl.expr.Expression], hl.expr.NumericExpression]]
    """

    def annotate_related_pairs(related_pairs: hl.Table, index_col: str) -> hl.Table:
        """
        Annotates each sample in a related pair with mean depth

        Returns Table with mean depth annotation and unkeys Table

        :param Table related_pairs: Table of samples with relatedness results from pc_relate
        :param str index_col: Name of column containing samples to be annotated
        :return: Table annotated with sample mean depth
        :rtype: hl.Table
        """
        related_pairs = related_pairs.key_by(**related_pairs[index_col])
        return related_pairs.annotate(
            **{
                index_col: related_pairs[index_col].annotate(
                    dp_mean=hl.or_else(
                        qc_ht[related_pairs.key].sample_qc.dp_stats.mean, -1.0
                    )
                )
            }
        ).key_by()

    relatedness_ht = annotate_related_pairs(relatedness_ht, "i")
    relatedness_ht = annotate_related_pairs(relatedness_ht, "j")

    def tie_breaker(
        l: hl.expr.Float64Expression, r: hl.expr.Float64Expression
    ) -> hl.expr.Float64Expression:
        """
        Tie breaker given to maximal independent set (MIS) that determines which related sample to remove.

        .. note::

            This tie breaker is r - l (instead of l - r) because we want to keep the sample with the greater depth.
            MIS removes the largest node.
            MIS removes the right node when the tie breaker returns a negative value.
            MIS removes the left node when the tie breaker returns a positive value.
            By using r - l, we make the sample with the smaller depth the larger node.
            This removes the sample with lower depth.

        :param hl.expr.Float64Expression l: Mean depth of a sample (left) in the related pair
        :param hl.expr.Float64Expression r: Mean depth of the second sample (right) in the related pair
        :return: The result of subtracting the left depth from the right depth.
            Returns a negative number with l > r. This mean MIS will remove the right (lower depth) node.
            Returns a positive number when l < r. This means MIS will remove the left (lower depth) node.
        :rtype: hl.expr.Float64Expression
        """
        return r.dp_mean - l.dp_mean

    return relatedness_ht, tie_breaker


def filter_related_samples(
    relatedness_ht: hl.Table, qc_ht: hl.Table, relationship: str,
) -> hl.Table:
    """
    Filters Table to keep only samples with input relationship.

    .. note::
        Input Table must be annotated with relationship (from `get_relationship_expr`).

    :param Table relatedness_ht: Table of samples with relatedness results from `pc_relate` annotated with `get_relationship_expr`.
    :param Table qc_ht: Table of samples and their sample QC metrics, including mean depth.
    :param str relationship: Desired relationship type to filter.
    :return: Filtered Table of samples with kinship above input cutoff
    :rtype: hl.Table
    """
    relatedness_ht = relatedness_ht.filter(relatedness_ht.relationship == relationship)
    related_pairs_ht, related_pairs_tie_breaker = rank_related_samples(
        relatedness_ht, qc_ht,
    )

    related_samples_to_drop_ht = hl.maximal_independent_set(
        related_pairs_ht.i,
        related_pairs_ht.j,
        keep=False,
        tie_breaker=related_pairs_tie_breaker,
    )
    related_samples_to_drop_ht = related_samples_to_drop_ht.key_by()
    related_samples_to_drop_ht = related_samples_to_drop_ht.select(
        **related_samples_to_drop_ht.node
    )
    related_samples_to_drop_ht = related_samples_to_drop_ht.key_by("s")
    return related_samples_to_drop_ht


def main(args):
    hl.init(
        log="/relatedness.log", default_reference="GRCh38",
    )

    data_source = "broad"
    freeze = args.freeze

    try:
        if not args.skip_pc_relate:
            logger.info("Running PCA for PC-Relate...")
            pruned_qc_mt = remove_hard_filter_samples(
                data_source,
                freeze,
                hl.read_matrix_table(qc_mt_path(data_source, freeze, ld_pruned=True)),
                gt_field="GT",
            ).unfilter_entries()
            eig, scores, _ = hl.hwe_normalized_pca(
                pruned_qc_mt.GT, k=10, compute_loadings=False
            )
            scores.write(
                relatedness_pca_scores_ht_path(data_source, freeze), args.overwrite
            )

            logger.info("Running PC-Relate...")
            # NOTE: This needs SSDs on your workers (for the temp files) and no preemptible workers while the BlockMatrix writes
            relatedness_ht = hl.pc_relate(
                pruned_qc_mt.GT,
                min_individual_maf=args.min_individual_maf,
                scores_expr=scores[pruned_qc_mt.col_key].scores,
                block_size=4096,
                min_kinship=args.min_emission_kinship,
                statistics="all",
            )
            relatedness_ht = relatedness_ht.checkpoint(
                get_checkpoint_path(data_source, freeze, name="pc_relate_temp")
            )

            logger.info("Annotating PC relate results with relationships...")
            relatedness_ht = relatedness_ht.annotate(
                relationship=get_relationship_expr(
                    kin_expr=relatedness_ht.kin,
                    ibd0_expr=relatedness_ht.ibd0,
                    ibd1_expr=relatedness_ht.ibd1,
                    ibd2_expr=relatedness_ht.ibd2,
                    first_degree_kin_thresholds=args.first_degree_kin_thresholds,
                    second_degree_min_kin=args.second_degree_kin_cutoff,
                    ibd0_0_max=args.ibd0_0_max,
                )
            )
            relatedness_ht = relatedness_ht.annotate_globals(
                min_individual_maf=args.min_individual_maf,
                min_emission_kinship=args.min_emission_kinship,
                ibd0_0_max=args.ibd0_0_max,
                second_degree_kin_cutoff=args.second_degree_kin_cutoff,
                first_degree_kin_thresholds=args.first_degree_kin_thresholds,
            )
            relatedness_ht.write(
                relatedness_ht_path(data_source, freeze), args.overwrite
            )

        if not args.skip_filter_dups:
            logger.info("Filtering duplicate samples...")
            sample_qc_ht = remove_hard_filter_samples(
                data_source, freeze, hl.read_table(qc_ht_path(data_source, freeze))
            )
            samples_rankings_ht = sample_qc_ht.select(
                rank=(-1 * sample_qc_ht.sample_qc.dp_stats.mean)
            )
            relatedness_ht = hl.read_table(relatedness_ht_path(data_source, freeze))
            dup_samples = get_duplicated_samples(relatedness_ht)
            if len(dup_samples) > 0:
                dups_ht = get_duplicated_samples_ht(dup_samples, samples_rankings_ht)
                dups_ht.write(
                    duplicates_ht_path(data_source, freeze, dup_sets=True),
                    overwrite=args.overwrite,
                )
                dups_ht = explode_duplicate_samples_ht(dups_ht)
                dups_ht = sample_qc_ht.select(
                    duplicate=hl.is_defined(dups_ht[sample_qc_ht.key].dup_filtered)
                    & dups_ht[sample_qc_ht.key].dup_filtered,
                )
            else:
                dups_ht = sample_qc_ht.select(duplicate=False)
            dups_ht.write(
                duplicates_ht_path(data_source, freeze), overwrite=args.overwrite
            )

        if not args.skip_infer_families:
            logger.info("Inferring families...")
            dups_ht = hl.read_table(duplicates_ht_path(data_source, freeze))
            ped = infer_families(
                hl.read_table(relatedness_ht_path(data_source, freeze)),
                hl.read_table(sex_ht_path(data_source, freeze)),
                dups_ht,
            )
            ped.write(inferred_ped_path(data_source, freeze))

        if not args.skip_filter_related_samples:
            logger.info("Filtering related samples...")
            relatedness_ht = hl.read_table(relatedness_ht_path(data_source, freeze))
            qc_ht = hl.read_table(qc_ht_path(data_source, freeze))

            # Filter second degree samples
            related_samples_to_drop_second_deg_ht = filter_related_samples(
                relatedness_ht, qc_ht, SECOND_DEGREE_RELATIVES
            )

            # Filter duplicate samples
            related_samples_to_drop_dup_ht = filter_related_samples(
                relatedness_ht, qc_ht, DUPLICATE_OR_TWINS
            )

            # Filter first degree samples
            related_samples_to_drop_pc_ht = filter_related_samples(
                relatedness_ht, qc_ht, PARENT_CHILD
            )
            related_samples_to_drop_sib_ht = filter_related_samples(
                relatedness_ht, qc_ht, SIBLINGS
            )

            # Combine first, second, and duplicate sample tables and annotate cutoffs
            related_samples_to_drop_ht = related_samples_to_drop_second_deg_ht.union(
                related_samples_to_drop_dup_ht,
                related_samples_to_drop_pc_ht,
                related_samples_to_drop_sib_ht,
            )
            related_samples_to_drop_ht = related_samples_to_drop_ht.annotate_globals(
                min_individual_maf=args.min_individual_maf,
                min_emission_kinship=args.min_emission_kinship,
                ibd0_0_max=args.ibd0_0_max,
                second_degree_kin_cutoff=args.second_degree_kin_cutoff,
                first_degree_kin_thresholds=args.first_degree_kin_thresholds,
            )

            related_samples_to_drop_ht = related_samples_to_drop_ht.checkpoint(
                related_drop_path(data_source, freeze), overwrite=args.overwrite
            )
            logger.info(
                f"{related_samples_to_drop_ht.filter(related_samples_to_drop_ht.relationship == SECOND_DEGREE_RELATIVES).count()}"
                "second degree samples flagged in callset using maximal independent set"
            )
            logger.info(
                f"{related_samples_to_drop_ht.filter(related_samples_to_drop_ht.relationship == SIBLINGS).count()}"
                "SIBLINGS flagged in callset using maximal independent set"
            )
            logger.info(
                f"{related_samples_to_drop_ht.filter(related_samples_to_drop_ht.relationship == PARENT_CHILD).count()}"
                "PARENT_CHILD samples flagged in callset using maximal independent set"
            )
            logger.info(
                f"{related_samples_to_drop_ht.filter(related_samples_to_drop_ht.relationship == DUPLICATE_OR_TWINS).count()}"
                "duplicate samples flagged in callset using maximal independent set"
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(data_source, freeze))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-s",
        "--data_source",
        help="Data source",
        choices=["regeneron", "broad"],
        default="broad",
    )
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )

    parser.add_argument(
        "--skip_pc_relate",
        help="Skip running PC-relate on all samples. \
        NOTE: This needs SSDs on your workers (for the temp files) and no pre-emptibles while the BlockMatrix writes",
        action="store_true",
    )
    parser.add_argument(
        "--min_emission_kinship",
        help="Minimum kinship threshold for emitting a pair of samples in PC relate and filtering related individuals.",
        default=0.05,
        type=float,
    )
    parser.add_argument(
        "--min_individual_maf",
        help="Minor allele frequency cutoff, must be greater that 0.001 \
        because the qc_mt used was already filtered to that maf.",
        default=0.05,
        type=float,
    )
    parser.add_argument(
        "--first_degree_kin_thresholds",
        help="First degree kinship threshold for filtering a pair of samples with a first degree relationship. \
        (Default = (0.1767767, 0.4); 1st degree relatives)",
        default="0.1767767,0.4",
    )
    parser.add_argument(
        "--second_degree_kin_cutoff",
        help="Minimum kinship threshold for filtering a pair of samples with a second degree relationship\
        in PC relate and filtering related individuals. (Default = 0.08838835)",
        default=0.08838835,
        type=float,
    )
    parser.add_argument(
        "--ibd0_0_max",
        help="IBD0 cutoff to determine parent offspring vs full sibling",
        default=0.05,
    )

    parser.add_argument(
        "--skip_filter_dups",
        help="Skip filtering duplicated samples",
        action="store_true",
    )
    parser.add_argument(
        "--skip_infer_families",
        help="Skip extracting duplicate samples and infers families samples based on PC-relate results",
        action="store_true",
    )
    parser.add_argument(
        "--skip_filter_related_samples",
        help="Skip Filter related samples \
        (based on the pairs present from running pc_relate and using the --second_degree_kin_cutoff value for that run)",
        action="store_true",
    )
    parser.add_argument(
        "--dup_kin_cutoff",
        help="Kinship threshold for filtering a pair of duplicate samples. (Default = 0.4)",
        default=0.4,
        type=float,
    )

    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
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
