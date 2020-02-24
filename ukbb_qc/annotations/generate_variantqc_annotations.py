from gnomad_hail import *
from gnomad_qc.annotations.generate_qc_annotations import *
from gnomad_hail.utils.relatedness import get_relationship_expr
from gnomad_hail.utils.annotations import (
    default_generate_sib_stats,
    generate_family_stats,
)
import gnomad_hail.resources.basics as gres
from gnomad_hail.utils.generic import vep_or_lookup_vep, vep_struct_to_csq
from ukbb_qc.resources.resources import *
from ukbb_qc.utils.utils import *
import argparse

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variantqc_annotations")
logger.setLevel(logging.INFO)


def annotate_truth_data(ht: hl.Table, truth_tables: Dict[str, hl.Table]) -> hl.Table:
    """
    Writes bi-allelic sites MT with the following annotations by default:
     - Omni
     - HapMap
     - 1KG high conf SNVs
     - Mills

    Additional truthset tables can be passed with `truth_tables`

    :param ht: Table to annotate
    :param truth_tables: A dictionary containing additional truth set hail tables keyed by the annotation to use in HT
    :return: Table with qc annotations
    """
    truth_ht = hl.read_table(truth_ht_path)
    ht = ht.join(truth_ht, how="left")

    ht = ht.annotate(
        **{root: hl.is_defined(truth[ht.key]) for root, truth in truth_tables.items()}
    )

    return ht


def main(args):
    hl.init(log="/generate_variantqc_annotations.log")
    data_source = args.data_source
    freeze = args.freeze

    if args.vep:
        logger.info(f"Running VEP on the MT...")
        ht = get_ukbb_data(data_source, freeze).rows().select()
        ht = vep_or_lookup_vep(ht, reference="GRCh38")
        ht = ht.annotate(vep_csq=vep_struct_to_csq(ht.vep))
        ht.write(var_annotations_ht_path(data_source, freeze, "vep"), args.overwrite)

    if (
        args.generate_allele_data
    ):  # TODO: I think this was moved to another location, if so I will remove from here
        mt = get_ukbb_data(data_source, freeze, split=False)
        generate_allele_data(mt).write(
            var_annotations_ht_path(data_source, freeze, "allele_data"),
            overwrite=args.overwrite,
        )

    if args.generate_qc_annotations:
        # Turn on spark speculation: --properties 'spark:spark.speculation=true,spark:spark.speculation.quantile=0.9,spark:spark.speculation.multiplier=3'
        mt = get_ukbb_data(data_source, freeze, non_refs_only=True, meta_root="meta")
        mt = generate_qc_annotations(
            mt,
            all_annotations=args.calculate_all_annotations,
            medians=args.calculate_medians,
        )
        mt.write(
            var_annotations_ht_path(data_source, freeze, "qc_stats"),
            stage_locally=True,
            overwrite=args.overwrite,
        )

    if args.generate_call_stats:  # TODO: This likely needs to be moved after densify
        mt = get_ukbb_data(data_source, freeze, meta_root="meta")
        ht = mt.annotate_rows(
            qc_callstats=hl.stuct(
                **{
                    group: hl.agg.filter(
                        mt[group], hl.agg.call_stats(mt.GT, mt.alleles)
                    )
                    for group in {
                        "qc_samples_raw": mt.meta.high_quality,
                        "all_samples_raw": True,
                    }
                }
            )
        ).rows()

        ht.write(
            var_annotations_ht_path(data_source, freeze, "call_stats"),
            overwrite=args.overwrite,
        )

    if args.generate_family_stats:  # TODO: Need to move this to be after a densify step
        mt = get_ukbb_data(data_source, freeze, meta_root="meta")
        ht, sample_table = generate_family_stats(
            mt, inferred_ped_path(data_source, freeze), args.include_adj_family_stats
        )
        ht.write(
            var_annotations_ht_path(data_source, freeze, "family_stats"),
            stage_locally=True,
            overwrite=args.overwrite,
        )
        sample_table.write(
            sample_annotations_table_path(data_source, freeze, "family_stats"),
            stage_locally=True,
            overwrite=args.overwrite,
        )

    if args.generate_sibling_singletons:
        relatedness_ht = hl.read_table(relatedness_ht_path(data_source, freeze))
        sex_ht = hl.read_table(sex_ht_path(data_source, freeze))
        mt = get_ukbb_data(data_source, freeze)

        relatedness_ht = relatedness_ht.annotate(
            relationship=get_relationship_expr(
                kin_expr=relatedness_ht.kin,
                ibd0_expr=relatedness_ht.ibd0,
                ibd1_expr=relatedness_ht.ibd1,
                ibd2_expr=relatedness_ht.ibd2,
            )
        )

        sib_stats_ht = default_generate_sib_stats(mt, relatedness_ht, sex_ht)

        sib_stats_ht.naive_coalesce(500).write(
            var_annotations_ht_path(data_source, freeze, "sib_stats"),
            overwrite=args.overwrite,
        )

    if args.generate_array_concordant_ht:
        variants = hl.read_table(array_variant_concordance_path(data_source, freeze))
        callrate_cutoff = variants.callrate_cutoff.take(1)[0]
        af_cutoff = variants.af_cutoff.take(1)[0]

        exome_mt = hl.read_matrix_table(
            get_mt_checkpoint_path(
                data_source,
                freeze,
                name=f"exome_subset_concordance_callrate_{callrate_cutoff}_af_{af_cutoff}",
            )
        )

        variants = variants.annotate(AF=exome_mt.rows()[variants.key].variant_qc.AF[1])
        variants = variants.filter(
            (variants.prop_gt_con_non_ref > args.concordance_cutoff)
            & (variants.AF > args.variant_qc_af_cutoff)
        )
        variants = variants.repartition(1000)
        variants.write(
            var_annotations_ht_path(
                data_source,
                freeze,
                f"array_con_con_{args.concordance_cutoff}_AF_{args.variant_qc_af_cutoff}",
            ),
            overwrite=args.overwrite,
        )

    if args.annotate_truth_data:
        ht = get_ukbb_data(data_source, freeze, meta_root=None).select_rows()
        truth_ht = annotate_truth_data(
            ht,
            {
                "sib_singletons": hl.read_table(
                    var_annotations_ht_path(data_source, freeze, "sib_singletons.train")
                ),
                "ukbb_array_con_common": hl.read_table(
                    var_annotations_ht_path(
                        data_source,
                        freeze,
                        "array_variant_concordance_callrate_0.95_non_ref_con_0.9_AF_0.001",
                    )
                ),
            },
        )

        truth_ht = truth_ht.checkpoint(
            var_annotations_ht_path(data_source, freeze, "truth_data"),
            overwrite=args.overwrite,
        )
        truth_ht.summarize()


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
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("-o", "--overwrite", help="Overwrite data", action="store_true")

    parser.add_argument("--vep", help="Runs VEP", action="store_true")
    parser.add_argument(
        "--generate_allele_data", help="Calculates allele data", action="store_true"
    )
    parser.add_argument(
        "--generate_qc_annotations",
        help="Calculates QC annotations",
        action="store_true",
    )
    parser.add_argument(
        "--generate_call_stats", help="Calculates call stats", action="store_true"
    )
    parser.add_argument(
        "--generate_family_stats", help="Calculates family stats", action="store_true"
    )
    parser.add_argument(
        "--generate_sibling_singletons",
        help="Creates a hail Table of variants that are sibling singletons",
        action="store_true",
    )
    parser.add_argument(
        "--generate_array_concordant_ht",
        help="Creates a hail Table of array concordant variants",
        action="store_true",
    )
    parser.add_argument(
        "--concordance_cutoff",
        help="Array exome concordance cutoff for variant QC HT.",
        type=float,
        default=0.9,
    )
    parser.add_argument(
        "--variant_qc_af_cutoff",
        help="Allele frequency cutoff used for variant QC HT creation, must be equal to or greater than af_cutoff.",
        type=float,
        default=0.001,
    )
    parser.add_argument(
        "--num_var_per_sibs_cutoff",
        help="Max number of sibling singletons in a pair for the pair to be included in truth set",
        default=40,
    )
    parser.add_argument(
        "--include_adj_family_stats",
        help="Also calculate family stats for adj genotypes",
        action="store_true",
    )
    parser.add_argument(
        "--test_train_split",
        help="Percentage of truth data to hold back for testing",
        default=0.2,
    )
    parser.add_argument(
        "--annotate_truth_data",
        help="Creates a HT of UKBB variants annotated with truth sites",
        action="store_true",
    )
    parser.add_argument(
        "--calculate_medians",
        help="Calculate metric medians (warning: slow)",
        action="store_true",
    )
    parser.add_argument(
        "--calculate_all_annotations",
        help="Calculation many more annotations (warning: slow)",
        action="store_true",
    )
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
