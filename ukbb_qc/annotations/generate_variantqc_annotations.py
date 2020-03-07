import argparse
import logging
from typing import Dict
import hail as hl
from gnomad.resources.grch38.reference_data import get_truth_ht
from gnomad.utils.annotations import (
    default_generate_sib_stats,
    filter_mt_to_trios,
    default_generate_trio_stats,
)
from gnomad.utils.generic import vep_or_lookup_vep, vep_struct_to_csq
from gnomad.utils.slack import try_slack
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.basics import get_ukbb_data
from ukbb_qc.resources.sample_qc import (
    array_variant_concordance_path,
    related_drop_path,
    sex_ht_path,
    inferred_ped_path,
)
from ukbb_qc.resources.variant_qc import var_annotations_ht_path, info_ht_path

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
    truth_ht = get_truth_ht()
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

    if args.generate_trio_stats:
        mt = get_ukbb_data(data_source, freeze, meta_root="meta")
        mt = mt.annotate_cols(is_female=mt.meta.is_female)

        ped_fp = inferred_ped_path(data_source, freeze)
        ped = hl.Pedigree.read(ped_fp, delimiter="\t")
        fam_ht = hl.import_fam(ped_fp, delimiter="\t")

        mt = filter_mt_to_trios(mt, fam_ht)
        mt = mt.select_entries("GT", "GQ", "AD", "END", "adj")

        mt = hl.experimental.densify(mt)

        mt = hl.trio_matrix(mt, pedigree=ped, complete_trios=True)
        trio_stats_ht = default_generate_trio_stats(mt)
        trio_stats_ht.write(
            var_annotations_ht_path(data_source, freeze, "trio_stats"),
            overwrite=args.overwrite,
        )

    if args.generate_sibling_stats:
        mt = get_ukbb_data(data_source, freeze, meta_root="meta")
        relatedness_ht = hl.read_table(related_drop_path(data_source, freeze))
        sex_ht = hl.read_table(sex_ht_path(data_source, freeze))

        sib_stats_ht = default_generate_sib_stats(mt, relatedness_ht, sex_ht)

        sib_stats_ht.naive_coalesce(500).write(
            var_annotations_ht_path(data_source, freeze, "sib_stats"),
            overwrite=args.overwrite,
        )

    if args.generate_array_concordant_ht:
        variants = hl.read_table(array_variant_concordance_path(data_source, freeze))
        info_ht = hl.read_table(info_ht_path(data_source, freeze, split=True))

        variants = variants.annotate(AF=info_ht[variants.key].ac_qc_samples_raw)
        variants = variants.filter(
            (variants.prop_gt_con_non_ref > args.concordance_cutoff)
            & (variants.AF > args.variant_qc_af_cutoff)
        )
        variants = variants.repartition(1000)
        variants = variants.annotate_globals(
            concordance_cutoff=args.concordance_cutoff,
            variant_qc_af_cutoff=args.variant_qc_af_cutoff,
        )
        variants.write(
            var_annotations_ht_path(
                data_source, freeze, "array_exome_concordant_variants",
            ),
            overwrite=args.overwrite,
        )

    if args.annotate_truth_data:
        ht = get_ukbb_data(data_source, freeze, meta_root=None).rows()
        truth_ht = annotate_truth_data(
            ht,
            {
                "ukbb_array_con_common": hl.read_table(
                    var_annotations_ht_path(
                        data_source,
                        freeze,
                        var_annotations_ht_path(
                            data_source, freeze, "array_exome_concordant_variants",
                        ),
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
        "--generate_trio_stats", help="Calculates trio stats", action="store_true"
    )
    parser.add_argument(
        "--generate_sibling_stats",
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
        "--annotate_truth_data",
        help="Creates a HT of UKBB variants annotated with truth sites",
        action="store_true",
    )
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
