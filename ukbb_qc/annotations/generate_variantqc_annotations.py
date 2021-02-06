import argparse
import logging

import hail as hl

from gnomad.resources.grch38.reference_data import get_truth_ht
from gnomad.sample_qc.relatedness import filter_mt_to_trios
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import VEP_CSQ_HEADER, vep_or_lookup_vep, vep_struct_to_csq
from gnomad.variant_qc.pipeline import generate_sib_stats, generate_trio_stats
from ukbb_qc.resources.basics import get_ukbb_data, logging_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import (
    array_concordance_results_path,
    array_concordance_sites_path,
    inferred_ped_path,
    relatedness_ht_path,
)
from ukbb_qc.resources.variant_qc import (
    get_true_positive_vcf_path,
    var_annotations_ht_path,
)
from ukbb_qc.slack_creds import slack_token

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variantqc_annotations")
logger.setLevel(logging.INFO)


def export_tp_vcf(
    data_source: str,
    freeze: int,
    transmitted_singletons: bool = True,
    sibling_singletons: bool = True,
    array_con_common: bool = True,
):
    """
    Export true positive variants to VCF for use in VQSR
    
    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param transmitted_singletons: Should transmitted singletons be included
    :param sibling_singletons: Should sibling singletons be included
    :param array_con_common: Should common variants that are concordant with array data be included
    :return: None
    """
    if not (transmitted_singletons | sibling_singletons | array_con_common):
        raise ValueError(
            "At least one of transmitted_singletons, sibling_singletons, or array_con_common must be set to True"
        )

    qc_ac_ht = hl.read_table(
        var_annotations_ht_path("allele_counts", data_source, freeze)
    )
    trio_stats_ht = hl.read_table(
        var_annotations_ht_path("trio_stats", data_source, freeze)
    )
    sib_stats_ht = hl.read_table(
        var_annotations_ht_path("sib_stats", data_source, freeze)
    )
    array_con_ht = hl.read_table(
        var_annotations_ht_path("array_exome_concordant_variants", data_source, freeze)
    )
    qc_ac_ht = qc_ac_ht.annotate(
        ukbb_array_con_common=hl.is_defined(array_con_ht[qc_ac_ht.key])
    )

    for transmission_confidence in ["raw", "adj"]:
        filter_expr = False
        true_positive_type = ""
        if transmitted_singletons:
            filter_expr = filter_expr | (
                trio_stats_ht[qc_ac_ht.key][f"n_transmitted_{transmission_confidence}"]
                == 1
            ) & (qc_ac_ht.ac_qc_samples_raw == 2)
            true_positive_type = true_positive_type + "ts_"

        if sibling_singletons:
            filter_expr = filter_expr | (
                sib_stats_ht[qc_ac_ht.key][
                    f"n_sib_shared_variants_{transmission_confidence}"
                ]
                == 1
            ) & (qc_ac_ht.ac_qc_samples_raw == 2)
            true_positive_type = true_positive_type + "ss_"

        if array_con_common:
            filter_expr = filter_expr | qc_ac_ht.ukbb_array_con_common
            true_positive_type = true_positive_type + "ac_"

        ht = qc_ac_ht.filter(filter_expr)
        mt = hl.MatrixTable.from_rows_table(ht)
        logger.info(
            f"Exporting {transmission_confidence} transmitted singleton VCF with {mt.count()} variants..."
        )
        hl.export_vcf(
            mt,
            get_true_positive_vcf_path(
                true_positive_type=true_positive_type,
                adj=(transmission_confidence == "adj"),
                data_source=data_source,
                freeze=freeze,
            ),
            tabix=True,
        )


def main(args):
    hl.init(log="/generate_variantqc_annotations.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze
    n_partitions = args.n_partitions
    overwrite = args.overwrite

    try:
        if args.vep:
            # Need to spin up a cluster with  --requester-pays-allow-buckets gs://hail-us-vep
            # (or --requester-pays-allows-all) in addition to --vep (build)
            logger.info("Running VEP on split hard call MT...")
            mt = get_ukbb_data(data_source, freeze)
            ht = mt.filter_rows(
                (hl.len(mt.alleles) > 1) & hl.agg.any(mt.GT.is_non_ref())
            ).rows()
            ht = vep_or_lookup_vep(ht)
            ht = ht.annotate(vep_csq=vep_struct_to_csq(ht.vep))
            ht = ht.annotate_globals(vep_csq_header=VEP_CSQ_HEADER)
            ht.naive_coalesce(n_partitions).write(
                var_annotations_ht_path("vep", data_source, freeze),
                overwrite=overwrite,
            )

        if args.generate_allele_counts:
            logger.info("Computing allele counts on split hard call MT...")
            mt = get_ukbb_data(data_source, freeze, meta_root="meta")
            mt = mt.filter_rows(
                (hl.len(mt.alleles) > 1) & hl.agg.any(mt.GT.is_non_ref())
            )
            ac_expr = {
                "ac_qc_samples_raw": mt.meta.sample_filters.high_quality
                | mt.meta.sample_filters.control,
                "ac_release_samples_raw": mt.meta.sample_filters.release,
                "ac_qc_samples_unrelated_raw": mt.meta.sample_filters.high_quality
                & ~mt.meta.sample_filters.related,
                "ac_qc_samples_adj": (
                    mt.meta.sample_filters.high_quality | mt.meta.sample_filters.control
                )
                & mt.adj,
                "ac_release_samples_adj": mt.meta.sample_filters.release & mt.adj,
                "ac_qc_samples_unrelated_adj": mt.meta.sample_filters.high_quality
                & ~mt.meta.sample_filters.related
                & mt.adj,
            }
            ht = mt.annotate_rows(
                **{
                    ac_name: hl.agg.filter(expr, hl.agg.sum(mt.GT.n_alt_alleles()))
                    for ac_name, expr in ac_expr.items()
                }
            ).rows()

            ht.naive_coalesce(n_partitions).write(
                var_annotations_ht_path("allele_counts", data_source, freeze),
                overwrite=overwrite,
            )

        if args.generate_trio_stats:
            logger.info("Generate trio statistics on split hard call MT...")
            mt = get_ukbb_data(data_source, freeze)

            ped_fp = inferred_ped_path(data_source, freeze)
            ped = hl.Pedigree.read(ped_fp, delimiter="\t")
            fam_ht = hl.import_fam(ped_fp, delimiter="\t")

            # Filter to autosomes to prevent unnecessary densify of the sex chromosomes
            mt = filter_to_autosomes(mt)
            mt = filter_mt_to_trios(mt, fam_ht)
            mt = hl.experimental.densify(mt)
            mt = mt.filter_rows(hl.len(mt.alleles) == 2)

            mt = hl.trio_matrix(mt, pedigree=ped, complete_trios=True)
            trio_stats_ht = generate_trio_stats(mt, bi_allelic_only=False)
            trio_stats_ht.naive_coalesce(n_partitions).write(
                var_annotations_ht_path("trio_stats", data_source, freeze),
                overwrite=overwrite,
            )

        if args.generate_sibling_stats:
            logger.info(
                "Generating sibling variant sharing statistics on split hard call MT..."
            )
            mt = get_ukbb_data(data_source, freeze)

            relatedness_ht = hl.read_table(relatedness_ht_path(data_source, freeze))
            sib_stats_ht = generate_sib_stats(mt, relatedness_ht, bi_allelic_only=False)
            sib_stats_ht.naive_coalesce(n_partitions).write(
                var_annotations_ht_path("sib_stats", data_source, freeze),
                overwrite=overwrite,
            )

        if args.generate_array_concordant_ht:
            logger.info(
                "Filtering array concordance data based on allele frequency and concordance..."
            )
            variants_ht = hl.read_table(
                array_concordance_results_path(data_source, freeze, sample=False)
            )

            # Get the tranche 2 (200K) allele frequency from the tranche 2 array concordance sites
            sites_ht = hl.read_table(array_concordance_sites_path())
            variants_ht = variants_ht.annotate(AF=sites_ht[variants_ht.key].AF)

            variants_ht = variants_ht.filter(
                (variants_ht.prop_gt_con_non_ref > args.concordance_cutoff)
                & (variants_ht.AF > args.variant_qc_af_cutoff)
            )
            variants_ht = variants_ht.repartition(n_partitions)
            variants_ht = variants_ht.annotate_globals(
                concordance_cutoff=args.concordance_cutoff,
                variant_qc_af_cutoff=args.variant_qc_af_cutoff,
            )
            variants_ht.write(
                var_annotations_ht_path(
                    "array_exome_concordant_variants", data_source, freeze,
                ),
                overwrite=overwrite,
            )

        if args.annotate_truth_data:
            logger.info("Joining truth data annotations...")
            ht = get_ukbb_data(data_source, freeze).rows().select()
            truth_ht = get_truth_ht()
            ht = ht.join(truth_ht, how="left")
            array_con_ht = hl.read_table(
                var_annotations_ht_path(
                    "array_exome_concordant_variants", data_source, freeze
                )
            )

            ht = ht.annotate(ukbb_array_con_common=hl.is_defined(array_con_ht[ht.key]))
            ht = ht.annotate_globals(
                concordance_cutoff=array_con_ht.index_globals().concordance_cutoff,
                variant_qc_af_cutoff=array_con_ht.index_globals().variant_qc_af_cutoff,
            )
            ht = ht.checkpoint(
                var_annotations_ht_path("truth_data", data_source, freeze),
                overwrite=overwrite,
            )
            ht.summarize()

        if args.export_true_positive_vcfs:
            logger.info("Exporting true positive variants to VCFs...")
            export_tp_vcf(
                data_source,
                freeze,
                transmitted_singletons=args.transmitted_singletons,
                sibling_singletons=args.sibling_singletons,
                array_con_common=args.array_con_common,
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(data_source, freeze))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    slack_params = parser.add_argument_group("Slack parameters")
    slack_params.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("-o", "--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument(
        "--n_partitions",
        help="Number of partitions for output HTs. \
        NOTE: This argument will be used for ALL output HTs",
        default=5000,
        type=int,
    )
    parser.add_argument("--vep", help="Runs VEP", action="store_true")
    parser.add_argument(
        "--generate_allele_counts",
        help="Creates a hail Table of allele counts for variants",
        action="store_true",
    )
    parser.add_argument(
        "--generate_trio_stats", help="Calculates trio stats", action="store_true"
    )
    parser.add_argument(
        "--generate_sibling_stats",
        help="Calculated sibling variant sharing stats",
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
    parser.add_argument(
        "--export_true_positive_vcfs",
        help="Exports true positive variants to VCF files.",
        action="store_true",
    )
    parser.add_argument(
        "--transmitted_singletons",
        help="Include transmitted singletons in the exports of true positive variants to VCF files.",
        action="store_true",
    )
    parser.add_argument(
        "--sibling_singletons",
        help="Include sibling singletons in the exports of true positive variants to VCF files.",
        action="store_true",
    )
    training_params.add_argument(
        "--array_con_common",
        help="Include common concordant array variants in the exports of true positive variants to VCF files.",
        action="store_true",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
