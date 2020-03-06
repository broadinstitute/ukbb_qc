import argparse
import hail as hl
import logging
from gnomad.utils.gnomad_functions import adjust_sex_ploidy, annotate_adj
from ukbb_qc.resources.basics import CURRENT_FREEZE, get_ukbb_data, get_ukbb_data_path
from ukbb_qc.resources.variant_qc import var_annotations_ht_path
from ukbb_qc.sample_qc.call_sex import run_impute_ploidy, run_impute_sex


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("hardcalls")
logger.setLevel(logging.INFO)


def main(args):
    hl.init(log="/generate_hardcalls.log", default_reference="GRCh38")

    data_source = args.data_source
    freeze = args.freeze

    if args.impute_sex:
        logger.info("Imputing sex...")
        mt = get_ukbb_data(
            data_source, freeze, split=False, raw=True, key_by_locus_and_alleles=True
        )

        # NOTE: Broad callset does not have filter annotations
        if data_source == "regeneron":
            mt = mt.filter_rows(hl.is_missing(mt.filters))

        sex_ht = default_annotate_sex(
            mt,
            included_intervals=capture_ht_path(data_source, freeze),
            sites_ht=hl.read_table(f_stats_sites_path()),
            aaf_expr="AF",
            gt_expr="LGT",
        )
        sex_ht = sex_ht.checkpoint(sex_ht_path(data_source, freeze), overwrite=True)
    # NOTE: check distributions here before continuing with hardcalls

    if args.write_hardcalls:
        logger.info("Generating split hardcalls...")
        mt = get_ukbb_data(
            data_source, freeze, split=False, raw=True, key_by_locus_and_alleles=True
        )
        sex_ht = hl.read_table(sex_ht_path(data_source, freeze))

        # Add allele data to mt
        allele_data = hl.struct(
            original_alleles=mt.alleles, has_star=hl.any(lambda a: a == "*", mt.alleles)
        )
        mt = mt.annotate_rows(
            allele_data=allele_data.annotate(**add_variant_type(mt.alleles))
        )

        # Split mt, add adj annotation, and adjust sex ploidies
        mt = hl.experimental.sparse_split_multi(mt)
        mt = annotate_adj(mt.select_cols(sex_karyotype=ht[mt.col_key].sex_karyotype))
        mt = mt.select_entries(
            GT=mt.GT, adj=mt.adj
        )  # Note: this is different from gnomAD hardcalls file because no PGT or PID
        mt = adjust_sex_ploidy(mt, mt.sex_karyotype, male_str="XY", female_str="XX")
        mt = mt.select_cols().naive_coalesce(args.n_partitions)

        # Filter out star alleles and checkpoint
        mt = mt.filter_rows(mt.alleles[1] != "*")
        mt = mt.checkpoint(
            get_ukbb_data_path(data_source, freeze, hardcalls=True, split=False),
            args.overwrite,
        )

        # Finish generating allele data
        ht = mt.rows().select("allele_data")
        allele_type = (
            hl.case()
            .when(hl.is_snp(ht.alleles[0], ht.alleles[1]), "snv")
            .when(hl.is_insertion(ht.alleles[0], ht.alleles[1]), "ins")
            .when(hl.is_deletion(ht.alleles[0], ht.alleles[1]), "del")
            .default("complex")
        )
        ht = ht.annotate(
            allele_data=ht.allele_data.annotate(
                allele_type=allele_type,
                was_mixed=ht.allele_data.variant_type == "mixed",
            )
        )
        ht.write(
            var_annotations_ht_path(data_source, freeze, "allele_data"), args.overwrite
        )

    if args.write_nonrefs:
        logger.info("Creating split sparse MT with only non-ref genotypes...")
        mt = get_ukbb_data(
            data_source, freeze, split=False, raw=True, key_by_locus_and_alleles=True
        ).select_cols()
        mt = mt.drop("gvcf_info")
        mt = mt.annotate_entries(is_missing=hl.is_missing(mt.LGT))
        mt = mt.filter_entries(mt.is_missing | mt.LGT.is_non_ref())

        if data_source == "regeneron":
            mt = mt.drop(
                "PL"
            )  # Note: I guess we need to check if PL has issues rather than just always dropping
        mt = hl.experimental.sparse_split_multi(mt)
        mt = mt.filter_entries(mt.is_missing | mt.GT.is_non_ref())
        mt = annotate_adj(mt)
        mt = mt.naive_coalesce(args.n_partitions)
        mt.write(
            get_ukbb_data_path(data_source, freeze, split=True, non_refs_only=True),
            args.overwrite,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
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
        "--n_partitions", help="Desired number of partitions for output", type=int
    )

    parser.add_argument(
        "--impute_sex",
        help="Impute sex on raw MT (prerequisite for creating hardcalls)",
        action="store_true",
    )
    parser.add_argument(
        "--write_hardcalls", help="Creates a split hardcalls mt", action="store_true"
    )
    parser.add_argument(
        "--write_nonrefs",
        help="Creates a split sparse mt with only non-ref genotypes",
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
