import argparse
import logging
import hail as hl
from gnomad.utils.slack import try_slack
from gnomad.utils.gnomad_functions import (
    adjust_sex_ploidy,
    annotate_adj,
    add_variant_type,
)
from gnomad.utils.sample_qc import default_annotate_sex
from ukbb_qc.resources.basics import (
    get_ukbb_data,
    get_ukbb_data_path,
    capture_ht_path,
    logging_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import f_stat_sites_path, sex_ht_path
from ukbb_qc.resources.variant_qc import var_annotations_ht_path


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("hardcalls")
logger.setLevel(logging.INFO)


def main(args):

    hl.init(log="/generate_hardcalls.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze

    try:
        if args.impute_sex:
            logger.info("Imputing sex...")
            # NOTE: Remove key_by_locus_and_alleles when MTs are keyed by locus and alleles by default
            # currently only keyed by locus
            mt = get_ukbb_data(
                data_source,
                freeze,
                key_by_locus_and_alleles=True,
                split=False,
                raw=True,
                repartition=args.repartition,
                n_partitions=args.n_partitions,
            )

            sex_ht = default_annotate_sex(
                mt,
                included_intervals=hl.read_table(capture_ht_path(data_source)),
                sites_ht=hl.read_table(f_stat_sites_path()),
                aaf_expr="AF",
                gt_expr="LGT",
            )
            sex_ht.write(sex_ht_path(data_source, freeze), overwrite=True)
        # NOTE: check distributions here before continuing with hardcalls

        if args.write_hardcalls:
            logger.info("Generating split hardcalls...")
            mt = get_ukbb_data(
                data_source,
                freeze,
                key_by_locus_and_alleles=True,
                split=False,
                raw=True,
                repartition=args.repartition,
                n_partitions=args.n_partitions,
            )
            sex_ht = hl.read_table(sex_ht_path(data_source, freeze))

            # Add allele data to mt
            allele_data = hl.struct(
                original_alleles=mt.alleles,
                has_star=hl.any(lambda a: a == "*", mt.alleles),
            )
            mt = mt.annotate_rows(
                allele_data=allele_data.annotate(**add_variant_type(mt.alleles))
            )

            # Split mt, add adj annotation, and adjust sex ploidies
            mt = hl.experimental.sparse_split_multi(mt)
            mt = annotate_adj(
                mt.select_cols(sex_karyotype=sex_ht[mt.col_key].sex_karyotype),
                haploid_adj_dp=5,
            )
            mt = mt.select_entries(
                GT=mt.GT, adj=mt.adj
            )  # Note: this is different from gnomAD hardcalls file because no PGT or PID
            mt = adjust_sex_ploidy(mt, mt.sex_karyotype, male_str="XY", female_str="XX")
            mt = mt.select_cols().naive_coalesce(args.n_partitions)

            # Filter out star alleles and checkpoint
            mt = mt.filter_rows(mt.alleles[1] != "*")
            mt = mt.checkpoint(
                get_ukbb_data_path(data_source, freeze, hardcalls=True), args.overwrite,
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
                var_annotations_ht_path(data_source, freeze, "allele_data"),
                args.overwrite,
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(data_source, freeze))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to.",
    )
    parser.add_argument(
        "--repartition",
        help="Repartition raw MT on read. Needs to be true for tranche 3/freeze 6/300K",
        action="store_true",
    )
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int,
    )
    parser.add_argument(
        "--n_partitions",
        help="Desired number of partitions for output. Also used to repartition raw MT on read (necessary only for tranche 3/freeze 6/300K!",
        type=int,
    )
    parser.add_argument(
        "--impute_sex",
        help="Impute sex on raw MT (prerequisite for creating hardcalls)",
        action="store_true",
    )
    parser.add_argument(
        "--write_hardcalls", help="Creates a split hardcalls mt", action="store_true",
    )
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
