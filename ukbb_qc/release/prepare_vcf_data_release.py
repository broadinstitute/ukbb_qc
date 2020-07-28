import argparse
import logging
import pickle
import sys
from typing import Dict, List, Union

import hail as hl

from gnomad.resources.grch38.reference_data import dbsnp
from gnomad.sample_qc.ancestry import POP_NAMES, SUBPOPS
from gnomad.utils.generic import get_reference_genome
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import (
    add_as_vcf_info_dict,
    ALLELE_TYPE_FIELDS,
    AS_FIELDS,
    ENTRIES,
    FAF_POPS,
    FORMAT_DICT,
    generic_field_check,
    GROUPS,
    HISTS,
    INFO_DICT,
    make_filters_sanity_check_expr,
    make_hist_bin_edges_expr,
    make_hist_dict,
    make_info_dict,
    make_label_combos,
    make_vcf_filter_dict,
    REGION_FLAG_FIELDS,
    RF_FIELDS,
    SITE_FIELDS,
    set_female_y_metrics_to_na,
    SEXES,
    SORT_ORDER,
    SPARSE_ENTRIES,
    VQSR_FIELDS,
)
from ukbb_qc.assessment.sanity_checks import sample_sum_check, sanity_check_release_mt
from ukbb_qc.resources.basics import (
    get_checkpoint_path,
    logging_path,
    release_header_path,
    release_mt_path,
    release_vcf_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.variant_qc import info_ht_path
from ukbb_qc.slack_creds import slack_token
from ukbb_qc.utils.constants import INTERVAL_QC_PARAMETERS
from ukbb_qc.utils.utils import make_index_dict


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("data_release")
logger.setLevel(logging.INFO)


# Add capture region, interval QC, and sibling singletons to vcf_info_dict
vcf_info_dict["in_capture_region"] = {
    "Description": "Variant falls within an exome capture region"
}
vcf_info_dict["fail_interval_qc"] = {
    "Description": f"Variant falls within a region where less than {INTERVAL_QC_PARAMETERS[0]}% of samples had a mean coverage of {INTERVAL_QC_PARAMETERS[1]}X"
}
vcf_info_dict["sibling_singleton"] = {
    "Description": "Variant was a callset-wide doubleton that was present only within a sibling pair"
}

# Add interval QC, capture region to REGION_FLAG_FIELDS and remove decoy, segdup
INTERVAL_FIELDS = ["fail_interval_qc", "in_capture_region"]
MISSING_REGION_FIELDS = ("decoy", "segdup")
REGION_FLAG_FIELDS = [field for field in REGION_FLAG_FIELDS if field not in MISSING_REGION_FIELDS]
REGION_FLAG_FIELDS.extend(INTERVAL_FIELDS)
# RF_FIELDS.append("interval_qc_pass")

# Add sibling singletons to SITE_FIELDS
SITE_FIELDS.append("sibling_singleton")

# Make subset list (used in properly filling out VCF header descriptions and naming VCF info fields)
SUBSET_LIST = ["", "gnomad_exomes", "gnomad_genomes"]  # empty for ukbb

# Get gnomAD subpop names
GNOMAD_NFE_SUBPOPS = map(lambda x: x.lower(), SUBPOPS["NFE"])
GNOMAD_EAS_SUBPOPS = map(lambda x: x.lower(), SUBPOPS["EAS"])


def make_info_expr(t: Union[hl.MatrixTable, hl.Table]) -> Dict[str, hl.expr.Expression]:
    """
    Makes Hail expression for variant annotations to be included in VCF INFO field.

    :param Table/MatrixTable t: Table/MatrixTable containing variant annotations to be reformatted for VCF export.
    :return: Dictionary containing Hail expressions for relevant INFO annotations.
    :rtype: Dict[str, hl.expr.Expression]
    """
    # Start info dict with region_flag and allele_info fields
    vcf_info_dict = {}
    for field in ALLELE_TYPE_FIELDS:
        vcf_info_dict[field] = t["allele_info"][f"{field}"]
    for field in REGION_FLAG_FIELDS:
        vcf_info_dict[field] = t["region_flag"][f"{field}"]

    # Add site-level annotations to vcf_info_dict
    for field in SITE_FIELDS:
        vcf_info_dict[field] = t["info"][f"{field}"]
    for field in RF_FIELDS:
        vcf_info_dict[field] = t["rf"][f"{field}"]
    for field in VQSR_FIELDS:
        vcf_info_dict[field] = t["vqsr"][f"{field}"]

    # Add AS annotations to info dict
    for field in AS_FIELDS:
        vcf_info_dict[field] = t["info"][f"{field}"]

    for hist in HISTS:
        for prefix in ["adj_qual_hists", "qual_hists"]:
            hist_name = hist
            if "adj" in prefix:
                hist_name = f"{hist}_adj"

            hist_dict = {
                f"{hist_name}_bin_freq": hl.delimit(
                    t[prefix][hist].bin_freq, delimiter="|"
                ),
                f"{hist_name}_bin_edges": hl.delimit(
                    t[prefix][hist].bin_edges, delimiter="|"
                ),
                f"{hist_name}_n_smaller": t[prefix][hist].n_smaller,
                f"{hist_name}_n_larger": t[prefix][hist].n_larger,
            }
            vcf_info_dict.update(hist_dict)
    return vcf_info_dict


def unfurl_nested_annotations(
    t: Union[hl.MatrixTable, hl.Table], gnomad: bool, genome: bool
) -> Dict[str, hl.expr.Expression]:
    """
    Create dictionary keyed by the variant annotation labels to be extracted from variant annotation arrays, where the values
    of the dictionary are Hail Expressions describing how to access the corresponding values.

    :param Table/MatrixTable t: Table/MatrixTable containing the nested variant annotation arrays to be unfurled.
    :param bool gnomad: Whether the annotations are from gnomAD.
    :param bool genome: Whether the annotations are from genome data (relevant only to gnomAD data).
    :return: Dictionary containing variant annotations and their corresponding values.
    :rtype: Dict[str, hl.expr.Expression]
    """
    expr_dict = dict()

    # Set variables to locate necessary fields, compute freq index dicts, and compute faf index dict for UKB
    if gnomad:
        if genome:
            faf = "gnomad_genomes_faf"
            freq = "gnomad_genomes_freq"
            faf_idx = hl.eval(t.gnomad_genomes_faf_index_dict)
            gnomad_prefix = "gnomad_genomes"
            popmax = "gnomad_genomes_popmax"
            freq_idx = make_index_dict(t, "gnomad_genomes_freq_meta")

        else:
            faf = "gnomad_exomes_faf"
            freq = "gnomad_exomes_freq"
            freq_idx = "gnomad_exomes_freq_index_dict"
            faf_idx = hl.eval(t.gnomad_exomes_faf_index_dict)
            gnomad_prefix = "gnomad_exomes"
            popmax = "gnomad_exomes_popmax"
            freq_idx = make_index_dict(t, "gnomad_exomes_freq_meta")

    else:
        faf = "faf"
        freq = "freq"
        faf_idx = make_index_dict(t, "faf_meta")
        popmax = "popmax"
        freq_idx = make_index_dict(t, "freq_meta")

    # Unfurl freq index dict
    for k, i in freq_idx.items():
        entry = k.split("_")

        # skip gnomad subsets
        if entry[0] == "non" or entry[0] == "controls":
            continue

        if gnomad:
            prefix = gnomad_prefix
        else:
            prefix = entry[0]

        if entry == ["raw"] or entry == ["gnomad", "raw"]:
            combo_fields = ["raw"]
        else:
            combo_fields = ["adj"] + entry[1:]

        combo = "_".join(combo_fields)
        combo_dict = {
            f"{prefix}_AC_{combo}": t[freq][i].AC,
            f"{prefix}_AN_{combo}": t[freq][i].AN,
            f"{prefix}_AF_{combo}": t[freq][i].AF,
            f"{prefix}_nhomalt_{combo}": t[freq][i].homozygote_count,
        }
        expr_dict.update(combo_dict)

    ## Unfurl FAF index dict
    for (
        k,
        i,
    ) in faf_idx.items():  # NOTE: faf annotations are all done on adj-only groupings
        entry = k.split("_")

        # skip gnomad subsets
        if entry[0] == "non" or entry[0] == "controls":
            continue

        combo_fields = ["adj"] + entry[1:]
        combo = "_".join(combo_fields)

        if gnomad:
            prefix = gnomad_prefix

            combo_dict = {
                f"{prefix}_faf95_{combo}": hl.or_missing(
                    hl.set(t[faf][i].meta.values()) == set(combo_fields),
                    t[faf][i].faf95,
                ),
                f"{prefix}_faf99_{combo}": hl.or_missing(
                    hl.set(t[faf][i].meta.values()) == set(combo_fields),
                    t[faf][i].faf99,
                ),
            }

        else:
            # NOTE: need to compute UKB separately because UKB no longer has faf meta bundled into faf
            combo_dict = {
                f"faf95_{combo}": hl.or_missing(
                    hl.set(hl.eval(t.faf_meta[i].values())) == set(entry),
                    t[faf][i].faf95,
                ),
                f"faf99_{combo}": hl.or_missing(
                    hl.set(hl.eval(t.faf_meta[i].values())) == set(entry),
                    t[faf][i].faf99,
                ),
            }
        expr_dict.update(combo_dict)

    ## Unfurl popmax
    if gnomad:
        prefix = gnomad_prefix + "_"
    else:
        prefix = ""
    if gnomad and not genome:
        idx = hl.eval(t.globals["gnomad_exomes_popmax_index_dict"]["gnomad"])
        combo_dict = {
            f"{prefix}popmax": t[popmax][idx].pop,
            f"{prefix}AC_popmax": t[popmax][idx].AC,
            f"{prefix}AN_popmax": t[popmax][idx].AN,
            f"{prefix}AF_popmax": t[popmax][idx].AF,
            f"{prefix}nhomalt_popmax": t[popmax][idx].homozygote_count,
        }
    else:
        combo_dict = {
            f"{prefix}popmax": t[popmax].pop,
            f"{prefix}AC_popmax": t[popmax].AC,
            f"{prefix}AN_popmax": t[popmax].AN,
            f"{prefix}AF_popmax": t[popmax].AF,
            f"{prefix}nhomalt_popmax": t[popmax].homozygote_count,
        }
    expr_dict.update(combo_dict)

    # Unfurl UKBB ages
    if not gnomad:
        age_hist_dict = {
            "age_hist_het_bin_freq": hl.delimit(t.age_hist_het.bin_freq, delimiter="|"),
            "age_hist_het_bin_edges": hl.delimit(
                t.age_hist_het.bin_edges, delimiter="|"
            ),
            "age_hist_het_n_smaller": t.age_hist_het.n_smaller,
            "age_hist_het_n_larger": t.age_hist_het.n_larger,
            "age_hist_hom_bin_freq": hl.delimit(t.age_hist_hom.bin_freq, delimiter="|"),
            "age_hist_hom_bin_edges": hl.delimit(
                t.age_hist_hom.bin_edges, delimiter="|"
            ),
            "age_hist_hom_n_smaller": t.age_hist_hom.n_smaller,
            "age_hist_hom_n_larger": t.age_hist_hom.n_larger,
        }
        expr_dict.update(age_hist_dict)

    return expr_dict


def main(args):

    hl.init(log="/release.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)

    try:

        if args.prepare_vcf_mt:
            logger.info("Starting VCF process...")
            logger.info("Getting raw MT and dropping all unnecessary entries...")

            # NOTE: reading in raw MatrixTable to be able to return all samples/variants
            mt = get_ukbb_data(
                data_source,
                freeze,
                key_by_locus_and_alleles=args.key_by_locus_and_alleles,
                split=False,
                raw=True,
                repartition=args.repartition,
                n_partitions=args.raw_partitions,
            ).select_entries(*SPARSE_ENTRIES)

            logger.info("Splitting raw MT...")
            mt = hl.experimental.sparse_split_multi(mt)
            mt = mt.select_entries(*ENTRIES)

            logger.info("Removing low QUAL variants...")
            info_ht = hl.read_table(info_ht_path(data_source, freeze))
            mt = mt.filter_rows(
                (hl.agg.any(mt.GT.is_non_ref() | hl.is_defined(mt.END)))
                & (~info_ht[mt.row_key].AS_lowqual)
            )

            logger.info("Reading in release HT and annotating onto raw MT...")
            ht = hl.read_matrix_table(release_ht_path(*tranche_data))
            mt = mt.annotate_rows(**ht[mt.row_key])
            mt = mt.annotate_globals(**ht.index_globals())

            logger.info("Dropping cohort frequencies (last index of freq)...")
            mt = mt.annotate_rows(freq=mt.freq[-1])
            mt = mt.annotate_globals(freq_meta=mt.freq_meta[-1])

            logger.info("Making histogram bin edges...")
            bin_edges = make_hist_bin_edges_expr(mt.rows(), prefix="")

            logger.info("Getting age hist data...")
            age_hist_data = hl.eval(mt.age_distribution)

            logger.info("Making INFO dict for VCF...")
            vcf_info_dict = INFO_DICT

            # Remove decoy and segdup from info dict
            vcf_info_dict.pop("decoy", None)
            vcf_info_dict.pop("segdup", None)
            vcf_info_dict.update(add_as_vcf_info_dict(vcf_info_dict))

            for subset in SUBSET_LIST:
                vcf_info_dict.update(make_info_dict(subset, dict(group=GROUPS)))
                vcf_info_dict.update(
                    make_info_dict(subset, dict(group=["adj"], sex=SEXES))
                )
                vcf_info_dict.update(
                    make_info_dict(subset, dict(group=["adj"]), faf=True)
                )
                vcf_info_dict.update(
                    make_info_dict(subset, dict(group=["adj"], sex=SEXES), faf=True)
                )
                vcf_info_dict.update(
                    make_info_dict(subset, dict(group=["adj"], pop=POP_NAMES))
                )
                vcf_info_dict.update(
                    make_info_dict(
                        subset, dict(group=["adj"], pop=POP_NAMES, sex=SEXES)
                    )
                )
                vcf_info_dict.update(
                    make_info_dict(subset, dict(group=["adj"], pop=FAF_POPS), faf=True)
                )
                vcf_info_dict.update(
                    make_info_dict(
                        subset, dict(group=["adj"], pop=FAF_POPS, sex=SEXES), faf=True,
                    )
                )

                if "gnomad" in subset:
                    description_text = " in gnomAD"
                    vcf_info_dict.update(
                        make_info_dict(
                            subset, popmax=True, description_text=description_text
                        )
                    )
                    vcf_info_dict.update(
                        make_info_dict(
                            subset,
                            dict(group=["adj"], pop=["nfe"], subpop=GNOMAD_NFE_SUBPOPS),
                            description_text=description_text,
                        )
                    )
                    vcf_info_dict.update(
                        make_info_dict(
                            subset,
                            dict(group=["adj"], pop=["eas"], subpop=GNOMAD_EAS_SUBPOPS),
                            description_text=description_text,
                        )
                    )

                else:
                    vcf_info_dict.update(
                        make_info_dict(
                            subset,
                            bin_edges=bin_edges,
                            popmax=True,
                            age_hist_data="|".join(str(x) for x in age_hist_data),
                        )
                    )
            vcf_info_dict.update(make_hist_dict(bin_edges, adj=False))

            # Adjust keys to remove adj tags before exporting to VCF
            new_vcf_info_dict = {
                i.replace("adj_", "").replace("_adj", "").replace("_adj_", ""): j
                for i, j in vcf_info_dict.items()
            }
            new_vcf_info_dict.update(make_hist_dict(bin_edges, adj=True))

            # Add non-PAR annotation
            mt = mt.annotate_rows(
                nonpar=(mt.locus.in_x_nonpar() | mt.locus.in_y_nonpar())
            )

            logger.info("Constructing INFO field")
            mt = mt.annotate_rows(info=hl.struct(**make_info_expr(mt)))
            mt = mt.annotate_rows(
                info=mt.info.annotate(
                    **unfurl_nested_annotations(mt, gnomad=False, genome=False)
                )
            )
            mt = mt.annotate_rows(
                info=mt.info.annotate(
                    **unfurl_nested_annotations(mt, gnomad=True, genome=True)
                )
            )
            mt = mt.annotate_rows(
                info=mt.info.annotate(
                    **unfurl_nested_annotations(mt, gnomad=True, genome=False)
                )
            )
            mt = mt.annotate_rows(**set_female_y_metrics_to_na(mt))

            # Select relevant fields for VCF export
            mt = mt.select_rows("info", "filters", "rsid", "qual")
            new_vcf_info_dict.update(
                {"vep": {"Description": hl.eval(mt.vep_csq_header)}}
            )
            header_dict = {
                "info": new_vcf_info_dict,
                "filter": make_vcf_filter_dict(
                    mt.rf_globals.snp_cutoff.min_score,
                    mt.rf_globals.indel_cutoff.min_score,
                    mt.rf_globals.inbreeding_cutoff,
                ),
                "format": FORMAT_DICT,
            }

            logger.info("Saving header dict to pickle...")
            with hl.hadoop_open(release_header_path(*tranche_data), "wb") as p:
                pickle.dump(header_dict, p, protocol=pickle.HIGHEST_PROTOCOL)
            mt = mt.annotate_rows(info=mt.info.annotate(vep=mt.vep_csq))
            mt.write(release_mt_path(*tranche_data), args.overwrite)

        if args.sanity_check:
            mt = hl.read_matrix_table(release_mt_path(*tranche_data, temp=True))
            mt = mt.annotate_rows(**set_female_y_metrics_to_na(mt))
            sanity_check_release_mt(
                mt, SUBSET_LIST, missingness_threshold=0.5, verbose=args.verbose
            )

        if args.prepare_release_vcf:

            logger.warning(
                "VCF export will densify! Make sure you have an autoscaling cluster."
            )
            if not args.per_chromosome and not args.parallelize:
                logger.error("Need to choose how to export the release VCF. Exiting...")
                sys.exit(1)

            mt = hl.read_matrix_table(release_mt_path(*tranche_data, temp=True))
            logger.info("Reading header dict from pickle...")
            with hl.hadoop_open(release_header_path(*tranche_data), "rb") as p:
                header_dict = pickle.load(p)

            # Reformat names to remove "adj" pre-export
            # All unlabled frequency information is assumed to be adj
            # All raw frequency information is labeled "_raw"
            row_annots = list(mt.row.info)
            new_row_annots = []
            for x in row_annots:
                if "hist" not in x:
                    x = (
                        x.replace("adj_", "",)
                        .replace("_adj", "")
                        .replace("_adj_", "")
                        .replace("raw_", "")
                    )
                new_row_annots.append(x)

            info_annot_mapping = dict(
                zip(new_row_annots, [mt.info[f"{x}"] for x in row_annots])
            )

            # Confirm all VCF fields/descriptions are present before exporting
            missing_fields = []
            missing_descriptions = []
            for item in ["info", "filter", "format"]:
                if item == "info":
                    annots = new_row_annots
                elif item == "format":
                    annots = list(mt.entry)
                else:
                    annot_mt = mt.explode_rows(mt.filters)
                    annots = list(
                        annot_mt.aggregate_rows(hl.agg.collect_as_set(annot_mt.filters))
                    )

                temp_missing_fields = []
                temp_missing_descriptions = []
                for field in annots:
                    try:
                        description = header_dict[item][field]
                        if len(description) == 0:
                            logger.warning(
                                f"{field} in mt info field has empty description in VCF header!"
                            )
                            temp_missing_descriptions.append(field)
                    except KeyError:
                        logger.warning(
                            f"{field} in mt info field does not exist in VCF header!"
                        )
                        if "hist" not in field:  # NOTE: some hists are not exported
                            temp_missing_fields.append(field)

                missing_fields.extend(temp_missing_fields)
                missing_descriptions.extend(temp_missing_descriptions)

            if len(missing_fields) != 0 or len(missing_descriptions) != 0:
                logger.error(
                    "Some fields are either missing or missing descriptions in the VCF header! Please reconcile."
                )
                return

            mt = mt.transmute_rows(info=hl.struct(**info_annot_mapping))

            # Rearrange INFO field in desired ordering
            drop_hists = (
                [x + "_n_smaller" for x in HISTS]
                + [x + "_bin_edges" for x in HISTS]
                + [x + "_n_larger" for x in HISTS if "dp_" not in x]
            )
            drop_hists.extend(
                [x + "_adj_n_smaller" for x in HISTS]
                + [x + "_adj_bin_edges" for x in HISTS]
                + [x + "_adj_n_larger" for x in HISTS if "dp_" not in x]
                + ["age_hist_hom_bin_edges", "age_hist_het_bin_edges"]
            )

            mt = mt.annotate_rows(
                info=mt.info.select(
                    "AC",
                    "AN",
                    "AF",
                    "rf_tp_probability",
                    *mt.info.drop("AC", "AN", "AF", "rf_tp_probability", *drop_hists),
                )
            )

            # Export VCFs by chromosome
            if args.per_chromosome:
                mt_path = get_checkpoint_path(*tranche_data, name="final_vcf", mt=True)
                mt = mt.write(mt_path)

                logger.info(f"VCF MT count: {mt.count()}")
                rg = get_reference_genome(mt.locus)
                contigs = rg.contigs[:24]  # autosomes + X/Y
                logger.info(f"Contigs: {contigs}")

                for contig in contigs:
                    # Faster way to filter to a contig
                    # TODO: Confirm with hail team if this is the fastest method for this
                    mt = hl.read_matrix_table(mt_path)
                    mt = hl.filter_intervals(mt, [hl.parse_locus_interval(contig)])
                    mt._calculate_new_partitions(10000)
                    intervals = [i for i in intervals if i.start.locus.contig == contig]
                    mt = hl.read_matrix_table(mt_path, _intervals=intervals)

                    logger.info("Densifying and exporting VCF...")
                    mt = hl.experimental.densify(mt)

                    hl.export_vcf(
                        contig_mt,
                        release_vcf_path(*tranche_data, contig=contig),
                        metadata=header_dict,
                    )

            # Export sharded VCF
            if args.parallelize:

                logger.info("Densifying...")
                mt = hl.experimental.densify(mt)
                mt = mt.naive_coalesce(args.n_shards)

                hl.export_vcf(
                    mt,
                    release_vcf_path(*tranche_data),
                    parallel="header_per_shard",
                    metadata=header_dict,
                )
    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(*tranche_data))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int
    )
    parser.add_argument(
        "--prepare_vcf_mt", help="Use release mt to create vcf mt", action="store_true"
    )
    parser.add_argument(
        "--sanity_check", help="Run sanity checks function", action="store_true"
    )
    parser.add_argument(
        "--verbose",
        help="Run sanity checks function with verbose output",
        action="store_true",
    )
    parser.add_argument(
        "--prepare_release_vcf", help="Prepare release VCF", action="store_true"
    )
    export_opts = parser.add_mutually_exclusive_group("Options for how to export VCF")
    export_opts.add_argument(
        "--per_chromosome",
        help="Prepare release VCFs per chromosome",
        action="store_true",
    )
    export_opts.add_argument(
        "--parallelize",
        help="Parallelize VCF export by exporting sharded VCF",
        action="store_true",
    )
    parser.add_argument(
        "--n_shards",
        help="Desired number of shards for output VCF (if --parallelize is set)",
        type=int,
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
