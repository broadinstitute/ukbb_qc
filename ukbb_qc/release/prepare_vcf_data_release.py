import argparse
import json
import logging
import pickle
import sys
from typing import Dict, List, Union

import hail as hl

from gnomad.resources.grch37.gnomad import SUBPOPS
from gnomad.sample_qc.ancestry import POP_NAMES
from gnomad.utils.generic import get_reference_genome
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import (
    add_as_vcf_info_dict,
    ALLELE_TYPE_FIELDS,
    AS_FIELDS,
    ENTRIES,
    FAF_POPS,
    FORMAT_DICT,
    GROUPS,
    HISTS,
    INFO_DICT,
    make_hist_bin_edges_expr,
    make_hist_dict,
    make_info_dict,
    make_vcf_filter_dict,
    REGION_FLAG_FIELDS,
    RF_FIELDS,
    SITE_FIELDS,
    set_female_y_metrics_to_na,
    SEXES,
    SPARSE_ENTRIES,
    VQSR_FIELDS,
)
from ukbb_qc.assessment.sanity_checks import (
    sanity_check_release_mt,
    vcf_field_check,
)
from ukbb_qc.resources.basics import (
    get_checkpoint_path,
    get_ukbb_data,
    logging_path,
    release_header_path,
    release_ht_path,
    release_mt_path,
    release_vcf_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.variant_qc import info_ht_path
from ukbb_qc.slack_creds import slack_token
from ukbb_qc.utils.utils import make_index_dict


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("vcf_release")
logger.setLevel(logging.INFO)


# Add capture region, interval QC, and sibling singletons to vcf_info_dict
VCF_INFO_DICT = INFO_DICT
VCF_INFO_DICT["in_capture_region"] = {
    "Description": "Variant falls within an exome capture region"
}
VCF_INFO_DICT["sibling_singleton"] = {
    "Description": "Variant was a callset-wide doubleton that was present only within a sibling pair"
}

# Add interval QC, capture region to REGION_FLAG_FIELDS and remove decoy, segdup
# NOTE: MISSING_REGION_FIELDS could change for 500K if we get hg38 files
INTERVAL_FIELDS = ["fail_interval_qc", "in_capture_region"]
MISSING_REGION_FIELDS = ("decoy", "segdup")
REGION_FLAG_FIELDS = [
    field for field in REGION_FLAG_FIELDS if field not in MISSING_REGION_FIELDS
]
REGION_FLAG_FIELDS.extend(INTERVAL_FIELDS)

# Add sibling singletons to AS_FIELDS
AS_FIELDS.append("sibling_singleton")

# Make subset list (used in properly filling out VCF header descriptions and naming VCF info fields)
SUBSET_LIST = ["", "gnomad_exomes", "gnomad_genomes"]  # empty for ukbb

# Get gnomAD subpop names
GNOMAD_NFE_SUBPOPS = map(lambda x: x.lower(), SUBPOPS["NFE"])
GNOMAD_EAS_SUBPOPS = map(lambda x: x.lower(), SUBPOPS["EAS"])


def populate_info_dict(
    subpops: Dict[str, List[str]],
    bin_edges: Dict[str, str],
    age_hist_data: str,
    info_dict: Dict[str, Dict[str, str]] = VCF_INFO_DICT,
    subset_list: List[str] = SUBSET_LIST,
    groups: List[str] = GROUPS,
    sexes: List[str] = SEXES,
    pops: List[str] = POP_NAMES,
    faf_pops: List[str] = FAF_POPS,
    gnomad_nfe_subpops: List[str] = GNOMAD_NFE_SUBPOPS,
    gnomad_eas_subpops: List[str] = GNOMAD_EAS_SUBPOPS,
) -> Dict[str, Dict[str, str]]:
    """
    Calls `make_info_dict` and `make_hist_dict` to populate INFO dictionary with specific sexes, population names, and filtering allele frequency (faf) pops.

    Used during VCF export.

    Creates:
        - INFO fields for age histograms (bin freq, n_smaller, and n_larger for heterozygous and homozygous variant carriers)
        - INFO fields for popmax AC, AN, AF, nhomalt, and popmax population
        - INFO fields for AC, AN, AF, nhomalt for each combination of sample population, sex, and subpopulation, both for adj and raw data
        - INFO fields for filtering allele frequency (faf) annotations 
        - INFO fields for variant histograms (hist_bin_freq, hist_n_smaller, hist_n_larger for each histogram)

    :param Dict[str, List[str]] subpops: Dictionary of global population names (keys)
        and all hybrid population cluster names associated with that global pop (values). 
    :param Dict[str, str] bin_edges: Dictionary of variant annotation histograms and their associated bin edges.
    :param str age_hist_data: Pipe-delimited string of age histograms, from `get_age_distributions`.
    :param Dict[str, Dict[str, str]] info_dict: INFO dict to be populated.
    :param List[str] subset_list: List of sample subsets in dataset. Default is SUBSET_LIST.
    :param List[str] groups: List of sample groups [adj, raw]. Default is GROUPS.
    :param List[str] sexes: List of sample sexes. Default is SEXES.
    :param List[str] pop_names: List of sample global population names. Default is POP_NAMES.
    :param List[str] faf_pops: List of faf population names. Default is FAF_POPS.
    :param List[str] gnomad_nfe_subpops: List of nfe subpopulations in gnomAD. Default is GNOMAD_NFE_SUBPOPS.
    :param List[str] gnomad_eas_subpops: List of eas subpopulations in gnomAD. Default is GNOMAD_EAS_SUBPOPS.
    :rtype: Dict[str, Dict[str, str]]
    """
    vcf_info_dict = info_dict

    # Remove MISSING_REGION_FIELDS from info dict
    for field in MISSING_REGION_FIELDS:
        vcf_info_dict.pop(field, None)
    vcf_info_dict.update(add_as_vcf_info_dict(vcf_info_dict))

    all_label_groups = [
        dict(group=["adj"], sex=sexes),
        dict(group=["adj"], pop=pops),
        dict(group=["adj"], pop=pops, sex=sexes),
    ]
    faf_label_groups = [
        dict(group=["adj"]),
        dict(group=["adj"], sex=sexes),
        dict(group=["adj"], pop=faf_pops),
        dict(group=["adj"], pop=faf_pops, sex=sexes),
    ]
    for subset in subset_list:
        vcf_info_dict.update(make_info_dict(subset, dict(group=groups)))

        for label_group in all_label_groups:
            vcf_info_dict.update(make_info_dict(subset, label_group))
        for label_group in faf_label_groups:
            vcf_info_dict.update(make_info_dict(subset, label_group, faf=True))

        if "gnomad" in subset:
            description_text = " in gnomAD"
            vcf_info_dict.update(
                make_info_dict(subset, popmax=True, description_text=description_text)
            )
            vcf_info_dict.update(
                make_info_dict(
                    subset,
                    dict(group=["adj"], pop=["nfe"], subpop=gnomad_nfe_subpops),
                    description_text=description_text,
                )
            )
            vcf_info_dict.update(
                make_info_dict(
                    subset,
                    dict(group=["adj"], pop=["eas"], subpop=gnomad_eas_subpops),
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

            for pop in subpops:
                vcf_info_dict.update(
                    make_info_dict(
                        subset,
                        dict(group=["adj"], pop=[pop], subpop=subpops[pop]),
                        description_text=description_text,
                    )
                )

    # Add variant quality histograms to info dict
    vcf_info_dict.update(make_hist_dict(bin_edges, adj=True))
    vcf_info_dict.update(make_hist_dict(bin_edges, adj=False))
    return vcf_info_dict


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
        for prefix in ["qual_hists", "raw_qual_hists"]:
            hist_name = hist
            if "raw" in prefix:
                hist_name = f"{hist}_raw"

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
    t: Union[hl.MatrixTable, hl.Table], gnomad: bool, genome: bool, subpops: List[str]
) -> Dict[str, hl.expr.Expression]:
    """
    Create dictionary keyed by the variant annotation labels to be extracted from variant annotation arrays, where the values
    of the dictionary are Hail Expressions describing how to access the corresponding values.

    :param Table/MatrixTable t: Table/MatrixTable containing the nested variant annotation arrays to be unfurled.
    :param bool gnomad: Whether the annotations are from gnomAD.
    :param bool genome: Whether the annotations are from genome data (relevant only to gnomAD data).
    :param List[str] subpops: List of all UKBB subpops (possible hybrid population cluster names). 
    :return: Dictionary containing variant annotations and their corresponding values.
    :rtype: Dict[str, hl.expr.Expression]
    """
    expr_dict = dict()

    # Set variables to locate necessary fields, compute freq index dicts, and compute faf index dict for UKBB
    if gnomad:
        data_type = "genomes" if genome else "exomes"
        gnomad_prefix = f"gnomad_{data_type}"
        popmax = f"{gnomad_prefix}_popmax"
        faf = f"{gnomad_prefix}_faf"
        freq = f"{gnomad_prefix}_freq"
        faf_idx = hl.eval(t.globals[f"{gnomad_prefix}_faf_index_dict"])
        freq_idx = make_index_dict(
            t, f"{gnomad_prefix}_freq_meta", [GNOMAD_NFE_SUBPOPS + GNOMAD_EAS_SUBPOPS]
        )

    else:
        faf = "faf"
        freq = "freq"
        faf_idx = make_index_dict(t, "faf_meta", subpops)
        popmax = "popmax"
        freq_idx = make_index_dict(t, "freq_meta", subpops)

    # Unfurl freq index dict
    # Cycles through each key and index (e.g., k=adj_afr, i=31)
    for k, i in freq_idx.items():
        prefix = ""
        if gnomad:
            prefix = f"{gnomad_prefix}_"

        # Split key to split group (adj, raw) and pop entry
        # e.g., split 'adj_afr' into ['adj', 'afr']
        entry = k.split("_")

        # Add 'adj' and 'raw' tags to combo_fields
        if entry == ["raw"]:
            combo_fields = entry
        else:
            # Add 'adj' tag AFTER all other entries
            # e.g., if entry is ['adj', 'afr'], set combo fields to ['afr', 'adj']
            combo_fields = entry[1:] + [entry[0]]

        combo = "_".join(combo_fields)
        combo_dict = {
            f"{prefix}AC_{combo}": t[freq][i].AC,
            f"{prefix}AN_{combo}": t[freq][i].AN,
            f"{prefix}AF_{combo}": t[freq][i].AF,
            f"{prefix}nhomalt_{combo}": t[freq][i].homozygote_count,
        }
        expr_dict.update(combo_dict)

    ## Unfurl FAF index dict
    for (
        k,
        i,
    ) in faf_idx.items():  # NOTE: faf annotations are all done on adj-only groupings
        entry = k.split("_")

        # Create combo_fields in same way as above ([pop, adj])
        combo_fields = entry[1:] + entry[0]
        combo = "_".join(combo_fields)

        if gnomad:
            prefix = f"{gnomad_prefix}_"
            combo_dict = {
                f"{prefix}faf95_{combo}": hl.or_missing(
                    hl.set(t[faf][i].meta.values()) == set(combo_fields),
                    t[faf][i].faf95,
                ),
                f"{prefix}faf99_{combo}": hl.or_missing(
                    hl.set(t[faf][i].meta.values()) == set(combo_fields),
                    t[faf][i].faf99,
                ),
            }

        else:
            # NOTE: need to compute UKBB separately because UKBB no longer has faf meta bundled into faf
            combo_dict = {
                f"faf95_{combo}": hl.or_missing(
                    hl.set(hl.eval(t.faf_meta[i].values())) == set(combo_fields),
                    t[faf][i].faf95,
                ),
                f"faf99_{combo}": hl.or_missing(
                    hl.set(hl.eval(t.faf_meta[i].values())) == set(combo_fields),
                    t[faf][i].faf99,
                ),
            }
        expr_dict.update(combo_dict)

    # Unfurl popmax
    if gnomad:
        prefix = f"{gnomad_prefix}_"
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

    hl.init(log="/vcf_release.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)
    hybrid_pop_map = args.hybrid_pop_map
    hybrid_pops = [pop for sublist in list(hybrid_pop_map.values()) for pop in sublist]

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

            logger.info("Reading in release HT and annotating onto raw MT...")
            ht = hl.read_matrix_table(release_ht_path(*tranche_data))
            mt = mt.annotate_rows(**ht[mt.row_key])
            mt = mt.annotate_globals(**ht.index_globals())

            logger.info(
                "Dropping cohort frequencies (necessary only for internal use; at last index of freq struct)..."
            )
            mt = mt.annotate_rows(freq=mt.freq[-1])
            mt = mt.annotate_globals(freq_meta=mt.freq_meta[-1])

            logger.info("Making histogram bin edges...")
            bin_edges = make_hist_bin_edges_expr(mt.rows(), prefix="")

            logger.info("Getting age hist data...")
            age_hist_data = hl.eval(mt.age_distribution)

            logger.info("Making INFO dict for VCF...")
            vcf_info_dict = populate_info_dict(
                subpops=hybrid_pop_map,
                bin_edges=bin_edges,
                age_hist_data=age_hist_data,
            )

            # Add interval QC parameters to INFO dict
            pct_samples = hl.eval(mt.rf_globals.interval_qc_cutoffs.pct_samples)
            autosome_cov = hl.eval(mt.rf_globals.interval_qc_cutoffs.autosome_cov)
            allosome_cov = hl.eval(mt.rf_globals.interval_qc_cutoffs.xy_cov)
            vcf_info_dict["fail_interval_qc"] = {
                "Description": f"Variant falls within a region where less than {pct_samples}% \
                                 of samples had a mean coverage of {autosome_cov}X on autosomes and \
                                 {allosome_cov}X on sex chromosomes"
            }

            # Adjust keys to remove adj tags before exporting to VCF
            new_vcf_info_dict = {
                i.replace("_adj", ""): j for i, j in vcf_info_dict.items()
            }

            # Add non-PAR annotation
            mt = mt.annotate_rows(
                nonpar=(mt.locus.in_x_nonpar() | mt.locus.in_y_nonpar())
            )

            logger.info("Constructing INFO field")
            # Add variant annotations to INFO field
            # This adds annotations from:
            #   RF struct, VQSR struct, allele_info struct,
            #   info struct (site and allele-specific annotations),
            #   region_flag struct, and
            #   raw_qual_hists/qual_hists structs.
            mt = mt.annotate_rows(info=hl.struct(**make_info_expr(mt)))
            # Unfurl nested UKBB frequency annotations and add to INFO field
            mt = mt.annotate_rows(
                info=mt.info.annotate(
                    **unfurl_nested_annotations(
                        mt, gnomad=False, genome=False, subpops=hybrid_pops
                    )
                )
            )
            # Unfurl nested gnomAD genome frequency annotations and add to info field
            mt = mt.annotate_rows(
                info=mt.info.annotate(
                    **unfurl_nested_annotations(
                        mt, gnomad=True, genome=True, subpops=hybrid_pops
                    )
                )
            )
            # Unfurl nested gnomAD exome frequency annotations and add to info field
            mt = mt.annotate_rows(
                info=mt.info.annotate(
                    **unfurl_nested_annotations(
                        mt, gnomad=True, genome=False, subpops=hybrid_pops
                    )
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
            mt = hl.read_matrix_table(release_mt_path(*tranche_data))
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

            mt = hl.read_matrix_table(release_mt_path(*tranche_data))
            logger.info("Reading header dict from pickle...")
            with hl.hadoop_open(release_header_path(*tranche_data), "rb") as p:
                header_dict = pickle.load(p)

            # Reformat names to remove "adj" pre-export
            # All unlabeled frequency information is assumed to be adj
            row_annots = list(mt.row.info)
            new_row_annots = []
            for x in row_annots:
                x = x.replace("_adj", "")
                new_row_annots.append(x)

            info_annot_mapping = dict(
                zip(new_row_annots, [mt.info[f"{x}"] for x in row_annots])
            )

            # Confirm all VCF fields and descriptions are present
            if not vcf_field_check(mt, header_dict, new_row_annots, list(mt.entry)):
                logger.error("Did not pass VCF field check.")
                return

            mt = mt.transmute_rows(info=hl.struct(**info_annot_mapping))

            # Rearrange INFO field in desired ordering
            # TODO: figure out which hists we want to export and only create those for 500K
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
                ht = mt.rows().checkpoint(
                    get_checkpoint_path(*tranche_data, name="flat_vcf_ready", mt=False),
                    overwrite=args.overwrite,
                )

                rg = get_reference_genome(mt.locus)
                contigs = rg.contigs[:24]  # autosomes + X/Y
                logger.info(f"Contigs: {contigs}")

                for contig in contigs:
                    # Faster way to filter to a contig
                    # TODO: Confirm with hail team if this is the fastest method for this
                    mt = hl.read_matrix_table(release_mt_path(*tranche_data))
                    mt = hl.filter_intervals(mt, [hl.parse_locus_interval(contig)])
                    intervals = mt._calculate_new_partitions(10000)
                    mt = hl.read_matrix_table(
                        release_mt_path(*tranche_data), _intervals=intervals
                    )
                    mt = mt.annotate_rows(**ht[mt.row_key])

                    logger.info("Densifying and exporting VCF...")
                    mt = hl.experimental.densify(mt)

                    logger.info("Removing low QUAL variants and * alleles...")
                    info_ht = hl.read_table(info_ht_path(data_source, freeze))
                    mt = mt.filter_rows(
                        (~info_ht[mt.row_key].AS_lowqual)
                        & ((hl.len(mt.alleles) > 1) & (mt.alleles[1] != "*"))
                    )

                    hl.export_vcf(
                        mt,
                        release_vcf_path(*tranche_data, contig=contig),
                        metadata=header_dict,
                    )

            # Export sharded VCF
            if args.parallelize:

                logger.info("Densifying...")
                mt = hl.experimental.densify(mt)

                logger.info("Removing low QUAL variants and * alleles...")
                info_ht = hl.read_table(info_ht_path(data_source, freeze))
                mt = mt.filter_rows(
                    (~info_ht[mt.row_key].AS_lowqual)
                    & ((hl.len(mt.alleles) > 1) & (mt.alleles[1] != "*"))
                )

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
        "--hybrid_pop_map",
        help='Dictionary mapping global populations (keys) to \
        list of all hybrid populations associated with that global population name (values).\n\
        e.g., \'{"nfe": [1,2,3], "afr": [4]}\'',
        type=json.loads,
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
