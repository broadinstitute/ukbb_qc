import argparse
import logging
import pickle
import sys
from typing import Dict, List, Union

import hail as hl

from gnomad.resources.grch38.reference_data import dbsnp
from gnomad.sample_qc.ancestry import POP_NAMES
from gnomad.utils.generic import get_reference_genome
from gnomad.utils.slack import slack_notifications
from gnomad_qc.v2.variant_qc.prepare_data_release import (
    FAF_POPS,
    GROUPS,
    HISTS,
    POPS,
    SEXES,
    SORT_ORDER,
    generic_field_check,
    make_filters_sanity_check_expr,
    make_label_combos,
)
from gnomad_qc.v2.variant_qc.prepare_data_release import (
    EAS_SUBPOPS as GNOMAD_EAS_SUBPOPS,
)
from gnomad_qc.v2.variant_qc.prepare_data_release import (
    NFE_SUBPOPS as GNOMAD_NFE_SUBPOPS,
)
from ukbb_qc.resources.basics import (
    get_checkpoint_path,
    logging_path,
    release_header_path,
    release_mt_path,
    release_vcf_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.slack_creds import slack_token
from ukbb_qc.utils.constants import FORMAT_DICT, INFO_DICT
from ukbb_qc.utils.utils import make_index_dict


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("data_release")
logger.setLevel(logging.INFO)


def make_info_expr(t: Union[hl.MatrixTable, hl.Table]) -> Dict[str, hl.expr.Expression]:
    """
    Makes Hail expression for variant annotations to be included in VCF INFO field.

    :param Table/MatrixTable t: Table/MatrixTable containing variant annotations to be reformatted for VCF export.
    :return: Dictionary containing Hail expressions for relevant INFO annotations.
    :rtype: Dict[str, hl.expr.Expression]
    """
    info_dict = {
        "FS": t.rf.info_FS,
        "InbreedingCoeff": t.rf.inbreeding_coeff,
        "MQ": t.rf.info_MQ,
        "MQRankSum": t.rf.info_MQRankSum,
        "QD": t.rf.info_QD,
        "ReadPosRankSum": t.rf.info_ReadPosRankSum,
        "SOR": t.rf.info_SOR,
        "VQSR_POSITIVE_TRAIN_SITE": t.vqsr.POSITIVE_TRAIN_SITE,
        "VQSR_NEGATIVE_TRAIN_SITE": t.vqsr.NEGATIVE_TRAIN_SITE,
        "AS_BaseQRankSum": t.vqsr.AS_BaseQRankSum,
        "VarDP": t.rf.info_VarDP,
        "AS_VarDP": t.vqsr.AS_VarDP,
        "AS_VQSLOD": t.vqsr.AS_VQSLOD,
        "AS_VQSR_culprit": t.vqsr.AS_culprit,
        "lcr": t.region_flag.lcr,
        "fail_interval_qc": t.region_flag.fail_interval_qc,
        "nonpar": t.nonpar,
        "rf_positive_label": t.rf.tp,
        "rf_negative_label": t.rf.fail_hard_filters,
        "rf_label": t.rf.rf_label,
        "rf_train": t.rf.rf_train,
        "rf_tp_probability": t.rf.rf_probability,
        "transmitted_singleton": t.rf.transmitted_singleton,
        "variant_type": t.allele_info.variant_type,
        "original_alleles": t.allele_info.original_alleles,
        "allele_type": t.allele_info.allele_type,
        "n_alt_alleles": t.allele_info.n_alt_alleles,
        "was_mixed": t.allele_info.was_mixed,
        "has_star": t.allele_info.has_star,
        "pab_max": t.rf.pab_max,
    }
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
            info_dict.update(hist_dict)
    return info_dict


def make_info_dict(
    prefix: str,
    label_groups: Dict[str, str] = None,
    bin_edges: Dict[str, str] = None,
    faf: bool = False,
    popmax: bool = False,
    age_hist_data: str = None,
) -> Dict[str, str]:
    """
    Generate dictionary of Number and Description attributes used for the VCF header during export.
    
    :param str prefix: gnomAD or empty (empty is for UKBB).
    :param dict label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"]).
    :param dict bin_edges: Dictionary keyed by annotation type, with values that reflect the bin edges corresponding to the annotation.
    :param bool faf: If True, use alternate logic to auto-populate dictionary values associated with filter allele frequency annotation.s
    :param bool popmax: If True, use alternate logic to auto-populate dictionary values associated with popmax annotations.
    :param str age_hist_data: Pipe-delimited string of age histograms, from get_age_distributions (required if popmax is set and prefix is empty).
    :return: Dictionary keyed by VCF INFO annotations, where values are dictionaries of Number and Description attributes.
    :rtype: Dict of str: (Dict of str: str)
    """
    info_dict = dict()

    if popmax:
        if "gnomad" in prefix:
            popmax_text = " in gnomAD"
        else:
            popmax_text = ""

            # Only create age hist dictionary for UKBB
            age_hist_dict = {
                "age_hist_het_bin_freq": {
                    "Number": "A",
                    "Description": f"Histogram of ages of heterozygous individuals; bin edges are: {bin_edges['het']}; total number of individuals of any genotype bin: {age_hist_data}",
                },
                "age_hist_het_n_smaller": {
                    "Number": "A",
                    "Description": "Count of age values falling below lowest histogram bin edge for heterozygous individuals",
                },
                "age_hist_het_n_larger": {
                    "Number": "A",
                    "Description": "Count of age values falling above highest histogram bin edge for heterozygous individuals",
                },
                "age_hist_hom_bin_freq": {
                    "Number": "A",
                    "Description": f"Histogram of ages of homozygous alternate individuals; bin edges are: {bin_edges['hom']}; total number of individuals of any genotype bin: {age_hist_data}",
                },
                "age_hist_hom_n_smaller": {
                    "Number": "A",
                    "Description": "Count of age values falling below lowest histogram bin edge for homozygous alternate individuals",
                },
                "age_hist_hom_n_larger": {
                    "Number": "A",
                    "Description": "Count of age values falling above highest histogram bin edge for homozygous alternate individuals",
                },
            }
            info_dict.update(age_hist_dict)

        popmax_dict = {
            f"{prefix}popmax": {
                "Number": "A",
                "Description": f"Population with maximum AF{popmax_text}",
            },
            f"{prefix}AC_popmax": {
                "Number": "A",
                "Description": f"Allele count in the population with the maximum AF{popmax_text}",
            },
            f"{prefix}AN_popmax": {
                "Number": "A",
                "Description": f"Total number of alleles in the population with the maximum AF{popmax_text}",
            },
            f"{prefix}AF_popmax": {
                "Number": "A",
                "Description": f"Maximum allele frequency across populations (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry){popmax_text}",
            },
            f"{prefix}nhomalt_popmax": {
                "Number": "A",
                "Description": f"Count of homozygous individuals in the population with the maximum allele frequency{popmax_text}",
            },
        }
        info_dict.update(popmax_dict)

    else:
        group_types = sorted(label_groups.keys(), key=lambda x: SORT_ORDER.index(x))
        combos = make_label_combos(label_groups)

        for combo in combos:
            combo_fields = combo.split("_")

            if not faf:
                combo_dict = {
                    f"{prefix}AC_{combo}": {
                        "Number": "A",
                        "Description": f"Alternate allele count{make_combo_header_text('for', group_types, combo_fields, prefix)}",
                    },
                    f"{prefix}AN_{combo}": {
                        "Number": "1",
                        "Description": f"Total number of alleles{make_combo_header_text('in', group_types, combo_fields, prefix)}",
                    },
                    f"{prefix}AF_{combo}": {
                        "Number": "A",
                        "Description": f"Alternate allele frequency{make_combo_header_text('in', group_types, combo_fields, prefix)}",
                    },
                    f"{prefix}nhomalt_{combo}": {
                        "Number": "A",
                        "Description": f"Count of homozygous individuals{make_combo_header_text('in', group_types, combo_fields, prefix)}",
                    },
                }
            else:
                combo_dict = {
                    f"{prefix}faf95_{combo}": {
                        "Number": "A",
                        "Description": f"Filtering allele frequency (using Poisson 95% CI) {make_combo_header_text('for', group_types, combo_fields, prefix, faf=True)}",
                    },
                    f"{prefix}faf99_{combo}": {
                        "Number": "A",
                        "Description": f"Filtering allele frequency (using Poisson 99% CI) {make_combo_header_text('for', group_types, combo_fields, prefix, faf=True)}",
                    },
                }
            info_dict.update(combo_dict)
    return info_dict


def make_hist_bin_edges_expr(ht: hl.Table) -> Dict[str, str]:
    """
    Create dictionarys containing variant histogram annotations and their associated bin edges, formatted into a string
    separated by pipe delimiters.

    :param Table ht: Table containing histogram variant annotations.
    :return: Dictionary keyed by histogram annotation name, with corresponding reformatted bin edges for values.
    :rtype: Dict[str, str]
    """
    edges_dict = {
        "het": "|".join(
            map(lambda x: f"{x:.1f}", ht.head(1).age_hist_het.collect()[0].bin_edges)
        ),
        "hom": "|".join(
            map(lambda x: f"{x:.1f}", ht.head(1).age_hist_hom.collect()[0].bin_edges)
        ),
    }
    for hist in HISTS:

        # Parse hists calculated on both raw and adj-filtered data
        for prefix in ["adj_qual_hists", "qual_hists"]:
            hist_name = hist
            if "adj" in prefix:
                hist_name = f"{hist}_adj"
            edges_dict[hist_name] = (
                "|".join(
                    map(lambda x: f"{x:.2f}", ht.take(1)[0][prefix][hist].bin_edges)
                )
                if "ab" in hist
                else "|".join(
                    map(lambda x: str(int(x)), ht.take(1)[0][prefix][hist].bin_edges)
                )
            )
    return edges_dict


def make_hist_dict(bin_edges: Dict, adj: bool) -> Dict[str, str]:
    """
    Generate dictionary of Number and Description attributes to be used in the VCF header, specifically for histogram annotations

    :param dict bin_edges: Dictionary keyed by histogram annotation name, with corresponding string-reformatted bin edges for values
    :param bool adj: Whether to create a header dict for raw or adj qual hists
    :return: Dictionary keyed by VCF INFO annotations, where values are Dictionaries of Number and Description attributes
    :rtype: Dict of str: (Dict of str: str)
    """
    header_hist_dict = {}
    for hist in HISTS:

        # get hists for raw and adj data
        if adj:
            hist = f"{hist}_adj"

        edges = bin_edges[hist]
        hist_fields = hist.split("_")
        hist_text = hist_fields[0].upper()

        if hist_fields[2] == "alt":
            hist_text = hist_text + " in heterozygous individuals"
        if adj:
            hist_text = hist_text + " calculated on high quality genotypes"

        hist_dict = {
            f"{hist}_bin_freq": {
                "Number": "A",
                "Description": f"Histogram for {hist_text}; bin edges are: {edges}",
            },
            f"{hist}_n_smaller": {
                "Number": "A",
                "Description": f"Count of {hist_fields[0].upper()} values falling below lowest histogram bin edge {hist_text}",
            },
            f"{hist}_n_larger": {
                "Number": "A",
                "Description": f"Count of {hist_fields[0].upper()} values falling above highest histogram bin edge {hist_text}",
            },
        }

        header_hist_dict.update(hist_dict)
    return header_hist_dict


def make_filter_dict(ht: hl.Table) -> Dict[str, str]:
    """
    Generates dictionary of Number and Description attributes to be used in the VCF header, specifically for FILTER annotations.

    :param Table ht: Table containing global annotations of the Random Forests SNP and indel cutoffs.
    :return: Dictionary keyed by VCF FILTER annotations, where values are Dictionaries of Number and Description attributes.
    :rtype: Dict[str, str]
    """
    snp_cutoff = hl.eval(ht.globals.rf_globals.rf_snv_cutoff)
    indel_cutoff = hl.eval(ht.globals.rf_globals.rf_indel_cutoff)
    filter_dict = {
        "AC0": {
            "Description": "Allele count is zero after filtering out low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)"
        },
        "InbreedingCoeff": {"Description": "InbreedingCoeff < -0.3"},
        "RF": {
            "Description": f"Failed random forest filtering thresholds of {np_cutoff.min_score}, {indel_cutoff.min_score} (probabilities of being a true positive variant) for SNPs, indels"
        },
        "PASS": {"Description": "Passed all variant filters"},
    }
    return filter_dict


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


def make_combo_header_text(
    preposition: str,
    group_types: List[str],
    combo_fields: List[str],
    prefix: str,
    faf: bool = False,
) -> str:
    """
    Programmatically generate text to populate the VCF header description for a given variant annotation with specific groupings and subset

    :param str preposition: Relevant preposition to precede automatically generated text.
    :param list of str group_types: List of grouping types, e.g. "sex" or "pop".
    :param list of str combo_fields: List of the specific values for each grouping type, for which the text is being generated.
    :param str prefix: gnomad_exomes or gnomad_genomes.
    :param bool faf: If True, use alternate logic to automatically populate descriptions for filter allele frequency annotations.
    :return: String with automatically generated description text for a given set of combo fields.
    :rtype: str
    """
    combo_dict = dict(zip(group_types, combo_fields))
    header_text = " " + preposition

    if "sex" in combo_dict.keys():
        header_text = header_text + " " + combo_dict["sex"]

    if faf:
        header_text = header_text + " samples"

    if "subpop" in combo_dict.keys():
        header_text = header_text + f" of {POP_NAMES[combo_dict['subpop']]} ancestry"
        combo_dict.pop("pop")
    if "pop" in combo_dict.keys():
        header_text = header_text + f" of {POP_NAMES[combo_dict['pop']]} ancestry"

    if "gnomad" in prefix:
        header_text = header_text + " in gnomAD"

    if "group" in group_types:
        if combo_dict["group"] == "raw":
            header_text = header_text + ", before removing low-confidence genotypes"

    if header_text == f" {preposition}":
        header_text = ""

    return header_text


def set_female_y_metrics_to_na(
    t: Union[hl.Table, hl.MatrixTable]
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Set AC, AN, and nhomalt Y variant annotations for females to NA (instead of 0).

    :param Table/MatrixTable t: Table/MatrixTable containing female variant annotations.
    :return: Table/MatrixTable with reset annotations.
    :rtype: hl.Table/hl.MatrixTable
    """
    metrics = list(t.row.info)
    female_metrics = [x for x in metrics if "_female" in x]
    faf_female_metrics = [x for x in female_metrics if "faf" in x]
    female_metrics = [
        x for x in female_metrics if ("nhomalt" in x) or ("AC" in x) or ("AN" in x)
    ]

    female_metrics_dict = {}
    for metric in female_metrics:
        female_metrics_dict.update(
            {
                f"{metric}": hl.cond(
                    t.locus.contig == "chrY", hl.null(hl.tint32), t.info[f"{metric}"]
                )
            }
        )
    for metric in faf_female_metrics:
        female_metrics_dict.update(
            {
                f"{metric}": hl.cond(
                    t.locus.contig == "chrY", hl.null(hl.tfloat32), t.info[f"{metric}"]
                )
            }
        )
    return (
        t.annotate(info=t.info.annotate(**female_metrics_dict))
        if isinstance(t, hl.Table)
        else t.annotate_rows(info=t.info.annotate(**female_metrics_dict))
    )


def sample_sum_check(
    ht: hl.Table,
    prefix: str,
    label_groups: Dict[str, List[str]],
    verbose: bool,
    subpop: bool = None,
) -> None:
    """
    Compute afresh the sum of annotations for a specified group of annotations, and compare to the annotated version;
    display results from checking the sum of the specified annotations in the terminal.

    :param Table ht: Table containing annotations to be summed.
    :param str prefix: Subset of gnomAD
    :param dict label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"]).
    :param bool verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param str subpop: Subpop abbreviation, supplied only if subpopulations are included in the annotation groups being checked.
    :rtype: None
    """
    combo_AC = [ht.info[f"{prefix}AC_{x}"] for x in make_label_combos(label_groups)]
    combo_AN = [ht.info[f"{prefix}AN_{x}"] for x in make_label_combos(label_groups)]
    combo_nhomalt = [
        ht.info[f"{prefix}nhomalt_{x}"] for x in make_label_combos(label_groups)
    ]

    group = label_groups.pop("group")[0]
    alt_groups = "_".join(
        sorted(label_groups.keys(), key=lambda x: SORT_ORDER.index(x))
    )

    annot_dict = {
        f"sum_AC_{group}_{alt_groups}": hl.sum(combo_AC),
        f"sum_AN_{group}_{alt_groups}": hl.sum(combo_AN),
        f"sum_nhomalt_{group}_{alt_groups}": hl.sum(combo_nhomalt),
    }

    ht = ht.annotate(**annot_dict)

    for subfield in ["AC", "AN", "nhomalt"]:
        if not subpop:
            generic_field_check(
                ht,
                (
                    ht.info[f"{prefix}{subfield}_{group}"]
                    != ht[f"sum_{subfield}_{group}_{alt_groups}"]
                ),
                f"adj_{prefix}{subfield}_{group} = sum({subfield}_{group}_{alt_groups})",
                [
                    f"info.{prefix}{subfield}_{group}",
                    f"sum_{subfield}_{group}_{alt_groups}",
                ],
                verbose,
            )
        else:
            generic_field_check(
                ht,
                (
                    ht.info[f"{prefix}{subfield}_{group}_{subpop}"]
                    != ht[f"sum_{subfield}_{group}_{alt_groups}"]
                ),
                f"{prefix}{subfield}_{group}_{subpop} = sum({subfield}_{group}_{alt_groups})",
                [
                    f"info.{prefix}{subfield}_{group}_{subpop}",
                    f"sum_{subfield}_{group}_{alt_groups}",
                ],
                verbose,
            )


def sanity_check_mt(
    mt: hl.MatrixTable,
    subsets: List[str],
    missingness_threshold: float = 0.5,
    verbose: bool = False,
) -> None:
    """
    Perform a battery of sanity checks on a specified group of subsets in a Hail Table containing variant annotations;
    includes summaries of % filter status for different partitions of variants; histogram outlier bin checks; checks on
    AC, AN, and AF annotations; checks that subgroup annotation values add up to the supergroup annotation values;
    checks on sex-chromosome annotations; and summaries of % missingness in variant annotations

    :param MatrixTable mt: MatrixTable containing variant annotations to check
    :param List[str] subsets: List of subsets to be checked
    :param bool verbose: If True, display top values of relevant annotations being checked, regardless of whether check
        conditions are violated; if False, display only top values of relevant annotations if check conditions are violated
    :return: Terminal display of results from the battery of sanity checks
    :rtype: None
    """
    n_samples = mt.count_cols()
    ht = mt.rows()
    n_sites = ht.count()
    contigs = ht.aggregate(hl.agg.collect_as_set(ht.locus.contig))
    logger.info(f"Found {n_sites} sites in contigs {contigs} in {n_samples} samples")
    info_metrics = list(ht.row.info)
    non_info_metrics = list(ht.row)
    non_info_metrics.remove("info")

    logger.info("VARIANT FILTER SUMMARIES:")
    ht_explode = ht.explode(ht.filters)
    logger.info(
        f"hl.agg.counter filters: {ht_explode.aggregate(hl.agg.counter(ht_explode.filters))}"
    )
    ht = ht.annotate(
        is_filtered=ht.filters.length() > 0,
        in_problematic_region=hl.any(
            lambda x: x, [ht.info.lcr, ht.info.fail_interval_qc]
        ),
    )

    # NOTE: make_filters_sanity_check_expr returns a dict with %ages of variants filtered
    ht_filter_check1 = (
        ht.group_by(ht.is_filtered)
        .aggregate(**make_filters_sanity_check_expr(ht))
        .order_by(hl.desc("n"))
    )
    ht_filter_check1.show()

    ht_filter_check2 = (
        ht.group_by(ht.info.allele_type)
        .aggregate(**make_filters_sanity_check_expr(ht))
        .order_by(hl.desc("n"))
    )
    ht_filter_check2.show()

    ht_filter_check3 = (
        ht.group_by(ht.info.allele_type, ht.in_problematic_region)
        .aggregate(**make_filters_sanity_check_expr(ht))
        .order_by(hl.desc("n"))
    )
    ht_filter_check3.show(50, 140)

    ht_filter_check4 = (
        ht.group_by(
            ht.info.allele_type, ht.in_problematic_region, ht.info.n_alt_alleles
        )
        .aggregate(**make_filters_sanity_check_expr(ht))
        .order_by(hl.desc("n"))
    )
    ht_filter_check4.show(50, 140)

    # NOTE: generic field check filters ht based on a certain condition and checks for the number of rows in ht that failed that condition
    # if the number of fails is 0, then the ht passes that check; otherwise it fails
    # its inputs are:
    #   ht, condition expression (the condition to be filtered to), check description (name of check),
    #   display fields (fields of ht to show), verbose (show values of checks that pass; default is to only show when they fail)
    logger.info("HISTOGRAM CHECKS:")
    for hist in HISTS:
        for suffix in ["", "adj"]:
            if suffix == "adj":
                logger.info("Checking adj qual hists")
                hist = f"{hist}_{suffix}"
            else:
                logger.info("Checking raw qual hists")

            # Check subfield == 0
            generic_field_check(
                ht,
                (ht.info[f"{hist}_n_smaller"] != 0),
                f"{hist}_n_smaller == 0",
                [f"info.{hist}_n_smaller"],
                verbose,
            )
            if hist not in [
                "dp_hist_alt",
                "dp_hist_all",
            ]:  # NOTE: DP hists can have nonzero values in n_larger bin
                generic_field_check(
                    ht,
                    (ht.info[f"{hist}_n_larger"] != 0),
                    f"{hist}_n_larger == 0",
                    [f"info.{hist}_n_larger"],
                    verbose,
                )

    logger.info("RAW AND ADJ CHECKS:")
    for subfield in ["AC", "AN", "AF"]:
        # Check AC, AN, AF > 0
        generic_field_check(
            ht,
            (ht.info[f"raw_{subfield}_raw"] <= 0),
            f"raw_{subfield}_raw > 0",
            [f"info.raw_{subfield}_raw"],
            verbose,
        )
        generic_field_check(
            ht,
            (ht.info[f"adj_{subfield}_adj"] < 0),
            f"adj_{subfield}_adj >= 0",
            [f"info.adj_{subfield}_adj", "filters"],
            verbose,
        )

    for subset in subsets:
        if subset == "":
            for subfield in ["AC", "AN", "nhomalt"]:
                # Check AC_raw >= AC adj
                generic_field_check(
                    ht,
                    (ht.info[f"raw_{subfield}_raw"] < ht.info[f"adj_{subfield}_adj"]),
                    f"raw_{subfield}_raw >= adj_{subfield}_adj",
                    [f"info.raw_{subfield}_raw", f"info.adj_{subfield}_adj"],
                    verbose,
                )
        else:
            for subfield in ["AC", "AN", "nhomalt"]:
                # Check AC_raw >= AC adj
                generic_field_check(
                    ht,
                    (
                        ht.info[f"{subset}{subfield}_raw"]
                        < ht.info[f"{subset}{subfield}_adj"]
                    ),
                    f"{subset}{subfield}_raw >= {subfield}_adj",
                    [f"info.{subset}{subfield}_raw", f"info.{subset}{subfield}_adj"],
                    verbose,
                )

    logger.info("FREQUENCY CHECKS:")
    for subset in subsets:
        if subset == "":
            logger.info("raw checks -- gnomad")
            for subfield in ["AC", "AN", "nhomalt"]:
                generic_field_check(
                    ht,
                    (
                        ht.info[f"gnomad_exomes_{subfield}_raw"]
                        == ht.info[f"gnomad_genomes_{subfield}_raw"]
                    ),
                    f"gnomad_exomes_{subfield}_raw == gnomad_genomes_{subfield}_raw",
                    [
                        f"info.gnomad_exomes_{subfield}_raw",
                        f"info.gnomad_genomes_{subfield}_raw",
                    ],
                    verbose,
                )
                logger.info("adj checks -- gnomad")
                generic_field_check(
                    ht,
                    (
                        ht.info[f"gnomad_exomes_{subfield}_adj"]
                        == ht.info[f"gnomad_genomes_{subfield}_adj"]
                    ),
                    f"gnomad_exomes_{subfield}_adj == gnomad_genomes_{subfield}_adj",
                    [
                        f"info.gnomad_exomes_{subfield}_adj",
                        f"info.gnomad_genomes_{subfield}_adj",
                    ],
                    verbose,
                )
        else:
            for subfield in ["AC", "AN", "nhomalt"]:
                logger.info("raw checks -- gnomad/ukb")
                generic_field_check(
                    ht,
                    (
                        ht.info[f"raw_{subfield}_raw"]
                        == ht.info[f"{subset}{subfield}_raw"]
                    ),
                    f"raw_{subfield}_raw == {subset}{subfield}_raw",
                    [f"info.raw_{subfield}_raw", f"info.{subset}{subfield}_raw"],
                    verbose,
                )
                logger.info("adj checks -- gnomad/ukb")
                generic_field_check(
                    ht,
                    (
                        ht.info[f"adj_{subfield}_adj"]
                        == ht.info[f"{subset}{subfield}_adj"]
                    ),
                    f"adj_{subfield}_adj == {subset}{subfield}_adj",
                    [f"info.adj_{subfield}_adj", f"info.{subset}{subfield}_adj"],
                    verbose,
                )

    freq_counts = ht.aggregate(
        hl.struct(
            total_defined_gnomad_wes_AC=hl.agg.count_where(
                hl.is_defined(ht.info.gnomad_exomes_AC_adj)
            ),
            total_defined_gnomad_wes_AC_raw=hl.agg.count_where(
                hl.is_defined(ht.info.gnomad_exomes_AC_raw)
            ),
            total_defined_gnomad_wgs_AC=hl.agg.count_where(
                hl.is_defined(ht.info.gnomad_genomes_AC_adj)
            ),
            total_defined_gnomad_wgs_AC_raw=hl.agg.count_where(
                hl.is_defined(ht.info.gnomad_genomes_AC_raw)
            ),
            total_defined_ukb_AC=hl.agg.count_where(hl.is_defined(ht.info.adj_AC_adj)),
            total_defined_ukb_AC_raw=hl.agg.count_where(
                hl.is_defined(ht.info.raw_AC_raw)
            ),
        )
    )
    logger.info(f"Frequency spot check counts: {freq_counts}")

    logger.info("SAMPLE SUM CHECKS:")
    for subset in subsets:
        # Check if pops are present
        if subset == "gnomad":
            pop_adjusted = list(
                set([x for x in info_metrics if subset in x and "raw" not in x])
            )
        else:
            pop_adjusted = list(
                set(
                    [
                        x
                        for x in info_metrics
                        if (
                            ("adj" in x or "raw" in x)
                            and ("gnomad" not in x)
                            and ("raw" not in x)
                        )
                    ]
                )
            )

        pop_adjusted = [
            i.replace("adj_", "").replace("_adj", "").replace("_adj_", "")
            for i in pop_adjusted
        ]
        found = []
        for i in pop_adjusted:
            if subset == "gnomad":
                for z in POP_NAMES:
                    if z in i:
                        found.append(z)
            else:
                for z in POPS:
                    if z in i:
                        found.append(z)
        missing_pops = set(POPS) - set(found)
        if len(missing_pops) != 0:
            logger.warning(f"Missing {missing_pops} pops in {subset} subset!")

        if subset == "gnomad":
            sample_sum_check(ht, subset, dict(group=["adj"], pop=POP_NAMES), verbose)
            sample_sum_check(ht, subset, dict(group=["adj"], sex=SEXES), verbose)
            sample_sum_check(
                ht, subset, dict(group=["adj"], pop=POP_NAMES, sex=SEXES), verbose
            )

            # Adjust subpops to those found in subset
            nfe_subpop_adjusted = list(
                set([x for x in pop_adjusted if "nfe_" in x and "male" not in x])
            )
            if nfe_subpop_adjusted != []:
                sample_sum_check(
                    ht,
                    subset,
                    dict(group=["adj"], pop=["nfe"], subpop=GNOMAD_NFE_SUBPOPS),
                    verbose,
                    subpop="nfe",
                )
            eas_subpop_adjusted = list(
                set(
                    [
                        x
                        for x in pop_adjusted
                        if subset in x and "eas_" in x and "male" not in x
                    ]
                )
            )
            if eas_subpop_adjusted != []:
                sample_sum_check(
                    ht,
                    subset,
                    dict(group=["adj"], pop=["eas"], subpop=GNOMAD_EAS_SUBPOPS),
                    verbose,
                    subpop="eas",
                )
        else:
            subset = "adj_"  # hacky add; UKB fields are weirdly named adj_AC_adj
            sample_sum_check(ht, subset, dict(group=["adj"], pop=POPS), verbose)
            sample_sum_check(ht, subset, dict(group=["adj"], sex=SEXES), verbose)
            sample_sum_check(
                ht, subset, dict(group=["adj"], pop=POPS, sex=SEXES), verbose
            )

    logger.info("SEX CHROMOSOME ANNOTATION CHECKS:")
    female_metrics = [x for x in info_metrics if "_female" in x]
    male_metrics = [x for x in info_metrics if "_male" in x]

    if "chrY" in contigs:
        logger.info("Check values of female metrics for Y variants are NA:")
        ht_y = hl.filter_intervals(ht, [hl.parse_locus_interval("chrY")])
        metrics_values = {}
        for metric in female_metrics:
            metrics_values[metric] = hl.agg.collect_as_set(ht_y.info[metric])
        output = ht_y.aggregate(hl.struct(**metrics_values))
        for metric, values in dict(output).items():
            if values == {None}:
                logger.info(f"PASSED {metric} = {None} check for Y variants")
            else:
                logger.info(f"FAILED Y check: Found {values} in {metric}")

    logger.info("Check values of male nhomalt metrics for X nonpar variants are 0:")
    ht_x = hl.filter_intervals(ht, [hl.parse_locus_interval("chrX")])
    ht_xnonpar = ht_x.filter(ht_x.locus.in_x_nonpar())
    n = ht_xnonpar.count()
    logger.info(
        f"Found {n} X nonpar sites"
    )  # Lots of values found in male X nonpar sites

    male_metrics = [x for x in male_metrics if "nhomalt" in x]
    metrics_values = {}
    for metric in male_metrics:
        metrics_values[metric] = hl.agg.collect_as_set(ht_xnonpar.info[metric])
    output = ht_xnonpar.aggregate(hl.struct(**metrics_values))
    for metric, values in dict(output).items():
        if values == {0}:
            logger.info(f"PASSED {metric} = 0 check for X nonpar variants")
        else:
            logger.info(f"FAILED X nonpar check: Found {values} in {metric}")

    logger.info("Check (nhomalt == nhomalt_female) for X nonpar variants:")
    female_metrics = [x for x in female_metrics if "nhomalt" in x]
    for metric in female_metrics:
        standard_field = metric.replace("_female", "")
        generic_field_check(
            ht_xnonpar,
            (ht_xnonpar.info[f"{metric}"] != ht_xnonpar.info[f"{standard_field}"]),
            f"{metric} == {standard_field}",
            [f"info.{metric}", f"info.{standard_field}"],
            verbose,
        )

    logger.info("MISSINGNESS CHECKS:")
    metrics_frac_missing = {}
    for x in info_metrics:
        metrics_frac_missing[x] = hl.agg.sum(hl.is_missing(ht.info[x])) / n_sites
    for x in non_info_metrics:
        metrics_frac_missing[x] = hl.agg.sum(hl.is_missing(ht[x])) / n_sites
    output = ht.aggregate(hl.struct(**metrics_frac_missing))

    n_fail = 0
    for metric, value in dict(output).items():
        message = f"missingness check for {metric}: {100 * value}%missing"
        if value > missingness_threshold:
            logger.info(f"FAILED {message}")
            n_fail += 1
        else:
            logger.info(f"Passed {message}")
    logger.info(f"{n_fail} missing metrics checks failed")


def main(args):

    hl.init(log="/release.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)

    try:

        if args.prepare_vcf_mt:
            logger.info("Starting VCF process...")
            mt = hl.read_matrix_table(release_mt_path(*tranche_data))

            logger.info("Dropping cohort frequencies (last index of freq)...")
            mt = mt.annotate_rows(freq=mt.freq[-1])
            mt = mt.annotate_globals(freq_meta=mt.freq_meta[-1])

            logger.info("Making histogram bin edges...")
            bin_edges = make_hist_bin_edges_expr(mt.rows())

            logger.info("Getting age hist data...")
            age_hist_data = hl.eval(mt.age_distribution)

            logger.info("Making INFO dict for VCF...")
            subset_list = ["", "gnomad_exomes_", "gnomad_genomes_"]  # empty for ukbb
            for subset in subset_list:
                INFO_DICT.update(make_info_dict(subset, dict(group=GROUPS)))
                INFO_DICT.update(make_info_dict(subset, dict(group=["adj"], sex=SEXES)))
                INFO_DICT.update(make_info_dict(subset, dict(group=["adj"]), faf=True))
                INFO_DICT.update(
                    make_info_dict(subset, dict(group=["adj"], sex=SEXES), faf=True)
                )
                INFO_DICT.update(make_info_dict(subset, dict(group=["adj"], pop=POPS)))
                INFO_DICT.update(
                    make_info_dict(subset, dict(group=["adj"], pop=POPS, sex=SEXES))
                )
                INFO_DICT.update(
                    make_info_dict(subset, dict(group=["adj"], pop=FAF_POPS), faf=True)
                )
                INFO_DICT.update(
                    make_info_dict(
                        subset, dict(group=["adj"], pop=FAF_POPS, sex=SEXES), faf=True,
                    )
                )

                if "gnomad" in subset:
                    INFO_DICT.update(make_info_dict(subset, popmax=True))
                    INFO_DICT.update(
                        make_info_dict(
                            subset,
                            dict(group=["adj"], pop=["nfe"], subpop=GNOMAD_NFE_SUBPOPS),
                        )
                    )
                    INFO_DICT.update(
                        make_info_dict(
                            subset,
                            dict(group=["adj"], pop=["eas"], subpop=GNOMAD_EAS_SUBPOPS),
                        )
                    )

                else:
                    INFO_DICT.update(
                        make_info_dict(
                            subset,
                            bin_edges=bin_edges,
                            popmax=True,
                            age_hist_data="|".join(str(x) for x in age_hist_data),
                        )
                    )
            INFO_DICT.update(make_hist_dict(bin_edges, adj=False))

            # Adjust keys to remove adj tags before exporting to VCF
            # new_info_dict = {i.replace('_adj', '').replace('adj_', '').replace('_adj_', '').replace('.', '_').replace('raw_', ''): j for i,j in INFO_DICT.items()}
            new_info_dict = {
                i.replace("adj_", "").replace("_adj", "").replace("_adj_", ""): j
                for i, j in INFO_DICT.items()
            }
            new_info_dict.update(make_hist_dict(bin_edges, adj=True))

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
            mt = set_female_y_metrics_to_na(mt)

            # Select relevant fields for VCF export
            mt = mt.select_rows("info", "filters", "rsid", "qual")
            new_info_dict.update({"vep": {"Description": hl.eval(mt.vep_csq_header)}})
            header_dict = {
                "info": new_info_dict,
                "filter": make_filter_dict(mt.rows()),
                "format": FORMAT_DICT,
            }

            logger.info("Saving header dict to pickle...")
            with hl.hadoop_open(release_header_path(*tranche_data), "wb") as p:
                pickle.dump(header_dict, p, protocol=pickle.HIGHEST_PROTOCOL)
            mt = mt.annotate_rows(info=mt.info.annotate(vep=mt.vep_csq))
            mt.write(release_mt_path(*tranche_data, temp=True), args.overwrite)

        if args.sanity_check:
            subset_list = ["", "gnomad_exomes_", "gnomad_genomes_"]
            mt = hl.read_matrix_table(release_mt_path(*tranche_data, temp=True))
            mt = set_female_y_metrics_to_na(mt)
            sanity_check_mt(
                mt, subset_list, missingness_threshold=0.5, verbose=args.verbose
            )

        if args.prepare_release_vcf:
            if not args.per_chromosome and not args.parallelize:
                logger.error("Need to choose how to export the release VCF. Exiting...")
                sys.exit(1)

            mt = hl.read_matrix_table(release_mt_path(*tranche_data, temp=True))
            logger.info("Reading header dict from pickle...")
            with hl.hadoop_open(release_header_path(*tranche_data), "rb") as p:
                header_dict = pickle.load(p)

            # Reformat names to remove "adj_" (this is a carryover from code that labeled everything "gnomad")
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
            for item in ["info", "filters", "format"]:
                if item == "info":
                    annots = new_row_annots
                elif item == "format":
                    annots = list(mt.entry)
                else:
                    item = "filter"
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
                    mt = hl.filter_intervals(
                        mt, [hl.parse_locus_interval(contig)]
                    )
                    mt._calculate_new_partitions(10000)
                    intervals = [i for i in intervals if i.start.locus.contig==contig]
                    mt = hl.read_matrix_table(mt_path, _intervals=intervals)

                    hl.export_vcf(
                        contig_mt,
                        release_vcf_path(*tranche_data, contig=contig),
                        metadata=header_dict,
                    )

            # Export sharded VCF
            if args.parallelize:
                mt = mt.naive_coalesce(5000)
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
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
