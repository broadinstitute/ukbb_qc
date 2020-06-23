import argparse
import logging
import pickle
from typing import Dict, List, Union

import hail as hl

from gnomad.resources.grch37.gnomad import (
    POP_NAMES,
)  # Note: this is in gnomad twice? Importing from larger dict
from gnomad.resources.grch38.reference_data import dbsnp, lcr_intervals
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
    index_globals,
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
    get_ukbb_data,
    logging_path,
    release_ht_path,
    release_mt_path,
    release_vcf_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.variant_qc import var_annotations_ht_path
from ukbb_qc.slack_creds import slack_token
from ukbb_qc.utils.constants import FORMAT_DICT, INFO_DICT
from ukbb_qc.utils.utils import get_age_ht


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("data_release")
logger.setLevel(logging.INFO)


def flag_problematic_regions(
    t: Union[hl.Table, hl.MatrixTable]
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Annotate HT/MT with `region_flag` struct. Struct contains flags for problematic regions.
    
    .. note:: 

        No hg38 resources for decoy, segdup, or self chain yet.

    :param Table/MatrixTable t: Input Table/MatrixTable.
    :return: Table/MatrixTable with `region_flag` struct row annotation.
    :rtype: Union[hl.Table, hl.MatrixTable]
    """
    lcr_ht = lcr_intervals.ht()

    if isinstance(t, hl.Table):
        t = t.annotate(
            region_flag=hl.struct(
                lcr=hl.is_defined(lcr_ht[t.locus]),
                fail_interval_qc=~(t.rf.interval_qc_pass),
            )
        )
    else:
        t = t.annotate_rows(
            region_flag=hl.struct(
                lcr=hl.is_defined(lcr_ht[t.locus]),
                fail_interval_qc=~(t.rf.interval_qc_pass),
            )
        )
    return t


def prepare_annotations(
    mt: hl.MatrixTable,
    freq_ht: hl.Table,
    rf_ht: hl.Table,
    vep_ht: hl.Table,
    index_dict: Dict,
    allele_ht: hl.Table,
    vqsr_ht: hl.Table,
) -> hl.MatrixTable:
    """
    Join all Tables with variant annotations, keeping only variants with non-zero AC.

    :param MatrixTable mt: MatrixTable to be annotated.
    :param Table freq_ht: Table with frequency annotations.
    :param Table rf_ht: Table with random forest variant annotations.
    :param Table vep_ht: Table with VEP variant annotations.
    :param dict index_dict: Dictionary containing index values for each entry in the frequency Table freq array, keyed by metadata label.
    :param Table allele_ht: Table containing allele annotations.
    :param Table vqsr_ht: Table containing VQSR annotations.
    :return: Table containing joined annotations.
    :rtype: hl.Table
    """
    logger.info("Removing unnecessary annotations from annotation Tables...")
    rf_ht = rf_ht.select(
        "info_FS",
        "inbreeding_coeff",
        "info_MQ",
        "info_MQRankSum",
        "info_QD",
        "info_ReadPosRankSum",
        "info_SOR",
        "tp",
        "fail_hard_filters",
        "rf_label",
        "rf_train",
        "rf_probability",
        "transmitted_singleton",
        "pab_max",
        "info_VarDP",
        "interval_qc_pass",
        "filters",
    )
    vep_ht = vep_ht.transmute(vep=vep_ht.vep.drop("colocated_variants"))

    logger.info(
        "Removing unnecessary fields from VQSR HT and splitting AS annotations..."
    )
    vqsr_ht = vqsr_ht.transmute(
        info=vqsr_ht.info.select(
            "NEGATIVE_TRAIN_SITE",
            "POSITIVE_TRAIN_SITE",
            "AS_VQSLOD",
            "AS_culprit",
            "AS_BaseQRankSum",
            "AS_VarDP",
        )
    ).drop("rsid", "was_split")
    vqsr_ht = vqsr_ht.annotate(
        info=vqsr_ht.info.annotate(
            **{
                f: [vqsr_ht.info[f][vqsr_ht.a_index - 1]]
                for f in vqsr_ht.info
                if f.startswith("AC")
                or (f.startswith("AS_") and not f == "AS_SB_TABLE")
            }
        )
    )
    vqsr_ht = vqsr_ht.annotate(info=vqsr_ht.info.annotate(qual=vqsr_ht.qual))
    freq_ht = freq_ht.drop(
        "InbreedingCoeff"
    )  # NOTE: Will need to drop InbreedingCoeff annotation here in 500K, but this isn't present in 300K

    logger.info("Filtering out low QUAL variants...")
    mt = mt.filter_rows(~vqsr_ht[mt.row_key].filters.contains("LowQual"))
    logger.info(f"Count after filtering out low QUAL variants {mt.count()}")

    logger.info(f"Frequency HT before filtering out chrM and AC < 1: {freq_ht.count()}")
    freq_ht = hl.filter_intervals(
        freq_ht, [hl.parse_locus_interval("chrM")], keep=False
    )
    raw_idx = index_dict["raw"]
    freq_ht = freq_ht.filter(freq_ht.freq[raw_idx].AC <= 0, keep=False)
    logger.info(
        f"Frequency HT count after filtering out chrM and AC < 1: {freq_ht.count()}"
    )

    logger.info(
        "Annotating MT with random forest HT and flagging problematic regions..."
    )
    mt = mt.annotate_rows(rf=rf_ht[mt.row_key])
    mt = flag_problematic_regions(mt)
    mt = mt.transmute_rows(rf=mt.rf.drop("interval_qc_pass"))

    logger.info("Annotating MT with frequency information...")
    mt = mt.annotate_rows(**freq_ht[mt.row_key])
    mt = mt.annotate_globals(**freq_ht.index_globals())

    dbsnp_ht = dbsnp.ht().select("rsid")
    mt = mt.annotate_rows(
        vep=vep_ht[mt.row_key].vep,
        allele_info=allele_ht[mt.row_key].allele_data,
        vqsr=vqsr_ht[mt.row_key].info,
        rsid=dbsnp_ht[mt.row_key],
    )
    mt = mt.annotate_globals(rf_globals=rf_ht.index_globals())
    logger.info(f"MT count: {mt.count()}")
    return mt


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


def make_freq_meta_index_dict(
    freq_meta: List[str], gnomad: bool, faf: bool
) -> Dict[str, int]:
    """
    Makes a dictionary of the entries in the frequency array annotation, where keys are the grouping combinations and the values
    are the 0-based integer indices.

    :param List[str] freq_meta: Ordered list containing string entries describing all the grouping combinations contained in the
        frequency array annotation.
    :param bool gnomad: Whether to index a list from gnomAD.
    :return: Dictionary keyed by grouping combinations in the frequency array, with values describing the corresponding index
        of each grouping entry in the frequency array
    :rtype: Dict[str, int]
    """
    index_dict = index_globals(freq_meta, dict(group=GROUPS))

    if faf:
        index_dict.update(index_globals(freq_meta, dict(group=GROUPS, pop=FAF_POPS)))

    else:
        index_dict.update(index_globals(freq_meta, dict(group=GROUPS, sex=SEXES)))
        if gnomad:
            index_dict.update(
                index_globals(freq_meta, dict(group=GROUPS, pop=POP_NAMES))
            )
            index_dict.update(
                index_globals(
                    freq_meta,
                    dict(group=GROUPS, pop=["nfe"], subpop=GNOMAD_NFE_SUBPOPS),
                )
            )
            index_dict.update(
                index_globals(
                    freq_meta,
                    dict(group=GROUPS, pop=["eas"], subpop=GNOMAD_EAS_SUBPOPS),
                )
            )
        else:
            index_dict.update(index_globals(freq_meta, dict(group=GROUPS, pop=POPS)))
            index_dict.update(
                index_globals(freq_meta, dict(group=GROUPS, pop=POPS, sex=SEXES))
            )

    return index_dict


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
                        "Description": "Alternate allele count{}".format(
                            make_combo_header_text(
                                "for", group_types, combo_fields, prefix
                            )
                        ),
                    },
                    f"{prefix}AN_{combo}": {
                        "Number": "1",
                        "Description": "Total number of alleles{}".format(
                            make_combo_header_text(
                                "in", group_types, combo_fields, prefix
                            )
                        ),
                    },
                    f"{prefix}AF_{combo}": {
                        "Number": "A",
                        "Description": "Alternate allele frequency{}".format(
                            make_combo_header_text(
                                "in", group_types, combo_fields, prefix
                            )
                        ),
                    },
                    f"{prefix}nhomalt_{combo}": {
                        "Number": "A",
                        "Description": "Count of homozygous individuals{}".format(
                            make_combo_header_text(
                                "in", group_types, combo_fields, prefix
                            )
                        ),
                    },
                }
            else:
                combo_dict = {
                    f"{prefix}faf95_{combo}": {
                        "Number": "A",
                        "Description": "Filtering allele frequency (using Poisson 95% CI) {}".format(
                            make_combo_header_text(
                                "for", group_types, combo_fields, prefix, faf=True
                            )
                        ),
                    },
                    f"{prefix}faf99_{combo}": {
                        "Number": "A",
                        "Description": "Filtering allele frequency (using Poisson 99% CI) {}".format(
                            make_combo_header_text(
                                "for", group_types, combo_fields, prefix, faf=True
                            )
                        ),
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


def make_index_dict(
    t: Union[hl.MatrixTable, hl.Table], freq_meta_str: str
) -> Dict[str, int]:
    """
    Create a look-up Dictionary for entries contained in the frequency annotation array.

    :param Table ht: Table or MatrixTable containing freq_meta global annotation to be indexed
    :param str freq_meta: freq_meta global annotation to be indexed (freq_meta, gnomad_exomes_freq_meta, or gnomad_genomes_freq_meta)
    :return: Dictionary keyed by grouping combinations in the frequency array, with values describing the corresponding index
        of each grouping entry in the frequency array
    :rtype: Dict of str: int
    """
    freq_meta = hl.eval(t.globals[freq_meta_str])
    # check if indexing gnomAD data
    if "gnomad" in freq_meta_str:
        index_dict = make_freq_meta_index_dict(freq_meta, gnomad=True, faf=False)
    else:
        index_dict = make_freq_meta_index_dict(freq_meta, gnomad=False, faf=False)
    return index_dict


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
            "Description": "Failed random forest filtering thresholds of {0}, {1} (probabilities of being a true positive variant) for SNPs, indels".format(
                snp_cutoff.min_score, indel_cutoff.min_score
            )
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


def get_age_distributions(ht: hl.Table) -> str:
    """
    Get background distribution of sample ages (using field 21022, age at recruitment).

    :param Table ht: Table containing samples and sample ages.
    :return: Pipe-delimited string with ages in pre-determined bins (<30, 30-35, ..., 75-80, 80+).
    :rtype: str
    """
    age_hist_data = ht.aggregate(hl.agg.hist(ht.age, 30, 80, 10))
    age_hist_data.bin_freq.insert(0, age_hist_data.n_smaller)
    age_hist_data.bin_freq.append(age_hist_data.n_larger)
    return age_hist_data.bin_freq


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
    logger.info(freq_counts)

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

    # tranche 2 dbsnp fix
    dbsnp_ht = dbsnp.ht().select("rsid")
    ht = ht.drop("rsid")
    ht = ht.annotate(**dbsnp_ht[ht.key])

    metrics_frac_missing = {}
    for x in info_metrics:
        metrics_frac_missing[x] = hl.agg.sum(hl.is_missing(ht.info[x])) / n_sites
    for x in non_info_metrics:
        metrics_frac_missing[x] = hl.agg.sum(hl.is_missing(ht[x])) / n_sites
    output = ht.aggregate(hl.struct(**metrics_frac_missing))

    n_fail = 0
    for metric, value in dict(output).items():
        if value > missingness_threshold:
            logger.info(
                "FAILED missing check for {}: {}% missing".format(metric, 100 * value)
            )
            n_fail += 1
        else:
            logger.info(
                "Missingness check for {}: {}% missing".format(metric, 100 * value)
            )
    logger.info("{} missing metrics checks failed".format(n_fail))


def main(args):

    hl.init(log="/release.log", default_reference="GRCh38")
    data_source = args.data_source
    freeze = args.freeze

    try:
        if args.prepare_internal_mt:
            logger.info("Getting age hist data...")
            age_hist_data = get_age_distributions(get_age_ht(freeze))

            logger.info("Getting raw MT and dropping all unnecessary entries...")
            # NOTE reading in raw to be able to return all samples/variants
            mt = get_ukbb_data(
                data_source,
                freeze,
                key_by_locus_and_alleles=args.key_by_locus_and_alleles,
                split=False,
                raw=True,
                repartition=args.repartition,
                n_partitions=args.raw_partitions,
            ).select_entries(
                "DP",
                "GQ",
                "LA",
                "LAD",
                "LGT",
                "LPGT",
                "LPL",
                "MIN_DP",
                "PID",
                "RGQ",
                "SB",
            )

            logger.info("Splitting raw MT...")
            mt = hl.experimental.sparse_split_multi(mt)
            mt = mt.select_entries(
                "GT", "GQ", "DP", "AD", "MIN_DP", "PGT", "PID", "PL", "SB"
            )

            logger.info("Reading in all variant annotation tables...")
            freq_ht = hl.read_table(
                var_annotations_ht_path(data_source, freeze, "join_freq")
            )
            rf_ht = hl.read_table(var_annotations_ht_path(data_source, freeze, "rf"))
            vep_ht = hl.read_table(var_annotations_ht_path(data_source, freeze, "vep"))
            allele_ht = hl.read_table(
                var_annotations_ht_path(data_source, freeze, "allele_data")
            )
            vqsr_ht = hl.read_table(
                var_annotations_ht_path(data_source, freeze, "vqsr")
            )

            logger.info("Adding annotations...")
            mt = prepare_annotations(
                mt,
                freq_ht,
                rf_ht,
                vep_ht,
                make_index_dict(freq_ht, "freq_meta"),
                allele_ht,
                vqsr_ht,
            )
            mt = mt.annotate_globals(age_distribution=age_hist_data)
            mt.write(release_mt_path(data_source, freeze), args.overwrite)

        if args.prepare_vcf_mt:
            logger.info("Starting VCF process...")
            mt = hl.read_matrix_table(release_mt_path(data_source, freeze))

            logger.info("Dropping cohort frequencies (last index of freq)...")
            mt = mt.annotate_rows(freq=mt.freq[-1])
            freq_meta = mt.freq_meta[-1]
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

                if "gnomad" in subset:
                    INFO_DICT.update(
                        make_info_dict(subset, dict(group=["adj"], pop=POP_NAMES))
                    )
                    INFO_DICT.update(make_info_dict(subset, popmax=True))
                    INFO_DICT.update(
                        make_info_dict(
                            subset, dict(group=["adj"], pop=POP_NAMES, sex=SEXES)
                        )
                    )
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
                    INFO_DICT.update(
                        make_info_dict(
                            subset, dict(group=["adj"], pop=FAF_POPS), faf=True
                        )
                    )
                    INFO_DICT.update(
                        make_info_dict(
                            subset,
                            dict(group=["adj"], pop=FAF_POPS, sex=SEXES),
                            faf=True,
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
                    INFO_DICT.update(
                        make_info_dict(subset, dict(group=["adj"], pop=POPS))
                    )
                    INFO_DICT.update(
                        make_info_dict(subset, dict(group=["adj"], pop=POPS, sex=SEXES))
                    )
                    INFO_DICT.update(
                        make_info_dict(
                            subset, dict(group=["adj"], pop=FAF_POPS), faf=True
                        )
                    )
                    INFO_DICT.update(
                        make_info_dict(
                            subset,
                            dict(group=["adj"], pop=FAF_POPS, sex=SEXES),
                            faf=True,
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
            logger.info(new_info_dict["faf95"])
            logger.info(new_info_dict["gnomad_exomes_faf95"])
            logger.info(new_info_dict["faf95_nfe"])
            logger.info(new_info_dict["faf95_nfe_male"])
            logger.info(new_info_dict["faf95_0"])
            logger.info(new_info_dict["gnomad_exomes_popmax"])
            logger.info(new_info_dict["gnomad_exomes_AC_popmax"])
            logger.info(new_info_dict["dp_hist_alt_adj_bin_freq"])
            logger.info(new_info_dict["dp_hist_alt_bin_freq"])
            logger.info(new_info_dict["AC"])
            logger.info(new_info_dict["AC_raw"])
            logger.info(new_info_dict["AC_nfe"])
            logger.info(new_info_dict["gnomad_exomes_AC"])

            # add non par annotation back
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

            # Add VEP annotations
            # TODO: After 200K release change this to use vep_csq annotation on the vep HT
            vep_csq_ht = hl.read_table(
                var_annotations_ht_path(data_source, freeze, "vep_csq")
            )
            new_info_dict.update(
                {"vep": {"Description": hl.eval(vep_csq_ht.globals.vep_csq_header)}}
            )
            header_dict = {
                "info": new_info_dict,
                "filter": make_filter_dict(mt.rows()),
                "format": FORMAT_DICT,
            }

            logger.info("Saving header dict to pickle")
            pickle_file = (
                f"gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/header_dict.pickle"
            )
            with hl.hadoop_open(pickle_file, "wb") as p:
                pickle.dump(header_dict, p, protocol=pickle.HIGHEST_PROTOCOL)
            mt = mt.annotate_rows(info=mt.info.annotate(vep=vep_csq_ht[mt.row_key].vep))
            mt.write(release_mt_path(data_source, freeze, temp=True), args.overwrite)

        if args.sanity_check:
            subset_list = ["", "gnomad_exomes_", "gnomad_genomes_"]
            mt = hl.read_matrix_table(release_mt_path(data_source, freeze, temp=True))
            mt = set_female_y_metrics_to_na(mt)
            sanity_check_mt(
                mt, subset_list, missingness_threshold=0.5, verbose=args.verbose
            )

        if args.prepare_browser_ht:
            mt = hl.read_matrix_table(release_mt_path(data_source, freeze))
            logger.info(f"mt count: {mt.count()}")

            logger.info("Adding missing filters field (tranche 2 fix)")
            rf_ht = hl.read_table(var_annotations_ht_path(data_source, freeze, "rf"))
            rf_ht = rf_ht.select("filters")
            mt = mt.annotate_rows(**rf_ht[mt.row_key])

            logger.info("Adding missing qual field (tranche 2 fix)")
            vqsr_ht = hl.read_table(
                var_annotations_ht_path(data_source, freeze, "vqsr")
            )
            vqsr_ht = vqsr_ht.select("qual")
            mt = mt.annotate_rows(**vqsr_ht[mt.row_key])

            logger.info("Adding rsids (tranche 2 fix)")
            dbsnp_ht = dbsnp.ht().select("rsid")
            mt = mt.drop("rsid")
            mt = mt.annotate_rows(**dbsnp_ht[mt.row_key])

            ht = mt.rows()
            logger.info(f"ht count: {ht.count()}")
            ht = ht.naive_coalesce(20000)
            ht.write(release_ht_path(data_source, freeze), args.overwrite)

        if args.prepare_release_vcf:

            mt = hl.read_matrix_table(release_mt_path(data_source, freeze, temp=True))
            logger.info("Reading header dict from pickle")
            pickle_file = (
                f"gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/header_dict.pickle"
            )
            with hl.hadoop_open(pickle_file, "rb") as p:
                header_dict = pickle.load(p)

            # tranche 2 set female faf on chrY to NA fix
            mt = set_female_y_metrics_to_na(mt)

            # tranche 2 dbsnp fix
            dbsnp_ht = dbsnp.ht().select("rsid")
            mt = mt.drop("rsid")
            mt = mt.annotate_rows(**dbsnp_ht[mt.row_key])

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
                        if "hist" not in field:  # some hists are not exported
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

            if args.per_chromosome:
                # Export VCFs by chromosome
                mt = mt.key_rows_by(
                    locus=hl.locus(mt.locus.contig, mt.locus.position),
                    alleles=mt.alleles,
                )
                mt.describe()

                logger.info(f"full mt count: {mt.count()}")
                # NOTE: need to run this on all workers
                rg = get_reference_genome(mt.locus)
                contigs = rg.contigs[:24]  # autosomes + X/Y
                logger.info(f"Contigs: {contigs}")

                for contig in contigs:
                    contig_mt = hl.filter_intervals(
                        mt, [hl.parse_locus_interval(contig)]
                    )
                    logger.info(f"{contig} mt count: {contig_mt.count()}")
                    hl.export_vcf(
                        contig_mt,
                        release_vcf_path(data_source, freeze, contig=contig),
                        metadata=header_dict,
                    )

            if args.h_per_shard:
                # Filter to autosomes + X/Y (remove chrM)
                rg = get_reference_genome(mt.locus)
                contigs = hl.parse_locus_interval(
                    f"{rg.contigs[0]}-{rg.contigs[23]}", reference_genome=rg,
                )
                mt = hl.filter_intervals(mt, [contigs])

                # Export sharded VCF
                mt = mt.naive_coalesce(5000)
                hl.export_vcf(
                    mt,
                    release_vcf_path(data_source, freeze),
                    parallel="header_per_shard",
                    metadata=header_dict,
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
        "--key_by_locus_and_alleles",
        help="Re-key raw MatrixTable by locus and alleles. REQUIRED only for tranche 3/300K dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--repartition",
        help="Repartition raw MatrixTable on read. REQUIRED only for tranche 3/300K dataset",
        action="store_true",
    )
    parser.add_argument(
        "--raw_partitions",
        help="Number of desired partitions for raw MatrixTable. Only used for tranche 3/300K dataset",
        default=30000,
        type=int,
    )
    parser.add_argument(
        "--n_partitions",
        help="Number of desired partitions for output Tables",
        default=10000,
        type=int,
    )
    parser.add_argument(
        "--prepare_internal_mt",
        help="Prepare internal MatrixTable",
        action="store_true",
    )
    parser.add_argument(
        "--prepare_vcf_mt", help="Use release mt to create vcf mt", action="store_true"
    )
    parser.add_argument(
        "--sanity_check", help="Run sanity checks function", action="store_true"
    )
    parser.add_argument(
        "--prepare_browser_ht", help="Prepare sites ht for browser", action="store_true"
    )
    parser.add_argument(
        "--verbose",
        help="Run sanity checks function with verbose output",
        action="store_true",
    )
    parser.add_argument(
        "--prepare_release_vcf", help="Prepare release VCF", action="store_true"
    )
    parser.add_argument(
        "--per_chromosome",
        help="Prepare release VCFs per chromosome",
        action="store_true",
    )
    parser.add_argument(
        "--h_per_shard",
        help="Parallelize release VCF export using header_per_shard",
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
