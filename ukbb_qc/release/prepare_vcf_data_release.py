import argparse
import logging
import pickle
import sys
from typing import Dict, List, Union

import hail as hl

from gnomad.resources.grch37.gnomad import SUBPOPS
from gnomad.resources.grch38.reference_data import dbsnp, seg_dup_intervals
from gnomad.resources.resource_utils import DataException
from gnomad.sample_qc.ancestry import POP_NAMES
from gnomad.sample_qc.sex import adjust_sex_ploidy
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import vep_struct_to_csq
from gnomad.utils.vcf import (
    add_as_info_dict,
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
    array_sample_map_ht_path,
    get_checkpoint_path,
    logging_path,
    raw_mt_path,
    release_header_path,
    release_ht_path,
    release_mt_path,
    release_vcf_path,
    release_vcf_ht_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import meta_ht_path
from ukbb_qc.resources.variant_qc import (
    info_ht_path,
    NA12878,
    SYNDIP,
    var_annotations_ht_path,
)
from ukbb_qc.slack_creds import slack_token
from ukbb_qc.utils.constants import SEXES_UKBB
from ukbb_qc.utils.utils import make_index_dict


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("vcf_release")
logger.setLevel(logging.INFO)


# Add END to entries
ENTRIES.append("END")

# Add capture region and sibling singletons to vcf_info_dict
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

# Remove BaseQRankSum from site and AS fields (doesn't exist in UKBB 300K)
SITE_FIELDS.remove("BaseQRankSum")
AS_FIELDS.remove("AS_BaseQRankSum")

# Add sibling singletons to AS_FIELDS
AS_FIELDS.append("sibling_singleton")

# Make subset list (used in properly filling out VCF header descriptions and naming VCF info fields)
SUBSET_LIST = ["", "gnomad_exomes", "gnomad_genomes"]  # empty for ukbb

# Get gnomAD subpop names
GNOMAD_NFE_SUBPOPS = list(map(lambda x: x.lower(), SUBPOPS["NFE"]))
GNOMAD_EAS_SUBPOPS = list(map(lambda x: x.lower(), SUBPOPS["EAS"]))

# Select populations to keep from the list of population names in POP_NAMES
# This removes pop names we don't want to use in the UKBB release
# (e.g., "uniform", "consanguineous") to reduce clutter
KEEP_POPS = ["afr", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"]
UKBB_POPS = {pop: POP_NAMES[pop] for pop in KEEP_POPS}

KEEP_POPS.extend(GNOMAD_NFE_SUBPOPS)
KEEP_POPS.extend(GNOMAD_EAS_SUBPOPS)

# Remove unnecessary pop names from pops dict
POPS = {pop: POP_NAMES[pop] for pop in KEEP_POPS}

# Separating gnomad exome/genome pops and adding 'ami' to gnomAD genomes pops
GNOMAD_EXOMES_POPS = POPS.copy()
GNOMAD_GENOMES_POPS = POPS.copy()
GNOMAD_GENOMES_POPS["ami"] = "Amish"


def populate_info_dict(
    bin_edges: Dict[str, str],
    age_hist_data: str,
    info_dict: Dict[str, Dict[str, str]] = VCF_INFO_DICT,
    subset_list: List[str] = SUBSET_LIST,
    groups: List[str] = GROUPS,
    ukbb_pops: Dict[str, str] = UKBB_POPS,
    gnomad_exomes_pops: Dict[str, str] = GNOMAD_EXOMES_POPS,
    gnomad_genomes_pops: Dict[str, str] = GNOMAD_GENOMES_POPS,
    faf_pops: List[str] = FAF_POPS,
    gnomad_sexes: List[str] = SEXES,
    ukbb_sexes: List[str] = SEXES_UKBB,
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

    :param Dict[str, str] bin_edges: Dictionary of variant annotation histograms and their associated bin edges.
    :param str age_hist_data: Pipe-delimited string of age histograms, from `get_age_distributions`.
    :param Dict[str, Dict[str, str]] info_dict: INFO dict to be populated.
    :param List[str] subset_list: List of sample subsets in dataset. Default is SUBSET_LIST.
    :param List[str] groups: List of sample groups [adj, raw]. Default is GROUPS.
    :param Dict[str, str] ukbb_pops: List of sample global population names for UKBB. Default is UKBB_POPS.
    :param Dict[str, str] gnomad_exomes_pops: List of sample global population names for gnomAD exomes. Default is GNOMAD_EXOMES_POPS.
    :param Dict[str, str] gnomad_genomes_pops: List of sample global population names for gnomAD genomes. Default is GNOMAD_GENOMES_POPS.
    :param List[str] faf_pops: List of faf population names. Default is FAF_POPS.
    :param List[str] gnomad_sexes: gnomAD sample sexes used in VCF export. Default is SEXES. 
    :param List[str] ukbb_sexes: UKBB sample sexes used in VCF export. Default is SEXES_UKBB.
    :param List[str] gnomad_nfe_subpops: List of nfe subpopulations in gnomAD. Default is GNOMAD_NFE_SUBPOPS.
    :param List[str] gnomad_eas_subpops: List of eas subpopulations in gnomAD. Default is GNOMAD_EAS_SUBPOPS.
    :rtype: Dict[str, Dict[str, str]]
    """
    vcf_info_dict = info_dict

    # Remove MISSING_REGION_FIELDS from info dict
    for field in MISSING_REGION_FIELDS:
        vcf_info_dict.pop(field, None)

    # Add allele-specific fields to info dict, including AS_VQSLOD and AS_culprit
    # NOTE: need to think about how to resolve AS VQSR fields to avoid having to make temp_AS_fields variable in the future
    temp_AS_fields = AS_FIELDS.copy()
    temp_AS_fields.extend(["AS_culprit", "AS_VQSLOD"])
    vcf_info_dict.update(
        add_as_info_dict(info_dict=vcf_info_dict, as_fields=temp_AS_fields)
    )

    def _create_label_groups(
        pops: Dict[str, str], sexes: List[str], group: List[str] = ["adj"],
    ) -> List[Dict[str, List[str]]]:
        """
        Generates list of label group dictionaries needed to populate info dictionary.

        Label dictionaries are passed as input to `make_info_dict`.

        :param Dict[str, str] pops: List of population names.
        :param List[str] sexes: List of sample sexes.
        :param List[str] group: List of data types (adj, raw). Default is ["adj"].
        :return: List of label group dictionaries.
        :rtype: List[Dict[str, List[str]]]
        """
        return [
            dict(group=groups),  # this is to capture raw fields
            dict(group=group, sex=sexes),
            dict(group=group, pop=pops),
            dict(group=group, pop=pops, sex=sexes),
        ]

    for subset in subset_list:

        if "gnomad" in subset:
            description_text = " in gnomAD"

            # faf labels are the same for exomes/genomes
            faf_label_groups = _create_label_groups(pops=faf_pops, sexes=gnomad_sexes)
            for label_group in faf_label_groups:
                vcf_info_dict.update(
                    make_info_dict(
                        prefix=subset,
                        pop_names=gnomad_exomes_pops,
                        label_groups=label_group,
                        faf=True,
                    )
                )

            if "exomes" in subset:
                gnomad_exomes_label_groups = _create_label_groups(
                    pops=gnomad_exomes_pops, sexes=gnomad_sexes
                )
                for label_group in gnomad_exomes_label_groups:
                    vcf_info_dict.update(
                        make_info_dict(
                            prefix=subset,
                            pop_names=gnomad_exomes_pops,
                            label_groups=label_group,
                            description_text=description_text,
                        )
                    )

                # Add popmax to info dict
                vcf_info_dict.update(
                    make_info_dict(
                        prefix=subset,
                        pop_names=gnomad_exomes_pops,
                        popmax=True,
                        description_text=description_text,
                    )
                )

                # Add gnomAD exome subpops to info dict
                vcf_info_dict.update(
                    make_info_dict(
                        prefix=subset,
                        pop_names=gnomad_exomes_pops,
                        label_groups=dict(
                            group=["adj"], pop=["nfe"], subpop=gnomad_nfe_subpops
                        ),
                        description_text=description_text,
                    )
                )
                vcf_info_dict.update(
                    make_info_dict(
                        prefix=subset,
                        pop_names=gnomad_exomes_pops,
                        label_groups=dict(
                            group=["adj"], pop=["eas"], subpop=gnomad_eas_subpops
                        ),
                        description_text=description_text,
                    )
                )

            else:
                gnomad_genomes_label_groups = _create_label_groups(
                    pops=gnomad_genomes_pops, sexes=gnomad_sexes
                )
                for label_group in gnomad_genomes_label_groups:
                    vcf_info_dict.update(
                        make_info_dict(
                            prefix=subset,
                            pop_names=gnomad_genomes_pops,
                            label_groups=label_group,
                            description_text=description_text,
                        )
                    )

                # Add popmax to info dict
                vcf_info_dict.update(
                    make_info_dict(
                        prefix=subset,
                        pop_names=gnomad_genomes_pops,
                        popmax=True,
                        description_text=description_text,
                    )
                )

        else:
            ukbb_label_groups = _create_label_groups(pops=ukbb_pops, sexes=ukbb_sexes)
            for label_group in ukbb_label_groups:
                vcf_info_dict.update(
                    make_info_dict(
                        prefix=subset, pop_names=ukbb_pops, label_groups=label_group,
                    )
                )

            faf_label_groups = _create_label_groups(pops=faf_pops, sexes=ukbb_sexes)
            for label_group in faf_label_groups:
                vcf_info_dict.update(
                    make_info_dict(
                        prefix=subset,
                        pop_names=ukbb_pops,
                        label_groups=label_group,
                        faf=True,
                    )
                )

            vcf_info_dict.update(
                make_info_dict(
                    prefix=subset,
                    bin_edges=bin_edges,
                    popmax=True,
                    age_hist_data="|".join(str(x) for x in age_hist_data),
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
    t: Union[hl.MatrixTable, hl.Table], gnomad: bool, genome: bool, pops: List[str],
) -> Dict[str, hl.expr.Expression]:
    """
    Create dictionary keyed by the variant annotation labels to be extracted from variant annotation arrays, where the values
    of the dictionary are Hail Expressions describing how to access the corresponding values.

    :param Table/MatrixTable t: Table/MatrixTable containing the nested variant annotation arrays to be unfurled.
    :param bool gnomad: Whether the annotations are from gnomAD.
    :param bool genome: Whether the annotations are from genome data (relevant only to gnomAD data).
    :param List[str] pops: List of global populations in frequency array. 
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
            t=t, freq_meta_str=f"{gnomad_prefix}_freq_meta", pops=pops,
        )

    else:
        faf = "faf"
        freq = "freq"
        faf_idx = make_index_dict(t=t, freq_meta_str="faf_meta", pops=pops)
        popmax = "popmax"
        freq_idx = make_index_dict(t=t, freq_meta_str="freq_meta", pops=pops)

    # Unfurl freq index dict
    # Cycles through each key and index (e.g., k=adj_afr, i=31)
    for k, i in freq_idx.items():
        prefix = ""
        if gnomad:
            prefix = f"{gnomad_prefix}_"

        # Set combination to key
        # e.g., set entry of 'afr_adj' to combo
        combo = k
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

        if gnomad:

            # Create combo_fields
            # NOTE: gnomad needs to be handled separately because gnomad faf index dict has different format
            # gnomad faf index dict has keys like 'adj_afr', but UKBB will index with keys like 'afr_adj'
            combo_fields = entry[1:] + [entry[0]]
            combo = "_".join(combo_fields)

            # Skip all gnomAD subsets
            # Note: this is relevant to the exomes only
            if not genome:
                subset_labels = ("controls", "non")
                if entry[0] in subset_labels:
                    continue

                # Re-create combo to make sure formatting is consistent with other pop labels
                # Manually create combo with "_adj"
                # NOTE: entry[0] in the gnomAD exomes faf meta is "gnomad",
                # but the value in the faf meta is still "adj"
                combo_fields = entry[1:] + ["adj"]
                combo = "_".join(combo_fields)

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
            # Set combo to equal entry
            combo_fields = entry
            combo = k

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

    try:
        if args.prepare_vcf_annotations:
            logger.info("Reading in release HT...")
            # NOTE: new_freq is an annotation created when fixing the missing 300K frequencies
            # it is a temp annotation and should have been dropped before writing HT
            # Also drop VEP here (this annotation uses the older VEP version with a bug),
            # so we want to add a new annotation from VEP version 101
            ht = hl.read_table(release_ht_path(*tranche_data)).drop("new_freq", "vep")

            logger.info(
                "Dropping frequencies stratified with subpops (using UKBB hybrid pops as the subpop)..."
            )
            logger.info(
                "Also dropping cohort frequencies (necessary only for internal use; at last four indices of freq struct)..."
            )
            # Keep:
            # Indices 0 and 1: adj, raw freqs
            # Indices 2-9: freq stratified by gnomAD pop label
            # Indices 10-11: freq stratified by sex
            # Indices 12-27: freq stratified by pop and sex
            # Cohort freq notes:
            # Cohort freq has 4 entries in freq and freq meta:
            # cohort (adj), cohort (raw), cohort (XX), and cohort (XY)
            ht = ht.annotate_rows(freq=ht.freq[:28])
            ht = ht.annotate_globals(freq_meta=ht.freq_meta[:28])

            logger.info("Making histogram bin edges...")
            # NOTE: using release HT here because age histograms aren't necessarily defined
            # in the first row of the raw MT (we may have filtered that row because it was low qual)
            bin_edges = make_hist_bin_edges_expr(ht, prefix="")

            logger.info("Getting age hist data...")
            age_hist_data = hl.eval(ht.age_distribution)

            # This was removed from the 455K but is necessary for the 300K
            from ukbb_qc.resources.basics import vqsr_sites_path
            from gnomad.utils.sparse_mt import split_info_annotation

            logger.info("Pulling AS_VarDP from VQSR sites HT (300K fix)...")
            vqsr_sites_ht = hl.read_table(vqsr_sites_path(*tranche_data))
            vqsr_sites_ht = hl.split_multi(vqsr_sites_ht)
            vqsr_sites_ht = vqsr_sites_ht.annotate(
                info=vqsr_sites_ht.info.annotate(
                    **split_info_annotation(vqsr_sites_ht.info, vqsr_sites_ht.a_index),
                ),
            )
            vqsr_sites_ht = vqsr_sites_ht.transmute(
                AS_VarDP=vqsr_sites_ht.info.AS_VarDP
            )
            ht = ht.annotate(
                info=ht.info.annotate(AS_VarDP=vqsr_sites_ht[ht.key].AS_VarDP)
            )

            logger.info("Making INFO dict for VCF...")
            vcf_info_dict = populate_info_dict(
                bin_edges=bin_edges, age_hist_data=age_hist_data,
            )

            # Add interval QC parameters to INFO dict
            pct_samples = hl.eval(ht.rf_globals.interval_qc_cutoffs.pct_samples) * 100
            autosome_cov = hl.eval(ht.rf_globals.interval_qc_cutoffs.autosome_cov)
            allosome_cov = hl.eval(ht.rf_globals.interval_qc_cutoffs.xy_cov)
            vcf_info_dict["fail_interval_qc"] = {
                "Description": f"Variant falls within a region where less than {pct_samples}% of samples had a mean coverage of {autosome_cov}X on autosomes and {allosome_cov}X on sex chromosomes"
            }

            # Adjust keys to remove adj tags before exporting to VCF
            new_vcf_info_dict = {
                i.replace("_adj", ""): j for i, j in vcf_info_dict.items()
            }

            # NOTE: We did not have a segdup file when preparing the 300K release,
            # which is why this annotation is missing
            logger.info("Updating region flag expr to include segdup intervals...")
            segdup_ht = seg_dup_intervals.ht()
            ht = ht.annotate_rows(
                region_flag=ht.region_flag.annotate(
                    # Also add non-PAR annotation
                    nonpar=(ht.locus.in_x_nonpar() | ht.locus.in_y_nonpar()),
                    segdup=hl.is_defined(segdup_ht[ht.locus]),
                )
            )

            # This was removed in the 455K tranche because the variant QC code was updated to have the correct names
            ht = ht.transmute(
                rf=ht.rf.annotate(
                    rf_positive_label=ht.rf.tp,
                    rf_negative_label=ht.rf.fail_hard_filters,
                )
            )
            ht = ht.annotate_rows(info=hl.struct(**make_info_expr(ht)))

            logger.info("Reformatting rsid...")
            # Update rsid annotation to be a set and not a string
            dbsnp_ht = dbsnp.ht().select("rsid")
            ht = ht.annotate(rsid=hl.str(";").join(dbsnp_ht[ht.key].rsid))

            logger.info("Constructing INFO field")
            # Add variant annotations to INFO field
            # This adds annotations from:
            #   RF struct, VQSR struct, allele_info struct,
            #   info struct (site and allele-specific annotations),
            #   region_flag struct, and
            #   raw_qual_hists/qual_hists structs.

            # Unfurl nested UKBB frequency annotations and add to INFO field
            ht = ht.annotate(
                info=ht.info.annotate(
                    **unfurl_nested_annotations(
                        ht, gnomad=False, genome=False, pops=UKBB_POPS,
                    )
                )
            )
            # Unfurl nested gnomAD genome frequency annotations and add to info field
            ht = ht.annotate(
                info=ht.info.annotate(
                    **unfurl_nested_annotations(
                        ht, gnomad=True, genome=True, pops=GNOMAD_GENOMES_POPS,
                    )
                )
            )
            # Unfurl nested gnomAD exome frequency annotations and add to info field
            ht = ht.annotate(
                info=ht.info.annotate(
                    **unfurl_nested_annotations(
                        ht, gnomad=True, genome=False, pops=GNOMAD_EXOMES_POPS,
                    )
                )
            )
            ht = ht.annotate(**set_female_y_metrics_to_na(ht))

            logger.info("Reading in VEP HT...")
            vep_ht = hl.read_table(var_annotations_ht_path("vep", *tranche_data))
            vep_ht = vep_ht.transmute(vep=vep_ht.vep.drop("colocated_variants"))

            # Reformat vep annotation
            vep_ht = vep_ht.annotate(vep=vep_struct_to_csq(vep_ht.vep))
            ht = ht.annotate(info=ht.info.annotate(vep=vep_ht[ht.key].vep))

            logger.info("Selecting fields and writing VCF HT...")
            ht = ht.select("info", "filters", "rsid", "qual")
            ht.write(release_vcf_ht_path(*tranche_data), args.overwrite)
            new_vcf_info_dict.update(
                {"vep": {"Description": hl.eval(mt.vep_csq_header)}}
            )

            # Make filter dict and add field for MonoAllelic filter
            filter_dict = make_vcf_filter_dict(
                hl.eval(mt.rf_globals.rf_snv_cutoff.min_score),
                hl.eval(mt.rf_globals.rf_indel_cutoff.min_score),
                hl.eval(mt.rf_globals.inbreeding_cutoff),
            )
            filter_dict["MonoAllelic"] = {
                "Description": "Samples are all homozygous reference or all homozygous alternate for the variant"
            }
            header_dict = {
                "info": new_vcf_info_dict,
                "filter": filter_dict,
                "format": FORMAT_DICT,
            }

            logger.info("Saving header dict to pickle...")
            with hl.hadoop_open(release_header_path(*tranche_data), "wb") as p:
                pickle.dump(header_dict, p, protocol=pickle.HIGHEST_PROTOCOL)

        if args.prepare_vcf_mt:
            logger.info("Preparing VCF MT...")
            logger.info("Getting raw MT and dropping all unnecessary entries...")

            # NOTE: Reading in raw MT from raw_mt_path to do custom sample filtering
            mt = hl.read_matrix_table(raw_mt_path(*tranche_data)).select_entries(
                *SPARSE_ENTRIES
            )
            # Add column index (to remove duplicate sample IDs)
            mt = mt.add_col_index()
            # Key MT by locus and alleles
            mt = hl.MatrixTable(
                hl.ir.MatrixKeyRowsBy(mt._mir, ["locus", "alleles"], is_sorted=True)
            )
            # Annotate meta HT information onto MT cols
            meta_ht = hl.read_table(meta_ht_path(*tranche_data))
            mt = mt.annotate_cols(**{"meta": meta_ht[mt.s]})
            mt = mt.transmute_cols(sex_karyotype=mt.meta.sex_imputation.sex_karyotype)

            logger.info("Removing duplicate sample IDs...")
            # UKB_4048554_0301608642 and UKB_1223807_0330880742 are both present in the 300K MT twice
            # Both versions of UKB_1223807_0330880742 (indices 202837 and 205124) are the same, so
            # code removes the version of the sample with a higher column index (205124)
            # UKB_4048554_0301608642 has column indices 81090 and 262414
            # The version of the sample with the higher index corresponds to the more recent version of the sample,
            # so code removes the sample with the lower column index
            original_sample_count = mt.count_cols()
            logger.info("Original sample count: %i", original_sample_count)
            mt = mt.annotate_cols(new_s=hl.format("%s_%s", mt.s, mt.col_idx))
            remove_ids = hl.literal(
                ["UKB_1223807_0330880742_205124", "UKB_4048554_0301608642_81090"]
            )
            mt = mt.filter_cols(~remove_ids.contains(mt.new_s)).drop("new_s", "col_idx")
            filtered_sample_count = mt.count_cols()
            logger.info(
                "Sample count after removing duplicate sample IDs: %i",
                filtered_sample_count,
            )
            if filtered_sample_count != original_sample_count - 2:
                raise DataException(
                    "Expected to remove two duplicate sample IDs. Please double check and rerun!"
                )

            logger.info(
                "Removing samples with withdrawn consent, control samples, and samples with undefined UKBB batch numbers..."
            )
            # File downloaded on 8/16/21
            withdrawn_ht = hl.import_table(
                "gs://broad-ukbb/resources/withdrawn_consents/w26041_20210809.csv",
                no_header=True,
            ).key_by("f0")
            sample_map_ht = hl.read_table(array_sample_map_ht_path(freeze))
            sample_map_ht = sample_map_ht.annotate(
                withdrawn_consent=hl.is_defined(
                    withdrawn_ht[sample_map_ht.ukbb_app_26041_id]
                )
            )
            withdrawn_ids = sample_map_ht.aggregate(
                hl.agg.count_where(sample_map_ht.withdrawn_consent)
            )
            logger.info(
                "Total number of IDs with withdrawn consents in sample map ht: %i",
                withdrawn_ids,
            )

            all_sample_count = mt.count_cols()
            logger.info(
                "MT sample count before removing samples with withdrawn consent, control samples, and samples with defined UKBB batch: %i",
                all_sample_count,
            )
            controls = hl.literal([NA12878, SYNDIP])
            mt = mt.filter_cols(
                hl.is_defined(sample_map_ht[mt.col_key].batch)
                & ~sample_map_ht[mt.col_key].withdrawn_consent
                & ~controls.contains(mt.col_key)
            )
            final_sample_count = mt.count_cols()
            logger.info(
                "MT sample count after removing withdrawn and control samples: %i",
                final_sample_count,
            )

            if (final_sample_count - all_sample_count) < withdrawn_ids:
                raise DataException(
                    "Number of removed samples is less than total number of samples with withdrawn consents in sample map HT. Please double check and rerun!"
                )

            logger.info("Removing chrM...")
            mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chrM")], keep=False)

            if args.test:
                logger.info("Filtering to chr20 and chrX (for tests only)...")
                # Using filter intervals to keep all the work done by get_ukbb_data
                # (removing sample with withdrawn consent/their ref blocks/variants,
                # also keeping meta col annotations)
                # Using chr20 to test a small autosome and chrX to test a sex chromosome
                # Some annotations (like FAF) are 100% missing on autosomes
                mt = hl.filter_intervals(
                    mt,
                    [hl.parse_locus_interval("chr20"), hl.parse_locus_interval("chrX")],
                )

            logger.info("Splitting raw MT...")
            mt = hl.experimental.sparse_split_multi(mt)
            mt = mt.select_entries(*ENTRIES)

            logger.info("Reading in VCF HT and annotating onto raw MT...")
            ht = hl.read_table(release_vcf_ht_path(*tranche_data))
            mt = mt.annotate_rows(**ht[mt.row_key])
            mt = mt.annotate_globals(**ht.index_globals())

            logger.info(
                "Coalescing MT to number of shards desired for output (necessary to get start/end positions for each shard correctly)..."
            )
            mt = mt.naive_coalesce(args.n_shards)
            mt.write(release_mt_path(*tranche_data), args.overwrite)

        if args.sanity_check:
            mt = hl.read_matrix_table(release_mt_path(*tranche_data))
            # NOTE: removing lowqual and star alleles here to avoid having additional failed missingness checks
            info_ht = hl.read_table(info_ht_path(data_source, freeze))
            mt = mt.filter_rows(
                (~info_ht[mt.row_key].AS_lowqual)
                & ((hl.len(mt.alleles) > 1) & (mt.alleles[1] != "*"))
            )
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

            # Drop unnecessary histograms
            # TODO: figure out which hists we want to export and only create those for 500K
            drop_hists = (
                [x + "_n_smaller" for x in HISTS]
                + [x + "_bin_edges" for x in HISTS]
                + [x + "_n_larger" for x in HISTS if "dp_" not in x]
            )
            drop_hists.extend(
                [x + "_raw_n_smaller" for x in HISTS]
                + [x + "_raw_bin_edges" for x in HISTS]
                + [x + "_raw_n_larger" for x in HISTS if "dp_" not in x]
                + ["age_hist_hom_bin_edges", "age_hist_het_bin_edges"]
            )
            mt = mt.annotate_rows(info=mt.info.drop(*drop_hists))

            # Reformat names to remove "adj" pre-export
            # e.g, renaming "AC_adj" to "AC"
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
            mt = mt.annotate_rows(
                info=mt.info.select(
                    "AC",
                    "AN",
                    "AF",
                    "rf_tp_probability",
                    *mt.info.drop("AC", "AN", "AF", "rf_tp_probability"),
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
                    # Checked with Hail team about the fastest way to filter to a contig
                    # This method shouldn't be any faster than `filter_intervals`: the same amount of data is read in both cases
                    # `_calculate_new_partitions` might give us more parallelism downstream
                    # Decided to stick with `_calculate_new_partitions` method because it felt much faster on the 300K tranche
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

                    logger.info("Adjusting sex ploidy...")
                    mt = adjust_sex_ploidy(
                        mt, mt.sex_karyotype, male_str="XY", female_str="XX"
                    )

                    hl.export_vcf(
                        mt.select_cols(),
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

                logger.info("Adjusting sex ploidy...")
                mt = adjust_sex_ploidy(
                    mt, mt.sex_karyotype, male_str="XY", female_str="XX"
                )
                mt = mt.select_cols()

                logger.info("Exporting VCF...")
                hl.export_vcf(
                    mt,
                    release_vcf_path(*tranche_data),
                    parallel="header_per_shard",
                    metadata=header_dict,
                    append_to_header="gs://broad-ukbb/broad.freeze_7/release/vcf/append_to_vcf_header.tsv",
                    tabix=True,
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
        "--key_by_locus_and_alleles",
        help="Whether to key raw MT by locus and alleles. REQUIRED only for tranche 3/freeze 6/300K",
        action="store_true",
    )
    parser.add_argument(
        "--repartition",
        help="Repartition raw MT on read. Needs to be true for tranche 3/freeze 6/300K.",
        action="store_true",
    )
    parser.add_argument(
        "--raw_partitions",
        help="Number of desired partitions for the raw MT. Necessary only for tranche 3/freeze 6/300K. Used only if --repartition is also specified",
        default=30000,
        type=int,
    )
    parser.add_argument(
        "--test",
        help="Create release files using only chr20 and chrX for testing purposes",
        action="store_true",
    )
    parser.add_argument(
        "--prepare_vcf_annotations",
        help="Use release HT to reformat VCF annotations",
        action="store_true",
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
    export_opts = parser.add_mutually_exclusive_group()
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
