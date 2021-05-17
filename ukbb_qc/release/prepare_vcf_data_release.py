import argparse
import logging
import pickle
from typing import Dict, List, Optional, Union

import hail as hl

from gnomad.resources.grch38.gnomad import SEXES
from gnomad.sample_qc.ancestry import POP_NAMES
from gnomad.sample_qc.sex import adjust_sex_ploidy

from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import (
    add_as_info_dict,
    ALLELE_TYPE_FIELDS,
    AS_FIELDS,
    ENTRIES,
    FAF_POPS,
    GROUPS,
    HISTS,
    INFO_DICT,
    make_hist_dict,
    make_info_dict,
    REGION_FLAG_FIELDS,
    RF_FIELDS,
    SITE_FIELDS,
    set_female_y_metrics_to_na,
    SPARSE_ENTRIES,
    VQSR_FIELDS,
)
from gnomad.utils.vcf import SEXES as SEXES_STR
from ukbb_qc.assessment.sanity_checks import (
    sanity_check_release_mt,
    vcf_field_check,
)
from ukbb_qc.resources.basics import (
    get_checkpoint_path,
    get_ukbb_data,
    release_ht_path,
    release_vcf_ht_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.variant_qc import info_ht_path
from ukbb_qc.slack_creds import slack_token
from ukbb_qc.utils.constants import UKBB_POPS
from ukbb_qc.utils.utils import GNOMAD_EAS_SUBPOPS, GNOMAD_NFE_SUBPOPS, make_index_dict


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
VCF_INFO_DICT["outside_capture_region"] = {
    "Description": "Variant falls outside exome capture regions"
}
VCF_INFO_DICT["sibling_singleton"] = {
    "Description": "Variant was a callset-wide doubleton that was present only within a sibling pair"
}

# Add interval QC, capture region to REGION_FLAG_FIELDS and remove decoy
# NOTE: MISSING_REGION_FIELDS could change for 500K if we get hg38 files
INTERVAL_FIELDS = ["fail_interval_qc", "outside_capture_region"]
MISSING_REGION_FIELDS = ["decoy"]
REGION_FLAG_FIELDS = [
    field for field in REGION_FLAG_FIELDS if field not in MISSING_REGION_FIELDS
]
REGION_FLAG_FIELDS.extend(INTERVAL_FIELDS)

# Remove BaseQRankSum from site and allele-specific fields (this is a legacy annotation)
SITE_FIELDS.remove("BaseQRankSum")
AS_FIELDS.remove("AS_BaseQRankSum")

# Add sibling singletons to AS_FIELDS
AS_FIELDS.append("sibling_singleton")

# Make subset list (used in properly filling out VCF header descriptions and naming VCF info fields)
SUBSET_LIST = ["", "gnomad_exomes", "gnomad_genomes"]  # empty for ukbb

# Select populations to keep from the list of population names in POP_NAMES
# This removes pop names we don't want to use in the UKBB release
# (e.g., "uniform", "consanguineous") to reduce clutter
KEEP_POPS = ["afr", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"]

# Separating gnomad exome/genome pops and adding 'ami', 'mid' to gnomAD genomes pops
GNOMAD_GENOMES_POPS = {pop: POP_NAMES[pop] for pop in KEEP_POPS}
GNOMAD_GENOMES_POPS["ami"] = "Amish"
GNOMAD_GENOMES_POPS["mid"] = "Middle Eastern"

# Remove unnecessary pop names from pops dict
# Store this as GNOMAD_EXOMES_POPS
KEEP_POPS.extend(GNOMAD_NFE_SUBPOPS)
KEEP_POPS.extend(GNOMAD_EAS_SUBPOPS)
GNOMAD_EXOMES_POPS = {pop: POP_NAMES[pop] for pop in KEEP_POPS}


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
    gnomad_sexes: List[str] = SEXES_STR,
    ukbb_sexes: List[str] = SEXES,
    gnomad_nfe_subpops: List[str] = GNOMAD_NFE_SUBPOPS,
    gnomad_eas_subpops: List[str] = GNOMAD_EAS_SUBPOPS,
    subpops: Optional[Dict[str, List[str]]] = None,
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

    .. note::
        If `subpops` is specified, `ukbb_pops` MUST include a description for each population in `subpops`.

    :param Dict[str, str] bin_edges: Dictionary of variant annotation histograms and their associated bin edges.
    :param str age_hist_data: Pipe-delimited string of age histograms, from `get_age_distributions`.
    :param Dict[str, Dict[str, str]] info_dict: INFO dict to be populated.
    :param List[str] subset_list: List of sample subsets in dataset. Default is SUBSET_LIST.
    :param List[str] groups: List of sample groups [adj, raw]. Default is GROUPS.
    :param Dict[str, str] ukbb_pops: List of sample global population names for UKBB. Default is UKBB_POPS.
    :param Dict[str, str] gnomad_exomes_pops: List of sample global population names for gnomAD exomes. Default is GNOMAD_EXOMES_POPS.
    :param Dict[str, str] gnomad_genomes_pops: List of sample global population names for gnomAD genomes. Default is GNOMAD_GENOMES_POPS.
    :param List[str] faf_pops: List of faf population names. Default is FAF_POPS.
    :param List[str] gnomad_sexes: gnomAD v2 sample sexes ("male", "female") used in VCF export. Default is SEXES_STR. 
    :param List[str] ukbb_sexes: UKBB, gnomAD v3 sample sexes ("XX", "XY") used in VCF export. Default is SEXES.
    :param List[str] gnomad_nfe_subpops: List of nfe subpopulations in gnomAD. Default is GNOMAD_NFE_SUBPOPS.
    :param List[str] gnomad_eas_subpops: List of eas subpopulations in gnomAD. Default is GNOMAD_EAS_SUBPOPS.
    :param Optional[Dict[str, List[str]]] subpops: Dictionary of global population names (keys)
        and all hybrid population cluster names associated with that global pop (values). 
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
        pops: Union[Dict[str, str], List[str]],
        sexes: List[str],
        group: List[str] = ["adj"],
    ) -> List[Dict[str, List[str]]]:
        """
        Generates list of label group dictionaries needed to populate info dictionary.

        Label dictionaries are passed as input to `make_info_dict`.

        :param Union[Dict[str, str], List[str]] pops: Dict or list of population names.
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

            if "exomes" in subset:
                faf_label_groups = _create_label_groups(
                    pops=faf_pops, sexes=gnomad_sexes
                )
                for label_group in faf_label_groups:
                    vcf_info_dict.update(
                        make_info_dict(
                            prefix=subset,
                            pop_names=gnomad_exomes_pops,
                            label_groups=label_group,
                            faf=True,
                            description_text=description_text,
                        )
                    )

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
                faf_label_groups = _create_label_groups(pops=faf_pops, sexes=ukbb_sexes)
                for label_group in faf_label_groups:
                    vcf_info_dict.update(
                        make_info_dict(
                            prefix=subset,
                            pop_names=gnomad_genomes_pops,
                            label_groups=label_group,
                            faf=True,
                            description_text=description_text,
                        )
                    )
                gnomad_genomes_label_groups = _create_label_groups(
                    pops=gnomad_genomes_pops, sexes=ukbb_sexes
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

            # NOTE: Using `ukbb_pops` here because all frequency calculations were run on pan-ancestry labels this tranche
            # `faf_pops` contains only gnomAD population labels
            faf_label_groups = _create_label_groups(pops=ukbb_pops, sexes=ukbb_sexes)
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

    if subpops:
        logger.info(
            "Adding UKBB subpops (hybrid pops) to UKBB population description dict..."
        )
        for pop in subpops:
            for cluster in subpops[pop]:
                ukbb_pops[
                    f"{cluster}"
                ] = f"{ukbb_pops[pop]} and hybrid population cluster {cluster}"
        for pop in subpops:
            vcf_info_dict.update(
                make_info_dict(
                    prefix="",
                    pop_names=ukbb_pops,
                    label_groups=dict(group=["adj"], pop=[pop], subpop=subpops[pop]),
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

    # Histograms to export are:
    # gq_hist_alt, gq_hist_all, dp_hist_alt, dp_hist_all, ab_hist_alt
    # We previously dropped:
    # _n_smaller for all hists
    # _bin_edges for all hists
    # _n_larger for all hists EXCEPT DP hists
    for hist in HISTS:
        for prefix in ["qual_hists", "raw_qual_hists"]:
            hist_name = hist
            if "raw" in prefix:
                hist_name = f"{hist}_raw"

            hist_dict = {
                f"{hist_name}_bin_freq": hl.delimit(
                    t[prefix][hist].bin_freq, delimiter="|"
                ),
            }
            vcf_info_dict.update(hist_dict)

            if "dp" in hist_name:
                vcf_info_dict.update(
                    {f"{hist_name}_n_larger": t[prefix][hist].n_larger},
                )
    return vcf_info_dict


def unfurl_nested_annotations(
    t: Union[hl.MatrixTable, hl.Table],
    gnomad: bool,
    genome: bool,
    pops: List[str],
    subpops: Optional[List[str]] = None,
) -> Dict[str, hl.expr.Expression]:
    """
    Create dictionary keyed by the variant annotation labels to be extracted from variant annotation arrays, where the values
    of the dictionary are Hail Expressions describing how to access the corresponding values.

    :param Table/MatrixTable t: Table/MatrixTable containing the nested variant annotation arrays to be unfurled.
    :param bool gnomad: Whether the annotations are from gnomAD.
    :param bool genome: Whether the annotations are from genome data (relevant only to gnomAD data).
    :param List[str] pops: List of global populations in frequency array. 
    :param List[str] subpops: List of all UKBB subpops (possible hybrid population cluster names). Default is None.
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
        if data_type != "genomes":
            subpops = [GNOMAD_NFE_SUBPOPS + GNOMAD_EAS_SUBPOPS]

        freq_idx = make_index_dict(
            t=t,
            freq_meta_str=f"{gnomad_prefix}_freq_meta",
            pops=pops,
            subpops=subpops,
            data_type=data_type,
        )

    else:
        faf = "faf"
        freq = "freq"
        faf_idx = make_index_dict(
            t=t, freq_meta_str="faf_meta", pops=pops, subpops=subpops
        )
        popmax = "popmax"
        freq_idx = make_index_dict(
            t=t, freq_meta_str="freq_meta", pops=pops, subpops=subpops
        )

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
                # NOTE: faf format in v3.1 changed; no longer has `meta` field
                prefix = f"{gnomad_prefix}_"
                combo_dict = {
                    f"{prefix}faf95_{combo}": hl.or_missing(
                        hl.set(hl.eval(t.gnomad_genomes_faf_meta[i].values()))
                        == set(combo_fields),
                        t[faf][i].faf95,
                    ),
                    f"{prefix}faf99_{combo}": hl.or_missing(
                        hl.set(hl.eval(t.gnomad_genomes_faf_meta[i].values()))
                        == set(combo_fields),
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
    # We previously dropped:
    # age_hist_hom_bin_edges, age_hist_het_bin_edges
    if not gnomad:
        age_hist_dict = {
            "age_hist_het_bin_freq": hl.delimit(t.age_hist_het.bin_freq, delimiter="|"),
            "age_hist_het_n_smaller": t.age_hist_het.n_smaller,
            "age_hist_het_n_larger": t.age_hist_het.n_larger,
            "age_hist_hom_bin_freq": hl.delimit(t.age_hist_hom.bin_freq, delimiter="|"),
            "age_hist_hom_n_smaller": t.age_hist_hom.n_smaller,
            "age_hist_hom_n_larger": t.age_hist_hom.n_larger,
        }
        expr_dict.update(age_hist_dict)

    return expr_dict


def main(args):

    hl.init(log="/vcf_release_patch.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)

    try:
        logger.info("Reading header dict from pickle...")
        # NOTE: I moved this file out of the temp bucket for permanent storage
        # We deleted the temp bucket to save on storage, and I wanted to keep this file
        with hl.hadoop_open(
            "gs://broad-ukbb/broad.freeze_7/release/vcf/header_dict.pickle", "rb"
        ) as p:
            header_dict = pickle.load(p)

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
            meta_root="meta",
        ).select_entries(*SPARSE_ENTRIES)

        logger.info("Loading het non ref sites to fix...")
        sites_ht = hl.read_matrix_table(args.het_non_ref).rows()

        logger.info("Densifying to het non ref sites only...")
        from gnomad.utils.sparse_mt import densify_sites
        from ukbb_qc.resources.basics import (
            last_END_positions_ht_path,
            get_release_path,
        )

        # NOTE: Set semi_join_rows to True here because setting it to False causes jobs to crash
        # NOTE: This densify required a highmem-32 master node
        # https://github.com/broadinstitute/gnomad_methods/blob/master/gnomad/utils/sparse_mt.py#L88
        mt = densify_sites(
            mt,
            sites_ht,
            hl.read_table(last_END_positions_ht_path(freeze)),
            semi_join_rows=True,
        )
        mt = mt.checkpoint(
            get_checkpoint_path(*tranche_data, name="het_non_ref_dense", mt=True)
        )

        logger.info("Adding het_non_ref annotation...")
        # Adding a Boolean for whether a sample had a heterozygous non-reference genotype
        # Need to add this prior to splitting MT to make sure these genotypes
        # are not adjusted by the homalt hotfix downstream
        mt = mt.annotate_entries(het_non_ref=mt.LGT.is_het_non_ref())

        logger.info("Splitting densified MT...")
        mt = hl.split_multi_hts(mt)
        mt = mt.select_entries(*ENTRIES)

        # NOTE: densify_sites operates using only the locus (not the alleles)
        logger.info("Filtering only to the impacted variants (filtering on alleles)...")
        logger.info("Removing low QUAL variants and * alleles...")
        info_ht = hl.read_table(info_ht_path(data_source, freeze))
        mt = mt.filter_rows(
            (~info_ht[mt.row_key].AS_lowqual)
            & ((hl.len(mt.alleles) > 1) & (mt.alleles[1] != "*"))
            & (hl.is_defined(sites_ht[mt.row_key]))
        )

        logger.info("Adjusting sex ploidy...")
        mt = mt.annotate_cols(sex_karyotype=mt.meta.sex_imputation.sex_karyotype)
        mt = adjust_sex_ploidy(mt, mt.sex_karyotype, male_str="XY", female_str="XX")

        # Temporary hotfix for depletion of homozygous alternate genotypes
        logger.info(
            "Setting het genotypes at sites with >1% AF and > 0.9 AB to homalt..."
        )
        # NOTE: Reading release HT here because frequency annotation was updated
        # Only `join_freq` and release HT have updated frequency annotation
        freq_ht = (
            hl.read_table(release_ht_path(*tranche_data))
            .select_globals()
            .select("freq")
        )
        freq_ht = freq_ht.select(AF=freq_ht.freq[0].AF)
        mt = mt.annotate_entries(
            GT=hl.if_else(
                mt.GT.is_het()
                # Skip adjusting genotypes if sample originally had a het nonref genotype
                & ~mt.het_non_ref
                & (freq_ht[mt.row_key].AF > 0.01)
                & (mt.AD[1] / mt.DP > 0.9),
                hl.call(1, 1),
                mt.GT,
            )
        )

        logger.info("Annotating release MT with HT annotations...")
        logger.info("Dropping frequency globals...")
        ht = hl.read_table(release_vcf_ht_path(*tranche_data)).drop("freq_meta")
        mt = mt.annotate_rows(**ht[mt.row_key])
        mt = mt.annotate_globals(**ht.index_globals())

        logger.info("Re-calculating frequency...")
        from gnomad.utils.annotations import (
            annotate_freq,
            faf_expr,
            pop_max_expr,
        )

        logger.info("Filtering related samples and their variants...")
        mt_filt = mt.select_rows().select_globals()
        mt_filt = mt_filt.filter_cols(~mt_filt.meta.sample_filters.related)
        mt_filt = mt_filt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
        mt_filt = annotate_freq(
            mt_filt,
            sex_expr=mt_filt.meta.sex_imputation.sex_karyotype,
            pop_expr=mt_filt.meta.pan_ancestry_meta.pop,
            additional_strata_expr=None,
        )
        faf, faf_meta = faf_expr(mt_filt.freq, mt_filt.freq_meta, mt_filt.locus)
        popmax = pop_max_expr(mt_filt.freq, mt_filt.freq_meta)
        logger.info("Calculating faf and popmax...")
        freq_ht = mt_filt.annotate_rows(faf=faf, popmax=popmax).rows()
        # TODO: Will need to update release HT frequencies at these het nonref sites
        freq_ht = freq_ht.checkpoint(
            get_checkpoint_path(*tranche_data, name="ukb_freq_het_non_ref")
        )

        logger.info("Annotating new frequency struct and globals onto MT...")
        mt = mt.annotate_rows(**freq_ht[mt.row_key])
        mt = mt.annotate_globals(**freq_ht.index_globals())

        logger.info(
            "Unfurling new UKBB frequency annotations and adding to INFO field..."
        )
        mt = mt.annotate_rows(
            info=mt.info.annotate(
                **unfurl_nested_annotations(
                    mt, gnomad=False, genome=False, pops=UKBB_POPS,
                )
            )
        )

        logger.info("Dropping column annotations...")
        mt = mt.select_cols()

        # NOTE: `qual` annotation is actually `QUALapprox` annotation in 455k tranche
        # Need to convert this field to a float because `export_vcf` won't export this field
        # if the type isn't float64
        mt = mt.annotate_rows(qual=hl.float(mt.qual))

        # NOTE: Fixing chrY metrics here for 455k tranche
        # because fix to `set_female_y_metrics_to_na` was pushed
        # after VCF MT generated
        # https://github.com/broadinstitute/gnomad_methods/pull/347
        mt = mt.annotate_rows(info=mt.info.annotate(**set_female_y_metrics_to_na(mt)))

        logger.info("Sanity checking release MT...")
        sanity_check_release_mt(
            mt,
            SUBSET_LIST,
            ukbb_pops=UKBB_POPS,
            gnomad_exomes_pops=GNOMAD_EXOMES_POPS,
            gnomad_genomes_pops=GNOMAD_GENOMES_POPS,
            missingness_threshold=0.5,
            verbose=args.verbose,
        )

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

        logger.info("Exporting VCF...")
        hl.export_vcf(
            dataset=mt,
            output=f"{get_release_path(*tranche_data)}/vcf/sharded_vcf/broad.freeze_7.patch.bgz",
            parallel=None,
            metadata=header_dict,
            # NOTE: I also moved this file out of temp when we were deleting temporary files
            append_to_header="gs://broad-ukbb/broad.freeze_7/release/vcf/append_to_vcf_header.tsv",
            tabix=True,
        )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(f"gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/")


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
        "--het_non_ref", help="Path to MT with het non ref sites to be fixed",
    )
    parser.add_argument(
        "--verbose",
        help="Run sanity checks function with verbose output",
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
