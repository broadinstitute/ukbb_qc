import argparse
import json
import logging
import pickle
import sys
from typing import Dict, List, Union

import hail as hl

from gnomad.resources.grch37.gnomad import SUBPOPS
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

# Remove BaseQRankSum from site fields (doesn't exist in UKBB 300K)
SITE_FIELDS.remove("BaseQRankSum")

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
KEEP_GNOMAD_POPS = ["afr", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"]
KEEP_GNOMAD_POPS.extend(GNOMAD_NFE_SUBPOPS)
KEEP_GNOMAD_POPS.extend(GNOMAD_EAS_SUBPOPS)


def populate_info_dict(
    subpops: Dict[str, List[str]],
    bin_edges: Dict[str, str],
    age_hist_data: str,
    info_dict: Dict[str, Dict[str, str]] = VCF_INFO_DICT,
    subset_list: List[str] = SUBSET_LIST,
    groups: List[str] = GROUPS,
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
    :param List[str] pop_names: List of sample global population names. Default is POP_NAMES.
    :param List[str] faf_pops: List of faf population names. Default is FAF_POPS.
    :param List[str] gnomad_nfe_subpops: List of nfe subpopulations in gnomAD. Default is GNOMAD_NFE_SUBPOPS.
    :param List[str] gnomad_eas_subpops: List of eas subpopulations in gnomAD. Default is GNOMAD_EAS_SUBPOPS.
    :rtype: Dict[str, Dict[str, str]]
    """
    vcf_info_dict = info_dict

    # Remove unnecessary pop names from pops dict
    for pop in pops.copy():
        if pop not in KEEP_GNOMAD_POPS:
            pops.pop(pop)

    logger.info("Removing gnomAD subpops from UKBB population description dict...")
    ukbb_pops = pops.copy()
    for pop in pops:
        if (pop in gnomad_nfe_subpops) or (pop in gnomad_eas_subpops):
            ukbb_pops.pop(pop)

    # Remove MISSING_REGION_FIELDS from info dict
    for field in MISSING_REGION_FIELDS:
        vcf_info_dict.pop(field, None)

    # Add allele-specific fields to info dict, including AS_VQSLOD and AS_culprit
    temp_AS_fields = AS_FIELDS.copy()
    temp_AS_fields.extend(["AS_culprit", "AS_VQSLOD"])
    vcf_info_dict.update(
        add_as_info_dict(info_dict=vcf_info_dict, as_fields=temp_AS_fields)
    )

    all_ukbb_label_groups = [
        dict(group=["adj"], sex=SEXES_UKBB),
        dict(group=["adj"], pop=ukbb_pops),
        dict(group=["adj"], pop=ukbb_pops, sex=SEXES_UKBB),
    ]
    all_gnomad_label_groups = [
        dict(group=["adj"], sex=SEXES),
        dict(group=["adj"], pop=pops),
        dict(group=["adj"], pop=pops, sex=SEXES),
    ]
    faf_label_groups = [
        dict(group=["adj"]),
        dict(group=["adj"], sex=SEXES_UKBB),
        dict(group=["adj"], pop=faf_pops),
        dict(group=["adj"], pop=faf_pops, sex=SEXES_UKBB),
    ]
    faf_gnomad_label_groups = [
        dict(group=["adj"]),
        dict(group=["adj"], sex=SEXES),
        dict(group=["adj"], pop=faf_pops),
        dict(group=["adj"], pop=faf_pops, sex=SEXES),
    ]

    for subset in subset_list:

        if "gnomad" in subset:
            description_text = " in gnomAD"

            vcf_info_dict.update(
                make_info_dict(
                    prefix=subset,
                    label_groups=dict(group=groups),
                    description_text=description_text,
                )
            )
            vcf_info_dict.update(
                make_info_dict(
                    prefix=subset, popmax=True, description_text=description_text
                )
            )

            for label_group in all_gnomad_label_groups:
                vcf_info_dict.update(
                    make_info_dict(prefix=subset, label_groups=label_group)
                )

            for label_group in faf_gnomad_label_groups:
                vcf_info_dict.update(
                    make_info_dict(prefix=subset, label_groups=label_group, faf=True)
                )

            vcf_info_dict.update(
                make_info_dict(
                    prefix=subset,
                    label_groups=dict(
                        group=["adj"], pop=["nfe"], subpop=gnomad_nfe_subpops
                    ),
                    description_text=description_text,
                )
            )
            vcf_info_dict.update(
                make_info_dict(
                    prefix=subset,
                    label_groups=dict(
                        group=["adj"], pop=["eas"], subpop=gnomad_eas_subpops
                    ),
                    description_text=description_text,
                )
            )
        else:
            vcf_info_dict.update(
                make_info_dict(prefix=subset, label_groups=dict(group=groups))
            )

            for label_group in all_ukbb_label_groups:
                vcf_info_dict.update(
                    make_info_dict(
                        prefix=subset, pop_names=ukbb_pops, label_groups=label_group
                    )
                )

            for label_group in faf_label_groups:
                vcf_info_dict.update(
                    make_info_dict(prefix=subset, label_groups=label_group, faf=True)
                )

            vcf_info_dict.update(
                make_info_dict(
                    prefix=subset,
                    bin_edges=bin_edges,
                    popmax=True,
                    age_hist_data="|".join(str(x) for x in age_hist_data),
                )
            )

    logger.info(
        "Adding UKBB subpops (hybrid pops) to UKBB population description dict..."
    )
    for pop in subpops:
        for cluster in subpops[pop]:
            ukbb_pops[
                f"{cluster}"
            ] = f"{pops[pop]} and hybrid population cluster {cluster}"
    for pop in subpops:
        vcf_info_dict.update(
            make_info_dict(
                prefix="",
                pop_names=ukbb_pops,
                label_groups=dict(group=["adj"], pop=[pop], subpop=subpops[pop]),
            )
        )

    # Add variant quality histograms to info dict
    # for key in vcf_info_dict:
    #    print(key, vcf_info_dict[key])
    # print(bin_edges)
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
        combo_fields = entry[1:] + [entry[0]]
        combo = "_".join(combo_fields)

        if gnomad:

            # Skip all gnomAD subsets
            # Note: this is relevant to the exomes only
            if not genome:
                subset_labels = ("controls", "non")
                if entry[0] in subset_labels:
                    continue

                # Re-create combo to make sure formatting is consistent with other pop lables
                # Manually create combo with "_adj"
                # NOTE: entry[0] in the gnomAD exomes faf meta is "gnomad"
                combo = "_".join(entry[1:] + ["adj"])

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

    try:

        if args.prepare_vcf_mt:
            hybrid_pop_map = args.hybrid_pop_map
            hybrid_pops = [
                pop for sublist in list(hybrid_pop_map.values()) for pop in sublist
            ]

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
                meta_root="meta",
            ).select_entries(*SPARSE_ENTRIES)
            mt = mt.transmute_cols(sex_karyotype=mt.meta.sex_imputation.sex_karyotype)

            logger.info("Filtering to chr20 (for tests only)...")
            # Using filter intervals to keep all the work done by get_ukbb_data
            # (removing sample with withdrawn consent/their ref blocks/variants,
            # also keeping meta col annotations)
            mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chr20")])

            logger.info("Splitting raw MT...")
            mt = hl.experimental.sparse_split_multi(mt)
            mt = mt.select_entries(*ENTRIES)

            logger.info("Reading in release HT and annotating onto raw MT...")
            ht = hl.read_table(release_ht_path(*tranche_data))
            mt = mt.annotate_rows(**ht[mt.row_key])
            mt = mt.annotate_globals(**ht.index_globals())

            logger.info(
                "Checkpointing release MT (chr20 only + sparse) for testing/pharmas..."
            )
            from ukbb_qc.resources.basics import get_release_path

            mt = mt.checkpoint(
                f"{get_release_path(*tranche_data)}/mt/{data_source}.freeze_{freeze}.release.sparse.mt",
                overwrite=args.overwrite,
            )
            logger.info(f"Release MT (sparse; chr20 only) count: {mt.count()}")
            # mt = hl.read_matrix_table(
            #    "gs://broad-ukbb/broad.freeze_6/release/mt/broad.freeze_6.release.sparse.mt"
            # )
            ht = hl.read_table(release_ht_path(*tranche_data))

            logger.info(
                "Dropping cohort frequencies (necessary only for internal use; at last four indices of freq struct)..."
            )
            # print(hl.eval(mt.freq_meta))
            mt = mt.annotate_rows(freq=mt.freq[:-4])
            mt = mt.annotate_globals(freq_meta=mt.freq_meta[:-4])
            # print(hl.eval(mt.freq_meta))

            logger.info("Making histogram bin edges...")
            # NOTE: using release HT here because age histograms aren't necessarily defined
            # in the first row of the raw MT (we may have filtered that row because it was low qual)
            bin_edges = make_hist_bin_edges_expr(ht, prefix="")

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
                region_flag=mt.region_flag.annotate(
                    nonpar=(mt.locus.in_x_nonpar() | mt.locus.in_y_nonpar())
                )
            )

            logger.info("Constructing INFO field")
            # Add variant annotations to INFO field
            # This adds annotations from:
            #   RF struct, VQSR struct, allele_info struct,
            #   info struct (site and allele-specific annotations),
            #   region_flag struct, and
            #   raw_qual_hists/qual_hists structs.

            # TODO: Remove this in the 500K, I modified the variant QC code to create the correct names
            mt = mt.transmute_rows(
                rf=mt.rf.annotate(
                    rf_positive_label=mt.rf.tp,
                    rf_negative_label=mt.rf.fail_hard_filters,
                )
            )
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

            # Reformat vep annotation
            mt = mt.annotate_rows(vep=vep_struct_to_csq(mt.vep))
            mt = mt.annotate_rows(info=mt.info.annotate(vep=mt.vep))

            # Select relevant fields for VCF export
            mt = mt.select_rows("info", "filters", "rsid", "qual")
            new_vcf_info_dict.update(
                {"vep": {"Description": hl.eval(mt.vep_csq_header)}}
            )

            # Make filter dict and add field for MonoAllelic filter
            filter_dict = make_vcf_filter_dict(
                mt.rf_globals.rf_snv_cutoff.min_score,
                mt.rf_globals.rf_indel_cutoff.min_score,
                mt.rf_globals.inbreeding_cutoff,
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

                mt.describe()
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
                mt = mt.naive_coalesce(args.n_shards)

                hl.export_vcf(
                    mt,
                    # release_vcf_path(*tranche_data),
                    "gs://broad-ukbb/broad.freeze_6/temp/chr20_release_vcf.bgz",
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
        "--key_by_locus_and_alleles",
        help="Whether to key raw MT by locus and alleles. REQUIRED only for the 300K tranche",
        action="store_true",
    )
    parser.add_argument(
        "--repartition",
        help="Repartition raw MT on read. Needs to be true for tranche 3/freeze 6/300K.",
        action="store_true",
    )
    parser.add_argument(
        "--raw_partitions",
        help="Number of desired partitions for the raw MT. Necessary only for 300K. Used only if --repartition is also specified",
        default=30000,
        type=int,
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
