import argparse
import logging
import pickle
from typing import Dict, List, Union

import hail as hl

from gnomad.resources.grch38.gnomad import SEXES
from gnomad.sample_qc.sex import adjust_sex_ploidy
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import vep_struct_to_csq
from gnomad.utils.vcf import (
    add_as_info_dict,
    ALLELE_TYPE_FIELDS,
    AS_FIELDS,
    ENTRIES,
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
    SPARSE_ENTRIES,
    VQSR_FIELDS,
)

from ukbb_qc.assessment.sanity_checks import (
    sanity_check_release_mt,
    vcf_field_check,
)
from ukbb_qc.resources.basics import (
    append_to_vcf_header_path,
    get_checkpoint_path,
    get_ukbb_data,
    logging_path,
    release_header_path,
    release_ht_path,
    release_vcf_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.variant_qc import info_ht_path
from ukbb_qc.slack_creds import slack_token
from ukbb_qc.utils.constants import UKBB_POPS
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
VCF_INFO_DICT["outside_capture_region"] = {
    "Description": "Variant falls outside exome capture regions"
}
VCF_INFO_DICT["sibling_singleton"] = {
    "Description": "Variant was a callset-wide doubleton that was present only within a sibling pair"
}

# Add AS_SB_TABLE and AS_QUALapprox to vcf_info_dict
# NOTE: This is necessary because they were added by default into `INFO_DICT` in a more
# current version of the gnomad repo than what is necessary to export UKBB data
VCF_INFO_DICT["AS_SB_TABLE"] = {
    "Number": ".",
    "Description": "Allele-specific forward/reverse read counts for strand bias tests",
}
VCF_INFO_DICT["AS_QUALapprox"] = {
    "Number": "1",
    "Description": "Sum of PL[0] values; used to approximate the QUAL score",
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
# Also remove SB from site fields since we do not release it
SITE_FIELDS.remove("BaseQRankSum")
SITE_FIELDS.remove("SB")
AS_FIELDS.remove("AS_BaseQRankSum")

# Add sibling singletons to AS_FIELDS
AS_FIELDS.append("sibling_singleton")


def populate_info_dict(
    bin_edges: Dict[str, str],
    age_hist_data: str,
    info_dict: Dict[str, Dict[str, str]] = VCF_INFO_DICT,
    groups: List[str] = GROUPS,
    ukbb_pops: Dict[str, str] = UKBB_POPS,
    ukbb_sexes: List[str] = SEXES,
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
    :param List[str] ukbb_sexes: UKBB, gnomAD v3 sample sexes ("XX", "XY") used in VCF export. Default is SEXES.
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

    # Set prefix to empty string (do not want to add extra subset description to autopopulated text for export)
    prefix = ""
    ukbb_label_groups = _create_label_groups(pops=ukbb_pops, sexes=ukbb_sexes)
    for label_group in ukbb_label_groups:
        vcf_info_dict.update(
            make_info_dict(
                prefix=prefix, pop_names=ukbb_pops, label_groups=label_group,
            )
        )

    # NOTE: Using `ukbb_pops` here because all frequency calculations were run on pan-ancestry labels
    # `faf_pops` contains only gnomAD population labels
    faf_label_groups = _create_label_groups(pops=ukbb_pops, sexes=ukbb_sexes)
    for label_group in faf_label_groups:
        vcf_info_dict.update(
            make_info_dict(
                prefix=prefix, pop_names=ukbb_pops, label_groups=label_group, faf=True,
            )
        )

    vcf_info_dict.update(
        make_info_dict(
            prefix=prefix,
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
    t: Union[hl.MatrixTable, hl.Table], pops: List[str],
) -> Dict[str, hl.expr.Expression]:
    """
    Create dictionary keyed by the variant annotation labels to be extracted from variant annotation arrays, where the values
    of the dictionary are Hail Expressions describing how to access the corresponding values.

    :param Table/MatrixTable t: Table/MatrixTable containing the nested variant annotation arrays to be unfurled.
    :param List[str] pops: List of global populations in frequency array. 
    :return: Dictionary containing variant annotations and their corresponding values.
    :rtype: Dict[str, hl.expr.Expression]
    """
    expr_dict = dict()

    # Set variables to locate necessary fields, compute freq index dicts, and compute faf index dict for UKBB
    faf = "faf"
    freq = "freq"
    prefix = ""
    faf_idx = make_index_dict(t=t, freq_meta_str="faf_meta", pops=pops)
    popmax = "popmax"
    freq_idx = make_index_dict(t=t, freq_meta_str="freq_meta", pops=pops)

    # Unfurl freq index dict
    # Cycles through each key and index (e.g., k=adj_afr, i=31)
    for k, i in freq_idx.items():
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

    hl.init(log="/vcf_release.log", default_reference="GRCh38")
    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)

    # NOTE: This isn't generally recommended, but these export steps appear to be particularly
    # slow if this flag isn't set
    logger.info("Setting hail flag (to try to speed up computations)...")
    hl._set_flags(no_whole_stage_codegen="1")

    try:

        if args.prepare_vcf_annotations:
            logger.info("Starting VCF process...")
            logger.info("Reading in release HT...")
            # Also drop all gnomAD annotations
            ht = hl.read_table(release_ht_path(*tranche_data)).drop(
                "gnomad_exomes_freq",
                "gnomad_exomes_popmax",
                "gnomad_exomes_faf",
                "gnomad_genomes_freq",
                "gnomad_genomes_popmax",
                "gnomad_genomes_faf",
                "gnomad_exomes_freq_meta",
                "gnomad_exomes_popmax_index_dict",
                "gnomad_exomes_faf_index_dict",
                "gnomad_exomes_freq_index_dict",
                "gnomad_genomes_freq_meta",
                "gnomad_genomes_faf_index_dict",
                "gnomad_genomes_faf_meta",
                "gnomad_genomes_freq_index_dict",
            )

            logger.info(
                "Dropping cohort frequencies (necessary only for internal use)..."
            )
            # Cohort freq has 22 entries in freq and freq meta:
            # cohort (adj), cohort (raw), cohort (pop), cohort (sex), cohort (pop and sex)
            # Two sexes: XX, XY
            # Six pops (pan-ancestry labels): CSA, MID, AFR, EAS, AMR, EUR
            ht = ht.annotate(freq=ht.freq[:22])
            ht = ht.annotate_globals(freq_meta=ht.freq_meta[:22])

            # Add AS_SB_TABLE, AS_QUALapprox here because it is missing from the release HT
            info_ht = hl.read_table(info_ht_path(data_source, freeze))
            ht = ht.annotate(
                info=ht.info.annotate(
                    AS_SB_TABLE=hl.array(
                        [
                            info_ht[ht.key].info.AS_SB_TABLE[:2],
                            info_ht[ht.key].info.AS_SB_TABLE[2:],
                        ]
                    ),
                    AS_QUALapprox=info_ht[ht.key].info.AS_QUALapprox,
                )
            )

            # Reformat AS_SB_TABLE for export
            # NOTE: Copied this function here because this doesn't exist in the version of the gnomad repo required for
            # UKB validity checks and export
            def _get_pipe_expr(
                array_expr: hl.expr.ArrayExpression,
            ) -> hl.expr.StringExpression:
                return hl.delimit(
                    array_expr.map(lambda x: hl.or_else(hl.str(x), "")), "|"
                )

            ht = ht.annotate(
                info=ht.info.annotate(
                    AS_SB_TABLE=_get_pipe_expr(
                        ht.info.AS_SB_TABLE.map(lambda x: hl.delimit(x, ","))
                    )
                )
            )

            logger.info("Reading in release patch frequencies...")
            patch_ht = (
                hl.read_matrix_table(
                    get_checkpoint_path(
                        *tranche_data,
                        name="release_patch_sites_dense_annot_no_hyphen.mt",
                        mt=True,
                    ),
                )
                .rows()
                .select_globals()
                .select("freq")
            )

            logger.info("Adding patch frequencies...")
            ht = ht.annotate(patch_freq=patch_ht[ht.key].freq)
            ht = ht.transmute(freq=hl.coalesce(ht.patch_freq, ht.freq))

            logger.info("Making histogram bin edges...")
            # NOTE: using release HT here because age histograms aren't necessarily defined
            # in the first row of the raw MT (we may have filtered that row because it was low qual)
            bin_edges = make_hist_bin_edges_expr(ht, prefix="")

            logger.info("Getting age hist data...")
            age_hist_data = hl.eval(ht.age_distribution)

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

            # Add non-PAR annotation
            ht = ht.annotate(
                region_flag=ht.region_flag.annotate(
                    nonpar=(ht.locus.in_x_nonpar() | ht.locus.in_y_nonpar())
                )
            )

            logger.info("Constructing INFO field")
            # Add variant annotations to INFO field
            # This adds annotations from:
            #   RF struct, VQSR struct, allele_info struct,
            #   info struct (site and allele-specific annotations),
            #   region_flag struct, and
            #   raw_qual_hists/qual_hists structs.

            ht = ht.annotate(info=hl.struct(**make_info_expr(ht)))

            # Unfurl nested UKBB frequency annotations and add to INFO field
            ht = ht.annotate(
                info=ht.info.annotate(**unfurl_nested_annotations(ht, pops=UKBB_POPS))
            )
            ht = ht.annotate(**set_female_y_metrics_to_na(ht))

            # Reformat vep annotation
            ht = ht.annotate(vep=vep_struct_to_csq(ht.vep))
            ht = ht.annotate(info=ht.info.annotate(vep=ht.vep))
            new_vcf_info_dict.update(
                {"vep": {"Description": hl.eval(ht.vep_csq_header)}}
            )

            # NOTE: rsid on release HT is a string but should be a set -- need to fix before releasing
            # Re-add rsid annotation to HT and reformat for export
            logger.info("Reformatting rsid...")
            # dbsnp might have multiple identifiers for one variant
            # thus, rsid is a set annotation, starting with version b154 for dbsnp resource:
            # https://github.com/broadinstitute/gnomad_methods/blob/master/gnomad/resources/grch38/reference_data.py#L136
            # `export_vcf` expects this field to be a string, and vcf specs
            # say this field may be delimited by a semi-colon:
            # https://samtools.github.io/hts-specs/VCFv4.2.pdf
            # dbsnp ht doesn't work in this version of gnomad methods
            dbsnp_ht = hl.read_table(
                "gs://gnomad-public-requester-pays/resources/grch38/dbsnp/dbsnp_b154_grch38_all_20200514.ht"
            ).select("rsid")
            ht = ht.annotate(rsid=dbsnp_ht[ht.key].rsid)
            ht = ht.annotate(rsid=hl.str(";").join(ht.rsid))

            logger.info(
                "Selecting relevant fields for VCF export and checkpointing HT..."
            )
            ht = ht.select("info", "filters", "rsid", "qual")
            ht.write(
                "gs://broad-ukbb/broad.freeze_7/release/ht/broad.freeze_7.release.vcf.ukb_official_export.ht",
                overwrite=args.overwrite,
            )

            # Make filter dict and add field for MonoAllelic filter
            filter_dict = make_vcf_filter_dict(
                hl.eval(ht.rf_globals.rf_snv_cutoff.min_score),
                hl.eval(ht.rf_globals.rf_indel_cutoff.min_score),
                hl.eval(ht.rf_globals.inbreeding_cutoff),
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
            logger.info("Getting raw MT and dropping all unnecessary entries...")
            # NOTE: reading in raw MatrixTable to be able to return all samples/variants
            mt = get_ukbb_data(
                *tranche_data,
                key_by_locus_and_alleles=args.key_by_locus_and_alleles,
                split=False,
                raw=True,
                ukbb_samples_only=True,
                repartition=args.repartition,
                n_partitions=args.raw_partitions,
                meta_root="meta",
            ).select_entries(*SPARSE_ENTRIES)
            mt = mt.annotate_cols(sex_karyotype=mt.meta.sex_imputation.sex_karyotype)

            logger.info("Removing chrM...")
            mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chrM")], keep=False)

            if args.test:
                logger.info("Filtering to chr20 and chrX (for tests only)...")
                # Using filter intervals to keep all the work done by get_ukbb_data
                # (removing sample with withdrawn consent/their ref blocks/variants,
                # also keeping meta col annotations)
                # Using chr20 to test a small autosome and chrX to test a sex chromosome
                # Some annotations (like FAF) are 100% missing on autosomes
                mt_chr20 = hl.filter_intervals(mt, [hl.parse_locus_interval("chr20")])
                mt_chr20 = mt_chr20._filter_partitions(range(2))

                mt_chrx = hl.filter_intervals(mt, [hl.parse_locus_interval("chrX")])
                mt_chrx = mt_chrx.filter_rows(mt_chrx.locus.in_x_nonpar())
                mt_chrx = mt_chrx._filter_partitions(range(2))
                mt = mt_chr20.union_rows(mt_chrx)

            logger.info("Adding het_non_ref annotation...")
            # Adding a Boolean for whether a sample had a heterozygous non-reference genotype
            # Need to add this prior to splitting MT to make sure these genotypes
            # are not adjusted by the homalt hotfix downstream
            mt = mt.annotate_entries(het_non_ref=mt.LGT.is_het_non_ref())

            # Add het_non_ref to ENTRIES (otherwise annotation gets accidentally dropped here)
            ENTRIES.append("het_non_ref")

            logger.info("Splitting raw MT...")
            mt = hl.experimental.sparse_split_multi(mt)
            mt = mt.select_entries(*ENTRIES)

            # Temporary hotfix for depletion of homozygous alternate genotypes
            logger.info(
                "Setting het genotypes at sites with >1% AF and > 0.9 AB to homalt..."
            )
            # NOTE: Reading release HT here because frequency annotation was updated
            # Only release HT has updated frequency annotation
            freq_ht = (
                hl.read_table(release_ht_path(data_source, freeze=6))
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

            logger.info("Changing sample IDs to UKBB IDs...")
            mt = mt.key_cols_by(s=mt.meta.ukbb_meta.ukbb_app_26041_id)

            logger.info("Annotating release MT with HT annotations...")
            ht = hl.read_table(
                "gs://broad-ukbb/broad.freeze_7/release/ht/broad.freeze_7.release.vcf.ukb_official_export.ht"
            )
            mt = mt.annotate_rows(**ht[mt.row_key])
            mt = mt.annotate_globals(**ht.index_globals())
            mt.write(
                "gs://broad-ukbb/broad.freeze_7/release/ht/broad.freeze_7.release.vcf.ukb_official_export.mt",
                args.overwrite,
            )

        if args.sanity_check:
            mt = hl.read_matrix_table(
                "gs://broad-ukbb/broad.freeze_7/release/ht/broad.freeze_7.release.vcf.ukb_official_export.mt"
            )

            # NOTE: removing lowqual and star alleles here to avoid having additional failed missingness checks
            info_ht = hl.read_table(info_ht_path(data_source, freeze))
            mt = mt.filter_rows(
                (~info_ht[mt.row_key].AS_lowqual)
                & ((hl.len(mt.alleles) > 1) & (mt.alleles[1] != "*"))
            )

            # NOTE: Fixing chrY metrics here because they were not correctly annotated into the info struct on the VCF HT
            mt = mt.annotate_rows(
                info=mt.info.annotate(**set_female_y_metrics_to_na(mt))
            )
            sanity_check_release_mt(
                mt,
                ukbb_pops=UKBB_POPS,
                missingness_threshold=0.5,
                verbose=args.verbose,
            )

        if args.prepare_release_vcf:

            logger.warning(
                "VCF export will densify! Make sure you have an autoscaling cluster."
            )
            mt = hl.read_matrix_table(
                "gs://broad-ukbb/broad.freeze_7/release/ht/broad.freeze_7.release.vcf.ukb_official_export.mt",
                _n_partitions=args.n_shards,
            )

            # NOTE: Fixing chrY metrics here because the code above previously annotated the fixed metrics onto the VCF HT
            # but added the metrics as top level annotations rather than adding them into the info struct
            # Line 488 should have been:
            # ht = ht.annotate(info=ht.info.annotate(**set_female_y_metrics_to_na(ht))
            mt = mt.annotate_rows(
                info=mt.info.annotate(**set_female_y_metrics_to_na(mt))
            )

            # NOTE: `qual` annotation is actually `QUALapprox` annotation in 455k tranche
            # Need to convert this field to a float because `export_vcf` won't export this field
            # if the type isn't float64
            mt = mt.annotate_rows(qual=hl.float(mt.qual))

            if args.test:
                mt = mt.filter_partitions(range(2))

            logger.info("Reading header dict from pickle...")
            with hl.hadoop_open(release_header_path(*tranche_data), "rb") as p:
                header_dict = pickle.load(p)

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
                    "nhomalt",
                    "rf_tp_probability",
                    *mt.info.drop("AC", "AN", "AF", "nhomalt", "rf_tp_probability"),
                )
            )

            ht = mt.rows().checkpoint(
                get_checkpoint_path(*tranche_data, name="flat_vcf_ready", mt=False),
                overwrite=args.overwrite,
            )

            # Export VCFs per chromosome
            rg = get_reference_genome(mt.locus)
            contigs = rg.contigs[:24]  # autosomes + X/Y
            logger.info(f"Contigs: {contigs}")

            for contig in contigs:
                # Read in MT and filter to contig
                mt = hl.read_matrix_table(
                    "gs://broad-ukbb/broad.freeze_7/release/ht/broad.freeze_7.release.vcf.ukb_official_export.mt"
                )
                mt = hl.filter_intervals(mt, [hl.parse_locus_interval(contig)])
                mt = mt.annotate_rows(**ht[mt.row_key])

                logger.info("Densifying and exporting VCF...")
                mt = hl.experimental.densify(mt)
                # Drop END and het non ref to avoid exporting
                mt = mt.drop("END", "het_non_ref")

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

                logger.info("Adjusting partitions...")
                mt = mt.naive_coalesce(args.n_shards)
                ht = mt.rows()
                # Unkey HT to avoid this error with map_partitions:
                # ValueError: Table._map_partitions must preserve key fields
                ht = ht.key_by()

                hl.export_vcf(
                    mt,
                    release_vcf_path(*tranche_data, contig=contig),
                    metadata=header_dict,
                    append_to_header=append_to_vcf_header_path(*tranche_data),
                    tabix=True,
                )

                logger.info("Getting start and stops per shard...")

                def part_min_and_max(part):
                    keys = part.map(lambda x: x.select("locus", "alleles"))
                    return hl.struct(start=keys[0], end=keys[-1])

                print(
                    ht._map_partitions(
                        lambda p: hl.array([part_min_and_max(p)])
                    ).collect()
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
        "--prepare_vcf_mt", help="Use release MT to create VCF MT", action="store_true"
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
    parser.add_argument(
        "--n_shards",
        help="Desired number of shards for output VCF. Can also be used to repartition raw MT on read.",
        type=int,
        default=500,
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
