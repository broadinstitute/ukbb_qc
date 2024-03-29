import json
import logging
from typing import Dict, List, Union

import hail as hl

from gnomad.resources.grch37.gnomad import SUBPOPS
from gnomad.resources.grch38.gnomad import SEXES
from gnomad.resources.grch38.reference_data import lcr_intervals
from gnomad.resources.resource_utils import DataException
from gnomad.sample_qc.relatedness import UNRELATED
from gnomad.utils.annotations import (
    age_hists_expr,
    qual_hist_expr,
)
from gnomad.utils.file_utils import file_exists
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.vcf import (
    GROUPS,
    make_label_combos,
)
from gnomad.utils.vcf import SEXES as SEXES_STR
from gnomad.variant_qc.random_forest import get_rf_runs

from ukbb_qc.load_data.utils import import_phenotype_ht
from ukbb_qc.resources.basics import (
    array_sample_map_ht_path,
    get_ukbb_data,
    phenotype_ht_path,
    raw_mt_path,
)
from ukbb_qc.resources.sample_qc import (
    f_stat_sites_path,
    interval_qc_path,
    hard_filters_ht_path,
    qc_mt_path,
    qc_sites_path,
)
from ukbb_qc.resources.variant_qc import vqsr_run_path


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("utils")
logger.setLevel(logging.INFO)


GNOMAD_NFE_SUBPOPS = list(map(lambda x: x.lower(), SUBPOPS["NFE"]))
"""
Dictionary of gnomAD v2 nfe subpop labels and description. 

Used in VCF export.
"""

GNOMAD_EAS_SUBPOPS = list(map(lambda x: x.lower(), SUBPOPS["EAS"]))
"""
Dictionary of gnomAD v2 eas subpop labels and description. 

Used in VCF export.
"""


# Generic
def compare_row_counts(ht1: hl.Table, ht2: hl.Table) -> bool:
    """
    Checks if row counts in two Tables are the same

    :param Table ht1: First Table to be checked
    :param Table ht2: Second Table to be checked
    :return: Whether the row counts are the same
    :rtype: bool
    """
    r_count1 = ht1.count()
    r_count2 = ht2.count()
    logger.info(f"{r_count1} rows in left table; {r_count2} rows in right table")
    return r_count1 == r_count2


def join_tables(
    left_ht: hl.Table, left_key: str, right_ht: hl.Table, right_key: str, join_type: str
) -> hl.Table:
    """
    Joins left and right tables using specified keys and join types and returns result.
    
    Also prints warning if sample counts are not the same.

    :param Table left_ht: Left Table to be joined
    :param str left_key: Key of left Table
    :param Table right_ht: Right Table to be joined
    :param str right_key: Key of right Table
    :param str join_type: Type of join
    :return: Table with annotations
    :rtype: Table
    """
    sample_count_match = compare_row_counts(left_ht, right_ht)
    if not sample_count_match:
        logger.warning("Sample counts in left and right tables do not match!")
    return left_ht.key_by(left_key).join(right_ht.key_by(right_key), how=join_type)


# Sample resources
def get_sample_ids(ids_file: str, header: bool = False) -> List[str]:
    """
    Open file with all sample IDs and stores IDs in a list.

    :param str ids_file: Path to file containing sample IDs with variants to extract for readviz.
    :param bool header: Whether IDs file has a header line. Default is False.
    :return: List of sample IDs
    :rtype: List[str]
    """
    sample_ids = []
    with hl.hadoop_open(ids_file) as i:
        if header:
            header = i.readline()
        for line in i:
            sample_ids.append(line.strip())
    return sample_ids


def remove_hard_filter_samples(
    data_source: str,
    freeze: int,
    t: Union[hl.MatrixTable, hl.Table],
    filter_rows: bool = True,
    non_refs_only: bool = True,
    gt_field: str = "LGT",
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Removes samples that failed hard filters from MatrixTable or Table. 

    Assumes MatrixTable's col key/Table's row key is sample.
    If non_refs_only is False, assumes input MatrixTable is sparse!

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param MatrixTable/Table t: Input MatrixTable or Table
    :param bool filter_rows: Whether to filter rows. Relevant only if input is a MatrixTable. Default is True
    :param bool non_refs_only: Whether to filter to non_reference sites only. 
        Relevant only if filter_rows is True and input is a MatrixTable.
    :param str gt_field: Field containing genotype. Relevant only if filter_rows is True.
        Default is LGT.
    :return: MatrixTable or Table with samples (and their variants, if non_refs_only is set) removed
    :rtype: MatrixTable or Table
    """
    ht = hl.read_table(hard_filters_ht_path(data_source, freeze)).select("hard_filters")
    hf_samples = ht.aggregate(hl.agg.count_where(ht.hard_filters.hard_filtered))

    logger.info(f"Removing {hf_samples} samples that failed hard filters")
    ht = ht.filter(~ht.hard_filters.hard_filtered)
    if isinstance(t, hl.MatrixTable):
        t = t.filter_cols(hl.is_defined(ht[t.col_key]))
        if filter_rows:
            if non_refs_only:
                t = t.filter_rows(hl.agg.any(t[gt_field].is_non_ref()))
            else:
                t = t.filter_rows(
                    hl.agg.any(t[gt_field].is_non_ref() | hl.is_defined(t.END))
                )
    else:
        t = t.filter(hl.is_defined(ht[t.key]))
    return t


# Variant resources
def annotate_interval_qc_filter(
    data_source: str,
    freeze: int,
    t: Union[hl.MatrixTable, hl.Table],
    cov_filter_field: str = "pct_samples_20x",
    xy_cov_filter_field: str = "pct_samples_10x",
    pct_samples: float = 0.85,
    autosomes_only: bool = False,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Removes poorly covered intervals from a MatrixTable or Table.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param MatrixTable/able t: Input MatrixTable or Table
    :param float cov_filter_field: Coverage field to use for filtering of autosomes/sex chr par/female X
    :param float xy_cov_filter_field: Coverage field to use for filtering of male X and Y
    :param float pct_samples: Percent of samples with coverage greater than `cov_filter_field` over the interval for filtering
    :param bool autosomes_only: Whether to annotate autosomes only. Default is False
    :return: MatrixTable or Table with variants within intervals excluded by coverage cutoffs removed
    :rtype: MatrixTable or Table
    """
    logger.warning(
        "This function looks only at variant position (doesn't check END fields)!\n"
        "That means this function only works on dense or non_ref only data."
    )
    interval_qc_ht = hl.read_table(interval_qc_path(data_source, freeze, "autosomes"))
    good_intervals_ht = interval_qc_ht.filter(
        interval_qc_ht[cov_filter_field] > pct_samples
    )

    if not autosomes_only:
        interval_qc_sex_ht = hl.read_table(
            interval_qc_path(data_source, freeze, "sex_chr")
        )

        good_sex_intervals_ht = interval_qc_sex_ht.filter(
            (
                interval_qc_sex_ht.interval.start.in_autosome_or_par()
                & (interval_qc_sex_ht[cov_filter_field]["XX"] > pct_samples)
                & (interval_qc_sex_ht[cov_filter_field]["XY"] > pct_samples)
            )
            | (
                (interval_qc_sex_ht.interval.start.contig == "chrX")
                & (interval_qc_sex_ht.interval.start.in_x_nonpar())
                & (interval_qc_sex_ht[cov_filter_field]["XX"] > pct_samples)
                & (interval_qc_sex_ht[xy_cov_filter_field]["XY"] > pct_samples)
            )
            | (
                (interval_qc_sex_ht.interval.start.contig == "chrY")
                & (interval_qc_sex_ht.interval.start.in_y_nonpar())
                & (interval_qc_sex_ht[xy_cov_filter_field]["XY"] > pct_samples)
            )
        )
        good_intervals_ht = good_intervals_ht.select().union(
            good_sex_intervals_ht.select()
        )

    if isinstance(t, hl.MatrixTable):
        t = t.annotate_rows(interval_qc_pass=hl.is_defined(good_intervals_ht[t.locus]))
    else:
        t = t.annotate(interval_qc_pass=hl.is_defined(good_intervals_ht[t.locus]))

    return t


def vqsr_run_check(data_source: str, freeze: int, vqsr_type: str) -> None:
    """
    Checks the run data log `vqsr_run_path` for the presence of `vqsr_type`.

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param int vqsr_type: Name of VQSR model
    :return: None
    """
    vqsr_runs = get_rf_runs(vqsr_run_path(data_source, freeze))
    if vqsr_type not in vqsr_runs:
        vqsr_runs_print = ""
        for run_hash, run_data in vqsr_runs.items():
            vqsr_runs_print = vqsr_runs_print + f"\n=== {run_hash} ==="
            vqsr_runs_print = vqsr_runs_print + json.dumps(
                run_data, sort_keys=True, indent=4, separators=(",", ": ")
            )

        raise DataException(
            f"Requested VQSR data doesn't exist: {vqsr_type}. Current loaded models: {vqsr_runs_print}"
        )


# Sites resources
def get_sites(
    af_cutoff: float,
    call_rate_cutoff: float,
    greater: bool = True,
    pass_interval_qc: bool = True,
    autosomes_only: bool = True,
) -> hl.Table:
    """
    Returns a Table of sites that meet af_cutoff from the 200K (tranche 2/freeze 5) callset. 
    
    :param float af_cutoff: Desired AF cutoff
    :param float call_rate_cutoff: Desired callrate cutoff
    :param bool greater: Whether to return variants with AFs higher than the cutoff.
    :param bool pass_interval_qc: Whether to only return sites that pass interval QC. 
    :param bool autosomes_only: Whether to return autosomal sites only. Also sets interval QC flag.
    :return: None
    :rtype: None
    """
    data_source = "broad"
    freeze = 5
    mt = get_ukbb_data(data_source, freeze, split=True)

    if autosomes_only:
        mt = filter_to_autosomes(mt)

    if pass_interval_qc:
        mt = annotate_interval_qc_filter("broad", 5, mt, autosomes_only=autosomes_only)
        mt = mt.filter_rows(mt.interval_qc_pass)

    # Add callstats (for AF) and filter on AF cutoff
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.call_rate > call_rate_cutoff)
    if greater:
        mt = mt.filter_rows(mt.variant_qc.AF[1] >= af_cutoff)
    else:
        mt = mt.filter_rows(mt.variant_qc.AF[1] < af_cutoff)

    ht = mt.rows()
    ht = ht.transmute(AF=ht.variant_qc.AF[1])
    return ht


def calculate_fstat_sites(call_rate: float = 0.99, af_cutoff: float = 0.001) -> None:
    """
    Writes a Table with high callrate, common, biallelic SNPs in regions that pass interval QC on chromosome X.

    This Table is designed to be used as an interval filter in sex imputation.
    NOTE: This function generates sites only for tranche 2/freeze 5, which is the last dataset that contains AF.

    :param float call_rate: Call rate cutoff used to filter sites
    :param float af_cutoff: Alternate allele frequency cutoff used to filter sites
    :return: None
    :rtype: None
    """
    data_source = "broad"
    freeze = 5
    mt = hl.read_matrix_table(raw_mt_path(data_source, freeze, densified=True))
    mt = mt.key_rows_by("locus", "alleles")
    mt = hl.filter_intervals(
        mt,
        [
            hl.parse_locus_interval(
                x_contig, reference_genome=get_reference_genome(mt.locus).name
            )
            for x_contig in get_reference_genome(mt.locus).x_contigs
        ],
    )
    mt = mt.filter_rows(
        (hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1])
    )
    mt = annotate_interval_qc_filter("broad", 5, mt, autosomes_only=False)
    mt = mt.filter_rows(mt.interval_qc_pass)
    mt = mt.transmute_entries(GT=hl.experimental.lgt_to_gt(mt.LGT, mt.LA))
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(
        (mt.variant_qc.call_rate > call_rate) & (mt.variant_qc.AF[1] > af_cutoff)
    )
    ht = mt.rows()
    ht = ht.transmute(AF=ht.variant_qc.AF[1])
    ht.write(f_stat_sites_path())


def get_qc_mt_sites() -> None:
    """
    Writes a Table with sites to use in QC MatrixTable generation. 

    Table includes the fields 'locus', 'alleles', and 'info'. 
    NOTE: This function generates sites based on the tranche 2/freeze 5 QC MatrixTable.

    :return: None
    :rtype: None
    """
    logger.info("Preparing to make QC MT sites HT")
    data_source = "broad"
    freeze = 5
    mt = hl.read_matrix_table(qc_mt_path(data_source, freeze, ld_pruned=True))
    ht = mt.rows()

    logger.info("Removing LCR intervals from QC sites")
    lcr = lcr_intervals.ht()
    ht = ht.filter(hl.is_missing(lcr[ht.locus]))
    ht.write(qc_sites_path())


# Sample-related resources
def get_phenotype_field(freeze: int, field: str) -> hl.Table:
    """
    Pull phenotype information using input field from UKBB phenotype file

    :param int freeze: One of the data freezes
    :param str field: Field containing desired phenotype information
    :return: Table with phenotype information per sample
    :rtype: hl.Table
    """
    if not file_exists(phenotype_ht_path()):
        import_phenotype_ht()
    ht = hl.read_table(phenotype_ht_path()).select(field)

    # Re-key phenotype table to UKBB ID using array sample map table and return
    sample_map_ht = hl.read_table(array_sample_map_ht_path(freeze))
    sample_map_ht = sample_map_ht.key_by("ukbb_app_26041_id")
    ht = ht.key_by(s=sample_map_ht[ht.key].s)
    return ht.select(field)


def get_age_ht(freeze: int, field: str = "f.21022.0.0") -> hl.Table:
    """
    Pull age information from UKBB phenotype file

    :param int freeze: One of the data freezes
    :param str field: Field containing age information. Default is f.21022.0.0.
    :return: Table with age at recruitment per sample
    :rtype: hl.Table
    """
    age_ht = get_phenotype_field(freeze, field)
    age_ht = age_ht.rename({field: "age"})
    return age_ht.select("age")


def get_array_sex_ht(freeze: int, field: str = "f.22001.0.0") -> hl.Table:
    """
    Pulls genetic sex information from UKBB phenotype file (Field: 22001)

    :param int freeze: One of the data freezes
    :param str field: Field containing age information. Default is f.22001.0.0.
    :return: Table with genetic sex per sample
    :rtype: Table
    """
    sex_ht = get_phenotype_field(freeze, field)
    return sex_ht.annotate(
        array_sex=hl.if_else(sex_ht["f.22001.0.0"] == 0, "XX", "XY")
    ).select("array_sex")


def get_relatedness_set_ht(relatedness_ht: hl.Table) -> hl.Table:
    """
    Parses relatedness Table to get every relationship (except UNRELATED) per sample.

    Returns Table keyed by sample with all sample relationships in a set.

    :param Table relatedness_ht: Table with inferred relationship information output by pc_relate. 
        Keyed by sample pair (i, j).
    :return: Table keyed by sample (s) with all relationships annotated as a set.
    :rtype: hl.Table
    """
    relatedness_ht = relatedness_ht.filter(relatedness_ht.relationship != UNRELATED)
    relatedness_ht = relatedness_ht.select("relationship", s=relatedness_ht.i.s).union(
        relatedness_ht.select("relationship", s=relatedness_ht.j.s)
    )
    relatedness_ht = relatedness_ht.group_by(relatedness_ht.s).aggregate(
        relationships=hl.agg.collect_as_set(relatedness_ht.relationship)
    )
    return relatedness_ht


def get_relationship_filter_expr(
    hard_filtered_expr: hl.expr.BooleanExpression,
    relationship: str,
    relationship_set: hl.expr.SetExpression,
) -> hl.expr.builders.CaseBuilder:
    """
    Returns case statement to populate relatedness filters in sample_filters struct

    :param hl.expr.BooleanExpression hard_filtered_expr: Boolean for whether sample was hard filtered. 
    :param str relationship: Relationship to check for. One of DUPLICATE_OR_TWINS, PARENT_CHILD, or SIBLINGS.
    :param hl.expr.SetExpression relationship_set: Set containing all possible relationship strings for sample.
    :return: Case statement used to population sample_filters related filter field.
    :rtype: hl.expr.builders.CaseBuilder
    """
    return (
        hl.case()
        .when(hard_filtered_expr, hl.null(hl.tbool))
        .when(hl.is_defined(relationship_set), relationship_set.contains(relationship))
        .default(False)
    )


# Release-related utils
def get_hists(mt: hl.MatrixTable, freeze: int) -> hl.MatrixTable:
    """
    Gets age (at recruitment; field 21022) and qual hists for UKBB.

    .. note::
        This function expects dense format data.

    :param MatrixTable mt: Input MatrixTable
    :param int freeze: One of the data freezes
    :return: MatrixTable with age and qual hists
    :rtype: MatrixTable
    """
    logger.warning(
        "This function expects dense data! Make sure you have run hl.experimental.densify."
    )
    logger.info("Getting age hists...")
    age_ht = get_age_ht(freeze)
    mt = mt.annotate_cols(age=age_ht[mt.col_key].age)
    mt = mt.annotate_rows(**age_hists_expr(mt.adj, mt.GT, mt.age))

    logger.info("Annotating with qual hists...")
    mt = mt.annotate_rows(qual_hists=qual_hist_expr(mt.GT, mt.GQ, mt.DP, mt.AD, mt.adj))

    logger.info(
        "Renaming hists (qual_hist_expr adds '_adj' suffix to adj hists) and returning..."
    )
    return mt.annotate_rows(
        qual_hists=hl.struct(
            **{
                i.replace("_adj", ""): mt.qual_hists[i]
                for i in mt.qual_hists
                if "_adj" in i
            }
        ),
        raw_qual_hists=hl.struct(
            **{i: mt.qual_hists[i] for i in mt.qual_hists if "_adj" not in i}
        ),
    )


def get_age_distributions(ht: hl.Table, bins: List[int] = [30, 80, 10]) -> str:
    """
    Get background distribution of sample ages (using field 21022, age at recruitment).

    :param Table ht: Table containing samples and sample ages.
    :param List[int]: Desired bin edges and bin size. Default: [30, 80, 10].
    :return: Pipe-delimited string with ages in bins (if using default bins: <30, 30-35, ..., 75-80, 80+).
    :rtype: str
    """
    age_hist_data = ht.aggregate(hl.agg.hist(ht.age, *bins))
    age_hist_data.bin_freq.insert(0, age_hist_data.n_smaller)
    age_hist_data.bin_freq.append(age_hist_data.n_larger)
    return age_hist_data.bin_freq


def index_globals(
    freq_meta: List[Dict[str, str]], label_groups: Dict
) -> Dict[str, int]:
    """
    Create a dictionary keyed by the specified label groupings with values describing the corresponding index of each grouping entry in the freq_meta array annotation.

    e.g., if freq_meta is [{'group': 'adj'}, {'group': 'raw'}, {'group': 'adj', 'pop': 'nfe'}], then this function will return 
    {'adj': 0, 'raw': 1, 'nfe_adj': 2}.

    :param List[Dict[str, str]] freq_meta: Ordered list containing dictionaries describing all the grouping combinations contained in the
        frequency row annotation.
    param Dict[str, List[str]] label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"]).
    :return: Dictionary keyed by specified label grouping combinations, with values describing the corresponding index
        of each grouping entry in the frequency array.
    :rtype: Dict[str, int]
    """
    combos = make_label_combos(label_groups)
    index_dict = {}

    for combo in combos:
        combo_fields = combo.split("_")
        for i, v in enumerate(freq_meta):
            if set(v.values()) == set(combo_fields):
                # NOTE: this is UKBB specific
                # This is to avoid updating the dictionary when the subpop is the same as the pop
                # E.g., this is to avoid updating the index for {'group': 'adj', 'pop': 'nfe'} with the
                # index for {'group': 'adj', 'pop': 'nfe', 'subpop': 'nfe'},
                if "subpop" in v and v["subpop"] == v["pop"]:
                    continue
                index_dict.update({f"{combo}": i})
    return index_dict


def make_freq_meta_index_dict(
    freq_meta: List[str],
    gnomad: bool,
    pops: List[str],
    subpops: List[str] = None,
    data_type: str = None,
    gnomad_nfe_subpops: List[str] = GNOMAD_NFE_SUBPOPS,
    gnomad_eas_subpops: List[str] = GNOMAD_EAS_SUBPOPS,
    groups: List[str] = GROUPS,
    gnomad_sexes: List[str] = SEXES_STR,
    ukbb_sexes: List[str] = SEXES,
) -> Dict[str, int]:
    """
    Makes a dictionary of the entries in the frequency array annotation, where keys are the grouping combinations and the values
    are the 0-based integer indices.

    :param List[str] freq_meta: Ordered list containing string entries describing all the grouping combinations contained in the
        frequency array annotation.
    :param bool gnomad: Whether to index a gnomAD sample freq_meta list.
    :param List[str] pops: List of global populations in frequency array. Used for both gnomAD and UKBB. 
        Can handle populations to unique to gnomAD/UKBB or a union of all population names.
    :param List[str] subpops: List of UKBB subpopulations in frequency array. Default is None.
    :param str data_type: Data type: one of "exomes" or "genomes". Required only for gnomAD data. Default is None.
    :param List[str] gnomad_nfe_subpops: List of nfe subpopulations in gnomAD. Default is GNOMAD_NFE_SUBPOPS.
    :param List[str] gnomad_eas_subpops: List of eas subpopulations in gnomAD. Default is GNOMAD_EAS_SUBPOPS.
    :param List[str] groups: Group names used to generate labels for high quality genotypes and all raw genotypes. Default is GROUPS.
    :param List[str] sexes: gnomAD v2 sample sexes ("male", "female") used in VCF export. Default is SEXES_STR. 
    :param List[str] ukbb_sexes: UKBB, gnomAD v3 sample sexes ("XX", "XY") used in VCF export. Default is SEXES.
    :return: Dictionary keyed by grouping combinations in the frequency array, with values describing the corresponding index
        of each grouping entry in the frequency array
    :rtype: Dict[str, int]
    """
    if gnomad:
        # v2 exomes used "male" and "female" labels
        if data_type == "exomes":
            sexes = gnomad_sexes
        # v3 genomes used "XX" and "XY" labels (same as UKBB)
        else:
            sexes = ukbb_sexes
    else:
        sexes = ukbb_sexes

    index_dict = index_globals(freq_meta, dict(group=groups))

    index_dict.update(index_globals(freq_meta, dict(group=groups, pop=pops)))
    index_dict.update(index_globals(freq_meta, dict(group=groups, sex=sexes)))
    index_dict.update(index_globals(freq_meta, dict(group=groups, pop=pops, sex=sexes)))

    if subpops:
        index_dict.update(
            index_globals(freq_meta, dict(group=groups, pop=pops, subpop=subpops))
        )

    if gnomad:
        index_dict.update(
            index_globals(
                freq_meta, dict(group=groups, pop=["nfe"], subpop=gnomad_nfe_subpops),
            )
        )
        index_dict.update(
            index_globals(
                freq_meta, dict(group=groups, pop=["eas"], subpop=gnomad_eas_subpops),
            )
        )

    return index_dict


def make_index_dict(
    t: Union[hl.MatrixTable, hl.Table],
    freq_meta_str: str,
    pops: List[str],
    subpops: List[str],
    data_type: str = None,
) -> Dict[str, int]:
    """
    Create a look-up Dictionary for entries contained in the frequency annotation array.

    :param Table ht: Table or MatrixTable containing freq_meta global annotation to be indexed
    :param str freq_meta: freq_meta global annotation to be indexed (freq_meta, gnomad_exomes_freq_meta, or gnomad_genomes_freq_meta)
    :param List[str] pops: List of global populations in frequency array. Used for both gnomAD and UKBB. 
        Can handle populations to unique to gnomAD/UKBB or a union of all population names.
    :param List[str] subpops: List of subpops in frequency array. Required for gnomAD populations.
    :param str data_type: Data type: one of "exomes" or "genomes". Required only for gnomAD data. Default is None.
    :return: Dictionary keyed by grouping combinations in the frequency array, with values describing the corresponding index
        of each grouping entry in the frequency array
    :rtype: Dict of str: int
    """
    freq_meta = hl.eval(t.globals[freq_meta_str])
    # check if indexing gnomAD data
    if "gnomad" in freq_meta_str:
        index_dict = make_freq_meta_index_dict(
            freq_meta, gnomad=True, pops=pops, subpops=subpops, data_type=data_type,
        )
    else:
        if subpops:
            index_dict = make_freq_meta_index_dict(
                freq_meta, gnomad=False, pops=pops, subpops=subpops
            )
        else:
            index_dict = make_freq_meta_index_dict(freq_meta, gnomad=False, pops=pops)
    return index_dict
