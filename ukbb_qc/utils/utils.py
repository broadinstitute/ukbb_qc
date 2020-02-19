import hail as hl
import logging
from typing import Union
from gnomad_hail.utils.generic import filter_to_autosomes, get_reference_genome
from gnomad_hail.resources.grch38.intervals import lcr
from ukbb_qc.resources.basics import get_ukbb_data, raw_mt_path
from ukbb_qc.resources.sample_qc import f_stat_sites_path, interval_qc_path, qc_temp_data_prefix, qc_sites_path


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("utils")
logger.setLevel(logging.INFO)


# Sample resources
def remove_hard_filter_samples(
    data_source: str, freeze: int, t: Union[hl.MatrixTable, hl.Table]
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Removes samples that failed hard filters from MatrixTable or Table. 
    Assumes MatrixTable's col key/Table's row key is sample.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param MatrixTable/able t: Input MatrixTable or Table
    :return: MatrixTable or Table with samples removed
    :rtype: MatrixTable or Table
    """

    # read in hard filters ht
    ht = hl.read_table(hard_filters_ht_path(data_source, freeze))

    # remove unnecessary fields from hard filters ht
    ht = ht.select("hard_filters")

    # get number of hard filtered samples
    hf_samples = ht.aggregate(hl.agg.count_where(hl.len(ht.hard_filters) > 0))
    logger.info(f"Removing {hf_samples} samples that failed hard filters")

    # remove hard filtered samples
    ht = ht.filter(hl.len(ht.hard_filters) == 0)
    if isinstance(t, hl.MatrixTable):
        t = t.filter_cols(hl.is_defined(ht[t.col_key]))
        t = t.annotate_rows(non_refs=hl.agg.count_where(t.GT.is_non_ref()))
        t = t.filter_rows(t.non_refs > 0).drop("non_refs")
    else:
        t.filter((hl.len(ht[t.key].hard_filters) == 0))
    return t


# Variant resources
def annotate_interval_qc_filter(
    data_source: str,
    freeze: int,
    t: Union[hl.MatrixTable, hl.Table],
    cov_filter_field: str = "pct_samples_20x",
    XY_cov_filter_field: str = "pct_samples_10x",
    pct_samples: float = 0.85,
    autosomes_only: bool = False
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Removes poorly covered intervals from a MatrixTable or Table.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param MatrixTable/able t: Input MatrixTable or Table
    :param float cov_filter_field: Coverage field to use for filtering of autosomes/sex chr par/female X
    :param float XY_cov_filter_field: Coverage field to use for filtering of male X and Y
    :param float pct_samples: Percent of samples with coverage greater than `cov_filter_field` over the interval for filtering
    :return: MatrixTable or Table with samples removed
    :rtype: MatrixTable or Table
    """
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
                & (interval_qc_sex_ht[cov_filter_field]["XX"] > pct_samples)
                & (interval_qc_sex_ht[XY_cov_filter_field]["XY"] > pct_samples)
            )
            | (
                (interval_qc_sex_ht.interval.start.contig == "chrY")
                & (interval_qc_sex_ht[XY_cov_filter_field]["XY"] > pct_samples)
            )
        )
        good_intervals_ht = good_intervals_ht.select().union(good_sex_intervals_ht.select())

    if isinstance(t, hl.MatrixTable):
        t = t.annotate_rows(
            interval_qc_pass=hl.is_defined(good_intervals_ht[t.locus])
        )
    else:
        t = t.annotate(
            interval_qc_pass=hl.is_defined(good_intervals_ht[t.locus])
        )

    return t


# Sites resources
def get_sites(
    af_cutoff: float, callrate_cutoff: float, greater: bool = True,  
    pass_interval_qc: bool = True, autosomes_only: bool = True
    ) -> hl.Table:
    """
    Returns a Table of sites that meet af_cutoff from the 200K (tranche 2/freeze 5) callset. 
    
    :param float af_cutoff: Desired AF cutoff
    :param float callrate_cutoff: Desired callrate cutoff
    :param bool greater: Whether to return variants with AFs higher than the cutoff.
    :param bool pass_interval_qc: Whether to only return sites that pass interval QC. 
    :param bool autosomes_only: Whether to return autosomal sites only. Also sets interval QC flag.
    :return: None
    :rtype: None
    """
    data_source = "broad"
    freeze = 5
    mt = hl.read_matrix_table(get_ukbb_data(data_source, freeze, split=True))

    if autosomes_only:
        mt = filter_to_autosomes(mt)

    if pass_interval_qc:
        mt = annotate_interval_qc_filter('broad', 5, mt, autosomes_only=autosomes_only)
        mt = mt.filter_rows(mt.interval_qc_pass)
    
    # Add callstats (for AF) and filter on AF cutoff
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.call_rate > callrate_cutoff)
    if greater:
        mt = mt.filter_rows(mt.variant_qc.AF[1] > af_cutoff)
    else:
        mt = mt.filter_rows(mt.variant_qc.AF[1] < af_cutoff)

    ht = mt.rows()
    ht = ht.transmute(AF=ht.variant_qc.AF[1])
    return ht


def calculate_fstat_sites() -> None:
    """
    Writes a Table with high callrate, common, biallelic SNPs in regions that pass interval QC on chromosome X.
    This Table is designed to be used as an interval filter in sex imputation.
    NOTE: This function generates sites only for tranche 2/freeze 5, which is the last dataset that contains AF.

    :return: None
    :rtype: None
    """
    data_source = "broad"
    freeze = 5
    mt = hl.read_matrix_table(raw_mt_path(data_source, freeze, densified=True))
    mt = mt.key_rows_by('locus', 'alleles')
    mt = hl.filter_intervals(mt,
                            [hl.parse_locus_interval(x_contig, reference_genome=get_reference_genome(mt.locus).name)
                            for x_contig in get_reference_genome(mt.locus).x_contigs])
    mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]))
    mt = annotate_interval_qc_filter('broad', 5, mt, autosomes_only=False)
    mt = mt.filter_rows(mt.interval_qc_pass)
    mt = mt.transmute_entries(GT=hl.experimental.lgt_to_gt(mt.LGT, mt.LA))
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows((hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99) & (mt.variant_qc.AF[1] > 0.001))
    ht = mt.rows()
    ht = ht.transmute(AF=ht.variant_qc.AF[1])
    ht = ht.write(f_stat_sites_path())


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
    lcr_intervals = lcr.ht()
    ht = ht.filter(
            hl.is_missing(lcr[ht.key])
    )

    logger.info("Adding info annotations to QC sites HT")
    # NOTE: info ht is hard coded from tranche 2/freeze 5
    info_ht = hl.read_matrix_table("gs://broad-ukbb/broad.freeze_5/temp/broad.freeze_5.sites.ht")
    ht = ht.annotate(info=info_ht[ht.key].info)
    ht.write(qc_sites_path())

