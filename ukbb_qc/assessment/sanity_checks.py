import logging
from typing import Dict, List, Optional, Tuple

import hail as hl

from gnomad.assessment.sanity_checks import (
    generic_field_check,
    make_filters_sanity_check_expr,
    sample_sum_check,
)
from gnomad.resources.grch38.gnomad import SEXES
from gnomad.utils.filtering import filter_to_adj
from gnomad.utils.vcf import HISTS
from gnomad.utils.vcf import SEXES as SEXES_STR

from ukbb_qc.utils.utils import GNOMAD_EAS_SUBPOPS, GNOMAD_NFE_SUBPOPS

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("sanity_checks")
logger.setLevel(logging.INFO)


# Resources to check raw MT upon loading
def summarize_mt(mt: hl.MatrixTable) -> hl.Struct:
    """
    Gets a summary of variants in a MatrixTable.
    Prints number of variants to stdout, and checks that each chromosome has variant calls

    :param MatrixTable mt: Raw MatrixTable to be checked
    :return: Struct of MatrixTable variant summary
    :rtype: Struct
    """

    var_summary = hl.summarize_variants(mt, show=False)
    logger.info(f"Dataset has {var_summary.n_variants} variants")

    # check that all contigs have variant calls
    for contig in var_summary.contigs:
        if var_summary.contigs[contig] == 0:
            logger.warning(f"{contig} has no variants called")

    return var_summary


def check_adj(
    mt: hl.MatrixTable, mt_adj: hl.MatrixTable, gt_expr: hl.expr.CallExpression
) -> bool:
    """
    Checks if MatrixTable has been filtered using adj criteria by checking allele counts pre and post adj filtration

    :param MatrixTable mt: MatrixTable to be checked
    :param MatrixTable mt_adj: MatrixTable filtered using adj criteria
    :param hl.expr.CallExpression gt_expr: Field containing genotype information 
    :return: Bool of whether MatrixTable has been adj filtered
    :rtype: bool
    """

    pre = mt.aggregate_entries(hl.agg.counter(mt[gt_expr].n_alt_alleles()))
    logger.info(f"\nAllele distribution pre adj filtration: {pre}")
    post = mt_adj.aggregate_entries(hl.agg.counter(mt_adj[gt_expr].n_alt_alleles()))
    logger.info(f"\nAllele distribution post adj filtration: {post}")

    adj = False
    if sum(pre.values()) != sum(post.values()):
        adj = True

    return adj


def sample_check(
    ht: hl.Table, exp_ht: hl.Table, show_mismatch: bool = True,
) -> Tuple[bool, bool]:
    """
    Checks for sample mismatch between samples in two Tables.
    If there is a sample mismatch, writes unique samples to output txt file
    Assumes the keys of the two tables match (uses anti_join).

    :param Table ht: Table containing samples to be checked
    :param Table exp_ht: Table with one column containing expected samples
    :param bool show_mismatch: Boolean whether to print sample mismatches to stdout. Default is True
    :return: Tuple of bools [whether there were missing samples, whether there were extra samples]
    :rtype: Tuple[bool, bool]
    """
    # bool to store whether there are samples missing from ht
    missing = False
    # bool to store whether ht contains samples not in exp_ht
    extra = False

    missing_samples = exp_ht.anti_join(ht).select()
    n_missing_samples = missing_samples.count()
    extra_samples = ht.anti_join(exp_ht).select()
    n_extra_samples = extra_samples.count()

    if n_missing_samples > 0:
        missing = True
        logger.info(
            f"Total number of IDs that are not in the sample HT: {n_missing_samples}..."
        )
        if show_mismatch:
            missing_samples.show(n_missing_samples)

    if n_extra_samples > 0:
        extra = True
        logger.info(f"Total number of extra IDs in the sample HT: {n_extra_samples}...")
        if show_mismatch:
            extra_samples.show(n_extra_samples)

    return (missing, extra)


# Resources to check MT upon VCF export
def filters_sanity_check(ht: hl.Table) -> None:
    """
    Summarizes variants filtered under various conditions in input Table.

    Summarizes counts for:
        - Total number of variants
        - Fraction of variants removed due to:
            - Any filter
            - Inbreeding coefficient filter in combination with any other filter
            - AC0 filter in combination with any other filter
            - Random forest filtering in combination with any other filter
            - Monoallelic filter in combination with any other filter
            - Only inbreeding coefficient filter
            - Only AC0 filter
            - Only random forest filtering
            - Only monoallelic filter

    :param hl.Table ht: Input Table.
    :return: None
    :rtype: None
    """
    ht_explode = ht.explode(ht.filters)
    logger.info(
        f"hl.agg.counter filters: {ht_explode.aggregate(hl.agg.counter(ht_explode.filters))}"
    )
    # NOTE: in_problematic_region check will need to be updated if we get hg38 decoy file
    ht = ht.annotate(
        is_filtered=hl.len(ht.filters) > 0,
        in_problematic_region=hl.any(
            lambda x: x, [ht.info.lcr, ht.info.fail_interval_qc, ht.info.segdup]
        ),
    )

    def _filter_agg_order(
        ht: hl.Table,
        group_exprs: Dict[str, hl.expr.Expression],
        n_rows: int = None,
        n_cols: int = None,
        extra_filter_checks: Optional[Dict[str, hl.expr.Expression]] = None,
    ) -> None:
        """
        Performs sanity checks to measure percentages of variants filtered under different conditions.

        :param hl.Table ht: Input Table.
        :param hl.expr.Expression group_exprs: Dictionary of expressions to group the table by.
        :param int n_rows: Number of rows to show.
        :param int n_cols: Number of columns to show.
        :return: None
        """
        # NOTE: make_filters_sanity_check_expr returns a dict with %ages of variants filtered
        ht.group_by(**group_exprs).aggregate(
            **make_filters_sanity_check_expr(ht, extra_filter_checks)
        ).order_by(hl.desc("n")).show(n_rows, n_cols)

    logger.info(
        "Checking distributions of filtered variants amongst variant filters..."
    )
    # Add extra check for monoallelic variants to make_filters_sanity_check_expr (currently UKBB-specific filter)
    monoallelic_dict = {
        "frac_monoallelic": hl.agg.fraction(ht.filters.contains("MonoAllelic")),
        "frac_monoallelic_only": hl.agg.fraction(
            ht.filters.contains("MonoAllelic") & (ht.filters.length() == 1)
        ),
    }
    _filter_agg_order(
        ht, {"is_filtered": ht.is_filtered}, extra_filter_checks=monoallelic_dict
    )

    logger.info("Checking distributions of variant type amongst variant filters...")
    _filter_agg_order(ht, {"allele_type": ht.info.allele_type})

    logger.info(
        "Checking distributions of variant type and region type amongst variant filters..."
    )
    _filter_agg_order(
        ht,
        {
            "allele_type": ht.info.allele_type,
            "in_problematic_region": ht.in_problematic_region,
        },
        50,
        140,
    )

    logger.info(
        "Checking distributions of variant type, region type, and number of alt alleles amongst variant filters..."
    )
    _filter_agg_order(
        ht,
        {
            "allele_type": ht.info.allele_type,
            "in_problematic_region": ht.in_problematic_region,
            "n_alt_alleles": ht.info.n_alt_alleles,
        },
        50,
        140,
    )


def histograms_sanity_check(
    ht: hl.Table, verbose: bool, hists: List[str] = HISTS
) -> None:
    """
    Checks the number of variants that have nonzero values in their n_smaller and n_larger bins of quality histograms (both raw and adj).

    :param hl.Table ht: Input Table.
    :param bool verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param List[str] hists: List of variant annotation histograms.
    :return: None
    :rtype: None
    """
    for hist in hists:
        for suffix in ["", "raw"]:
            if suffix == "raw":
                logger.info("Checking raw qual hists...")
                hist = f"{hist}_{suffix}"
            else:
                logger.info("Checking adj qual hists...")

            # Check subfield == 0
            generic_field_check(
                ht,
                cond_expr=(ht.info[f"{hist}_n_smaller"] != 0),
                check_description=f"{hist}_n_smaller == 0",
                display_fields=[f"info.{hist}_n_smaller"],
                verbose=verbose,
            )
            if hist not in [
                "dp_hist_alt",
                "dp_hist_all",
            ]:  # NOTE: DP hists can have nonzero values in n_larger bin
                generic_field_check(
                    ht,
                    cond_expr=(ht.info[f"{hist}_n_larger"] != 0),
                    check_description=f"{hist}_n_larger == 0",
                    display_fields=[f"info.{hist}_n_larger"],
                    verbose=verbose,
                )


def raw_and_adj_sanity_checks(ht: hl.Table, subsets: List[str], verbose: bool):
    """
    Performs sanity checks on raw and adj data in input Table.

    Checks that:
        - Raw AC, AN, AF are not 0
        - Adj AN is not 0 and AC and AF are not negative
        - Raw values for AC, AN, nhomalt in each sample subset are greater than or equal to their corresponding adj values

    :param hl.Table ht: Input Table.
    :param List[str] subsets: List of sample subsets.
    :param bool verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :return: None
    :rtype: None
    """
    for subfield in ["AC", "AF"]:
        # Check raw AC, AF > 0
        # NOTE: some sites should fail the raw AC > 0 for UKBB (but only in the low ten thousands range)
        # We generate the release MT based on the raw MT, but the frequency HT calculates frequency only on samples that pass QC
        generic_field_check(
            ht,
            cond_expr=(ht.info[f"{subfield}_raw"] <= 0),
            check_description=f"{subfield}_raw > 0",
            display_fields=[f"info.{subfield}_raw"],
            verbose=verbose,
        )
        # Check adj AC, AF >=0
        generic_field_check(
            ht,
            cond_expr=(ht.info[f"{subfield}_adj"] < 0),
            check_description=f"{subfield}_adj >= 0",
            display_fields=[f"info.{subfield}_adj", "filters"],
            verbose=verbose,
        )

    # Check raw AN > 0
    generic_field_check(
        ht,
        cond_expr=(ht.info.AN_raw <= 0),
        check_description="AN_raw > 0",
        display_fields=["info.AN_raw"],
        verbose=verbose,
    )

    # Check adj AN >= 0
    generic_field_check(
        ht,
        cond_expr=(ht.info.AN_adj < 0),
        check_description="AN_adj >= 0",
        display_fields=["info.AN_adj"],
        verbose=verbose,
    )

    for subset in subsets:
        if subset != "":
            subset += "_"
        for subfield in ["AC", "AN", "nhomalt"]:
            # Check AC_raw >= AC adj
            generic_field_check(
                ht,
                cond_expr=(
                    ht.info[f"{subset}{subfield}_raw"]
                    < ht.info[f"{subset}{subfield}_adj"]
                ),
                check_description=f"{subset}{subfield}_raw >= {subset}{subfield}_adj",
                display_fields=[
                    f"info.{subset}{subfield}_raw",
                    f"info.{subset}{subfield}_adj",
                ],
                verbose=verbose,
            )


def frequency_sanity_checks(ht: hl.Table, subsets: List[str], verbose: bool) -> None:
    """
    Performs sanity checks on frequency data in input Table.

    Checks:
        - Number of sites where gnomAD exome frequency is equal to the gnomAD genome frequency (both raw and adj)
        - Number of sites where the UKBB exome frequency is equal to the gnomAD exome frequency (both raw and adj)
        - Number of sites where the UKBB exome frequency is equal to the gnomAD genome frequency (both raw and adj)

    Also performs small spot checks:
        - Counts total number of sites where the gnomAD exome allele count annotation is defined (both raw and adj)
        - Counts total number of sites where the gnomAD genome allele count annotation is defined (both raw and adj)
        - Counts total number of sites where the UKBB exome allele count annotation is defined (both raw and adj)
        
    :param hl.Table ht: Input Table.
    :param List[str] subsets: List of sample subsets.
    :param bool verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :return: None
    :rtype: None
    """
    for subset in subsets:
        if subset == "":
            for subfield in ["AC", "AN", "nhomalt"]:

                for group_type in ["adj", "raw"]:
                    logger.info(f"{group_type} checks -- gnomAD/gnomAD")

                    generic_field_check(
                        ht,
                        cond_expr=(
                            ht.info[f"gnomad_exomes_{subfield}_{group_type}"]
                            == ht.info[f"gnomad_genomes_{subfield}_{group_type}"]
                        ),
                        check_description=f"gnomad_exomes_{subfield}_{group_type} != gnomad_genomes_{subfield}_{group_type}",
                        display_fields=[
                            f"info.gnomad_exomes_{subfield}_{group_type}",
                            f"info.gnomad_genomes_{subfield}_{group_type}",
                        ],
                        verbose=verbose,
                    )
        else:
            for subfield in ["AC", "AN", "nhomalt"]:
                logger.info("raw checks -- gnomAD/UKBB")
                generic_field_check(
                    ht,
                    cond_expr=(
                        ht.info[f"{subfield}_raw"]
                        == ht.info[f"{subset}_{subfield}_raw"]
                    ),
                    check_description=f"{subfield}_raw != {subset}_{subfield}_raw",
                    display_fields=[
                        f"info.{subfield}_raw",
                        f"info.{subset}_{subfield}_raw",
                    ],
                    verbose=verbose,
                    show_percent_sites=True,
                )
                logger.info("adj checks -- gnomAD/UKBB")
                generic_field_check(
                    ht,
                    cond_expr=(
                        ht.info[f"{subfield}_adj"]
                        == ht.info[f"{subset}_{subfield}_adj"]
                    ),
                    check_description=f"{subfield}_adj != {subset}_{subfield}_adj",
                    display_fields=[
                        f"info.{subfield}_adj",
                        f"info.{subset}_{subfield}_adj",
                    ],
                    verbose=verbose,
                    show_percent_sites=True,
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
            total_defined_ukb_AC=hl.agg.count_where(hl.is_defined(ht.info.AC_adj)),
            total_defined_ukb_AC_raw=hl.agg.count_where(hl.is_defined(ht.info.AC_raw)),
        )
    )
    logger.info(f"Frequency spot check counts: {freq_counts}")


def sample_sum_sanity_checks(
    ht: hl.Table,
    subsets: List[str],
    info_metrics: List[str],
    verbose: bool,
    ukbb_pops: Dict[str, str],
    gnomad_exomes_pops: Dict[str, str],
    gnomad_genomes_pops: Dict[str, str],
) -> None:
    """
    Performs sanity checks on sample sums in input Table.

    Computes afresh the sum of annotations for a specified group of annotations, and compare to the annotated version;
    displays results from checking the sum of the specified annotations in the terminal.

    Also checks that annotations for all expected sample populations are present (both for gnomAD and UKBB).

    :param hl.Table ht: Input Table.
    :param List[str] subsets: List of sample subsets.
    :param List[str] info_metrics: List of metrics in info struct of input Table.
    :param bool verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param Dict[str, str] ukbb_names: Dict with UKBB population names (keys) and population descriptions (values).
    :param Dict[str, str] gnomad_exomes_pops: Dict with gnomAD v2 exomes population names (keys) and population descriptions (values).
    :param Dict[str, str] gnomad_genomes_pops: Dict with gnomAD v3 genomes population names (keys) and population descriptions (values).
    :return: None
    :rtype: None
    """
    # Check if pops are present
    # Get list of all pops present in HT
    for subset in subsets:
        if "gnomad" in subset:
            if "exomes" in subset:
                sexes = SEXES_STR

                # Remove subpops here -- they have a different format in the info annotations
                # and are checked later in this function
                # genomes do not have subpops yet, so this only applies to exomes
                pops = gnomad_exomes_pops
                for subpop in GNOMAD_NFE_SUBPOPS + GNOMAD_EAS_SUBPOPS:
                    pops.pop(subpop, None)
            else:
                sexes = SEXES
                pops = gnomad_genomes_pops
            pop_adjusted = list(
                set([x for x in info_metrics if (subset in x) and ("raw" not in x)])
            )
        else:
            sexes = SEXES
            pops = ukbb_pops
            pop_adjusted = list(
                set(
                    [
                        x
                        for x in info_metrics
                        if (("adj" in x) and ("gnomad" not in x) and ("raw" not in x))
                    ]
                )
            )
        pop_adjusted = [i.replace("_adj", "") for i in pop_adjusted]

        # NOTE: Added filter to remove subset here
        # We don't retain any gnomAD subset info in the UKBB release files
        # This filter is a fix to prevent sample sum checks from
        # checking for HGDP/TGP population labels (gnomAD v3)
        # HGDP/TGP population labels are more granular than the global pops
        # and are not included in the UKBB release
        pop_found = ht[f"{subset + '_' if subset != '' else subset}freq_meta"].filter(
            lambda x: x.contains("pop") & ~x.contains("subset")
        )
        pop_found = list(hl.eval(pop_found.group_by(lambda x: x["pop"])).keys())
        for pop in pop_found:
            no_pop = True
            for i in pop_adjusted:
                if pop in i:
                    no_pop = False
            if no_pop:
                pop_found.remove(pop)
                logger.warning(
                    f"{pop} found in {subset} subset freq_meta but not in info_metrics!"
                )

        # Print any missing pops to terminal
        missing_pops = set(pops.keys()) - set(pop_found)
        if len(missing_pops) != 0:
            logger.warning(f"Missing {missing_pops} pops in {subset} subset!")

        # Perform sample sum checks
        sample_sum_check(
            ht, subset, dict(group=["adj"], pop=list(set(pop_found))), verbose
        )
        sample_sum_check(ht, subset, dict(group=["adj"], sex=sexes), verbose)
        sample_sum_check(
            ht,
            subset,
            dict(group=["adj"], pop=list(set(pop_found)), sex=sexes),
            verbose,
        )

        if "gnomad" in subset and "exomes" in subset:
            # Adjust subpops to those found in subset
            # This is checking v2 exomes subpops
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


def sex_chr_sanity_checks(
    ht: hl.Table, info_metrics: List[str], contigs: List[str], verbose: bool
) -> None:
    """
    Performs sanity checks for annotations on the sex chromosomes.

    Checks:
        - That metrics for chrY variants in female samples are NA and not 0
        - That nhomalt counts are equal to female nhomalt counts for all non-PAR chrX variants

    :param hl.Table ht: Input Table.
    :param List[str] info_metrics: List of metrics in info struct of input Table.
    :param List[str] contigs: List of contigs present in input Table.
    :param bool verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :return: None
    :rtype: None
    """
    female_metrics = [x for x in info_metrics if "_female" in x or "_XX" in x]

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

    ht_x = hl.filter_intervals(ht, [hl.parse_locus_interval("chrX")])
    ht_xnonpar = ht_x.filter(ht_x.locus.in_x_nonpar())
    n = ht_xnonpar.count()
    logger.info(f"Found {n} X nonpar sites")

    logger.info("Check (nhomalt == nhomalt_female) for X nonpar variants:")
    female_metrics = [x for x in female_metrics if "nhomalt" in x]
    for metric in female_metrics:
        standard_field = metric.replace("_female", "").replace("_XX", "")
        generic_field_check(
            ht_xnonpar,
            (ht_xnonpar.info[f"{metric}"] != ht_xnonpar.info[f"{standard_field}"]),
            f"{metric} == {standard_field}",
            [f"info.{metric}", f"info.{standard_field}"],
            verbose,
        )


def missingness_sanity_checks(
    ht: hl.Table,
    info_metrics: List[str],
    non_info_metrics: List[str],
    n_sites: int,
    missingness_threshold: float,
) -> None:
    """
    Checks amount of missingness in all row annotations.

    Prints metric to terminal if more than missingness_threshold% of annotations for that metric are missing.

    :param hl.Table ht: Input Table.
    :param List[str] info_metrics: List of metrics in info struct of input Table.
    :param List[str] non_info_metrics: List of row annotations minus info struct from input Table.
    :param int n_sites: Number of sites in input Table.
    :param float missingness_threshold: Upper cutoff for allowed amount of missingness.
    :return: None
    :rtype: None
    """
    logger.info(
        f"Missingness threshold (upper cutoff for what is allowed for missingness checks): {missingness_threshold}"
    )
    metrics_frac_missing = {}
    for x in info_metrics:
        metrics_frac_missing[x] = hl.agg.sum(hl.is_missing(ht.info[x])) / n_sites
    for x in non_info_metrics:
        metrics_frac_missing[x] = hl.agg.sum(hl.is_missing(ht[x])) / n_sites
    output = ht.aggregate(hl.struct(**metrics_frac_missing))

    n_fail = 0
    for metric, value in dict(output).items():
        message = f"missingness check for {metric}: {100 * value}% missing"
        if value > missingness_threshold:
            logger.info(f"FAILED {message}")
            n_fail += 1
        else:
            logger.info(f"Passed {message}")
    logger.info(f"{n_fail} missing metrics checks failed")


def sanity_check_release_mt(
    mt: hl.MatrixTable,
    subsets: List[str],
    ukbb_pops: Dict[str, str],
    gnomad_exomes_pops: Dict[str, str],
    gnomad_genomes_pops: Dict[str, str],
    missingness_threshold: float = 0.5,
    verbose: bool = False,
    n_alt_alleles_hist_start: int = 1,
    n_alt_alleles_hist_end: int = 1000,
    n_alt_alleles_hist_bins: int = 100,
) -> None:
    """
    Perform a battery of sanity checks on a specified group of subsets in a MatrixTable containing variant annotations.

    Includes:
    - Summaries of % filter status for different partitions of variants
    - Histogram outlier bin checks
    - Checks on AC, AN, and AF annotations
    - Checks that subgroup annotation values add up to the supergroup annotation values
    - Checks on sex-chromosome annotations; and summaries of % missingness in variant annotations

    :param MatrixTable mt: MatrixTable containing variant annotations to check.
    :param List[str] subsets: List of subsets to be checked.
    :param ukbb_pops: Dict with UKBB population names (keys) and population descriptions (values).
    :param gnomad_exomes_pops: Dict with gnomAD v2 exomes population names (keys) and population descriptions (values).
    :param gnomad_genomes_pops: Dict with gnomAD v3 genomes population names (keys) and population descriptions (values).
    :param float missingness_threshold: Upper cutoff for allowed amount of missingness. Default is 0.5
    :param bool verbose: If True, display top values of relevant annotations being checked, regardless of whether check
        conditions are violated; if False, display only top values of relevant annotations if check conditions are violated.
    :param int n_alt_alleles_hist_start: Start of range for `n_alt_alleles` histogram. Default is 1.
    :param int n_alt_alleles_hist_end: End of range for `n_alt_alleles` histogram. Default is 1000.
    :param int n_alt_alleles_hist_bins: Number of bins for `n_alt_alleles` histogram. Default is 100.
    :return: None (terminal display of results from the battery of sanity checks).
    :rtype: None
    """
    # Perform basic checks -- number of variants, number of contigs, number of samples,
    # number of alternate alleles histogram
    logger.info("BASIC SUMMARY OF INPUT TABLE:")
    n_samples = mt.count_cols()
    ht = mt.rows()
    n_sites = ht.count()
    contigs = ht.aggregate(hl.agg.collect_as_set(ht.locus.contig))
    n_alt_alleles_hist = ht.aggregate(
        hl.agg.hist(
            ht.info.n_alt_alleles,
            start=n_alt_alleles_hist_start,
            end=n_alt_alleles_hist_end,
            bins=n_alt_alleles_hist_bins,
        )
    )
    logger.info(f"Found {n_sites} sites in contigs {contigs} in {n_samples} samples")
    logger.info(f"n_alt_alleles histogram: {n_alt_alleles_hist}")

    logger.info("VARIANT FILTER SUMMARIES:")
    filters_sanity_check(ht)

    # NOTE: This check won't work in the 455k tranche
    # The `histograms_sanity_check` code checks `_n_smaller` for all hists and
    # `_n_larger` for every hist except the DP hists
    # ^ we dropped all of these annotations in the 455k tranche
    # logger.info("HISTOGRAM CHECKS:")
    # histograms_sanity_check(ht, verbose=verbose)

    logger.info("RAW AND ADJ CHECKS:")
    raw_and_adj_sanity_checks(ht, subsets, verbose)

    logger.info("FREQUENCY CHECKS:")
    frequency_sanity_checks(ht, subsets, verbose)

    # Pull row annotations from HT
    info_metrics = list(ht.row.info)
    non_info_metrics = list(ht.row)
    non_info_metrics.remove("info")

    logger.info("SAMPLE SUM CHECKS:")
    sample_sum_sanity_checks(
        ht,
        subsets,
        info_metrics,
        verbose,
        ukbb_pops,
        gnomad_exomes_pops,
        gnomad_genomes_pops,
    )

    logger.info("SEX CHROMOSOME ANNOTATION CHECKS:")
    sex_chr_sanity_checks(ht, info_metrics, contigs, verbose)

    logger.info("MISSINGNESS CHECKS:")
    missingness_sanity_checks(
        ht, info_metrics, non_info_metrics, n_sites, missingness_threshold
    )


def sanity_check_release_patch(
    mt: hl.MatrixTable, meta_ht: hl.Table, verbose: bool = False,
) -> None:
    """
    Perform sanity checks on release patch MatrixTable.

    Includes:
    - Total sample and variant count (should be 454699 samples, 72220 rows)
    - Homozygote count in VCF MT (created using raw MT) comparison to homozygote count produced when applying homalt hotfix
    - Frequency recalculation and comparison to VCF MT frequencies
    - Comparison of VCF MT frequencies and previously calculated, incorrectly adjusted 455K frequencies 
    (to confirm the frequencies have been updated)

    :param hl.MatrixTable mt: Release patch MatrixTable.
    :param hl.Table meta_ht: Table with sample metadata information.
    :param bool verbose: If True, display top values of relevant annotations being checked, regardless of whether check
        conditions are violated; if False, display only top values of relevant annotations if check conditions are violated.
    :return: None (terminal display of sanity checks results).
    :rtype: None
    """
    logger.info("CHECKING TOTAL COUNTS:")
    logger.info(mt.count())

    logger.info(
        "RE-APPLYING HOMALT HOTFIX (WITHOUT HETNONREF PATCH) AND COMPARING HOMOZYGOTE COUNTS..."
    )
    # Add homalt hotfix
    mt = mt.annotate_entries(
        GT_adj=hl.if_else(
            mt.GT.is_het() & (mt.info.AF_adj > 0.01) & (mt.AD[1] / mt.DP > 0.9),
            hl.call(1, 1),
            mt.GT,
        )
    )
    mt = mt.annotate_rows(homalt_stats=hl.agg.call_stats(mt.GT_adj, mt.alleles))

    logger.info("NOTE that this check should fail for a small number of rows (~80)...")
    # Check where homozygote count after applying hotfix is the same as the raw homozygote count
    # This check will fail for sites that had sample genotypes adjusted in ONLY the frequency code
    # or ONLY the previous release code
    # These sites needed either frequencies or genotypes adjusted and are unrelated to the het nonref issue
    ht = mt.rows()
    ht = ht.annotate(
        # Check the raw homozygote counts because these are the genotypes that will exist in the release VCF
        # (release VCF has no sample/variant filters and is not filtered to adj)
        het_nonref_hotfix_nhomalt=ht.info.nhomalt_raw,
        hotfix_nhomalt=ht.homalt_stats.homozygote_count[1],
    )
    generic_field_check(
        ht=ht,
        cond_expr=(ht.het_nonref_hotfix_nhomalt == ht.hotfix_nhomalt),
        check_description="homozygote_counts_homalt_hotfix != homozygote_counts_hetnonref_homalt_hotfix",
        display_fields=["het_nonref_hotfix_nhomalt", "hotfix_nhomalt"],
        verbose=verbose,
    )

    logger.info("RECALCULATING FREQUENCY AND COMPARING WITH FREQUENCIES IN MT...")
    mt.drop("GT_adj")
    mt = mt.annotate_cols(**meta_ht[mt.col_key])
    # Filter to unrelated, high quality samples WITH UKBB batch number defined ONLY
    mt = mt.filter_cols(
        (mt.sample_filters.high_quality)
        & hl.is_defined(mt.ukbb_meta.batch)
        & ~mt.sample_filters.related
    )
    # Filter to adj
    mt = filter_to_adj(mt)
    # Recalculate frequencies using `call_stats`
    mt = mt.annotate_rows(freq=hl.agg.call_stats(mt.GT, mt.alleles))
    ht = mt.rows()
    ht = ht.annotate(
        # Take 1st index of frequencies calculated using call_stats because 0th index is for ref allele
        recalc_AC=ht.freq.AC[1],
        recalc_AF=ht.freq.AF[1],
        recalc_AN=ht.freq.AN,
        recalc_nhomalt=ht.freq.homozygote_count[1],
        adj_AC=ht.info.AC_adj,
        adj_AF=ht.info.AF_adj,
        adj_AN=ht.info.AN_adj,
        adj_nhomalt=ht.info.nhomalt_adj,
    )
    for freq in ["AC", "AF", "AN", "nhomalt"]:
        generic_field_check(
            ht,
            cond_expr=(ht[f"recalc_{freq}"] != ht[f"adj_{freq}"]),
            check_description=f"recalculated_{freq} == adj_{freq}",
            display_fields=[f"recalc_{freq}", f"adj_{freq}",],
            verbose=verbose,
        )

    logger.info(
        "COMPARING VCF MT FREQUENCIES WITH PREVIOUSLY CALCULATED FREQUENCIES  (INCORRECT HOM ALT DEPLETION ADJUSTMENT ON HET NON REF VARAINTS)..."
    )
    # Overwrote HT at release HT path with frequencies recalculated from first patch release,
    # which is why this HT is in temp
    release_ht = (
        hl.read_table(
            "gs://broad-ukbb/broad.freeze_7/temp/broad.freeze_7.release.sites.ht"
        )
        .select_globals()
        .select("freq")
    )
    ht = ht.annotate(
        prev_AC=release_ht[ht.key].freq[1].AC,
        prev_AF=release_ht[ht.key].freq[1].AF,
        prev_AN=release_ht[ht.key].freq[1].AN,
        prev_nhomalt=release_ht[ht.key].freq[1].homozygote_count,
        raw_AC=ht.info.AC_raw,
        raw_AF=ht.info.AF_raw,
        raw_AN=ht.info.AN_raw,
        raw_nhomalt=ht.info.nhomalt_raw,
    )
    # NOTE: Expect most of these sites to pass
    # Not all variants included in the patch release required frequency adjustment
    # E.g., only ~18K of the ~61K variants should fail the recalc_AC != prev_AC check
    # Most of these 18K are het nonref impacted variants
    # A small subset (~10) of the 18K are variants common only in the 300K
    for freq in ["AC", "AF", "AN", "nhomalt"]:
        generic_field_check(
            ht,
            cond_expr=(ht[f"raw_{freq}"] != ht[f"prev_{freq}"]),
            check_description=f"recalculated_raw_{freq} == prev_raw_{freq}",
            display_fields=[f"raw_{freq}", f"prev_{freq}"],
            verbose=verbose,
        )


def vcf_field_check(
    mt: hl.MatrixTable,
    header_dict: Dict[str, Dict[str, Dict[str, str]]],
    row_annotations: List[str],
    entry_annotations: List[str],
    hists: List[str] = HISTS,
) -> bool:
    """
    Checks that all VCF fields and descriptions are present in input MatrixTable and VCF header dictionary.

    :param hl.MatrixTable mt: Input MatrixTable to be exported to VCF.
    :param Dict[str, Dict[str, Dict[str, str]]] header_dict: VCF header dictionary.
    :param List[str] row_annotations: List of row annotations in MatrixTable.
    :param List[str] entry_annotations: List of entry annotations in MatrixTable.
    :param List[str] hists: List of variant histogram annotations. Default is HISTS.
    :return: Bool with whether all expected fields and descriptions are present.
    :rtype: bool
    """
    # Confirm all VCF fields/descriptions are present before exporting
    hist_fields = []
    for hist in hists:
        hist_fields.extend(
            [
                f"{hist}_bin_freq",
                f"{hist}_n_smaller",
                f"{hist}_n_larger",
                f"{hist}_raw_bin_freq",
                f"{hist}_raw_n_smaller",
                f"{hist}_raw_n_larger",
            ],
        )

    missing_fields = []
    missing_descriptions = []
    for item in ["info", "format", "filter"]:
        if item == "info":
            annots = row_annotations
        elif item == "format":
            annots = entry_annotations
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
                        f"{field} in MT info field has empty description in VCF header!"
                    )
                    temp_missing_descriptions.append(field)
            except KeyError:
                logger.warning(
                    f"{field} in MT info field does not exist in VCF header!"
                )
                # NOTE: some hists are not exported, so ignoring here
                # END entry is also not exported (removed during densify)
                if (field not in hist_fields) and (field != "END"):
                    temp_missing_fields.append(field)

        missing_fields.extend(temp_missing_fields)
        missing_descriptions.extend(temp_missing_descriptions)

    if len(missing_fields) != 0 or len(missing_descriptions) != 0:
        logger.error(
            "Some fields are either missing or missing descriptions in the VCF header! Please reconcile."
        )
        logger.error(f"Missing fields: {missing_fields}")
        logger.error(f"Missing descriptions: {missing_descriptions}")
        return False

    logger.info("Passed VCF fields check!")
    return True
