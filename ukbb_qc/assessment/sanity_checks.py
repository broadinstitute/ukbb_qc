import logging
import hail as hl
from typing import Tuple


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
def sanity_check_release_mt(
    mt: hl.MatrixTable,
    subsets: List[str],
    missingness_threshold: float = 0.5,
    verbose: bool = False,
) -> None:
    """
    Perform a battery of sanity checks on a specified group of subsets in a MatrixTable containing variant annotations.

    Includes:
    - Summaries of % filter status for different partitions of variants
    - Histogram outlier bin checks
    - Checks on AC, AN, and AF annotations
    - Checks that subgroup annotation values add up to the supergroup annotation values
    - Checks on sex-chromosome annotations; and summaries of % missingness in variant annotations

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
        is_filtered=hl.len(ht.filters) > 0,
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
                for z in POPS_NAMES:
                    if z in i:
                        found.append(z)
        missing_pops = set(POP_NAMES) - set(found)
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
            sample_sum_check(ht, subset, dict(group=["adj"], pop=POP_NAMES), verbose)
            sample_sum_check(ht, subset, dict(group=["adj"], sex=SEXES), verbose)
            sample_sum_check(
                ht, subset, dict(group=["adj"], pop=POP_NAMES, sex=SEXES), verbose
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
