import argparse
import hail as hl
import logging
from typing import Union
from gnomad_hail.resources.grch38.intervals import lcr
from gnomad_hail.utils.generic import (
                                bi_allelic_expr, filter_low_conf_regions, filter_to_autosomes
                                )
from gnomad_hail.utils.sample_qc import compute_stratified_metrics_filter
from gnomad_hail.utils.slack import try_slack
from ukbb_qc.resources.basics import array_sample_map_ht_path, get_ukbb_data
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import ancestry_hybrid_ht_path, platform_pop_outlier_ht_path, qc_temp_data_prefix
from ukbb_qc.utils import annotate_interval_qc_filter, remove_hard_filter_samples


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("outlier_filter")
logger.setLevel(logging.INFO)


def filter_low_conf_regions(
    t: Union[hl.MatrixTable, hl.Table], filter_lcr: bool = True
    ) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filters low-confidence regions

    :param MatrixTable or Table mt: MatrixTable or Table to filter
    :param bool filter_lcr: Whether to filter LCR regions
    :param bool filter_decoy: Whether to filter decoy regions
    :param bool filter_segdup: Whether to filter Segdup regions
    :param bool filter_exome_low_coverage_regions: Whether to filter exome low confidence regions
    :param list of str high_conf_regions: Paths to set of high confidence regions to restrict to (union of regions)
    :return: MatrixTable or Table with low confidence regions removed
    :rtype: MatrixTable or Table
    """
    criteria = []
    if filter_lcr:
        lcr_intervals = lcr.ht()
        lcr = hl.import_locus_intervals(lcr_intervals_path, reference_genome="GRCh38")
        criteria.append(hl.is_missing(lcr[t.locus]))

    if criteria:
        filter_criteria = functools.reduce(operator.iand, criteria)
        t = t.filter_rows(filter_criteria) if isinstance(t, hl.MatrixTable) else t.filter(filter_criteria)

    return t


def run_sample_qc(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Filter input MatrixTable to bi-allelic sites, removes problematic intervals, and computes sample QC metrics

    :param MatrixTable mt: Input MatrixTable on which sample QC metrics need to be computed
    :return: MatrixTable filtered to autosomes and high-confidence regions, with computed sample QC column annotations
    :rtype: hl.MatrixTable
    """
    # Filtering to bi-allelic sites outside of low confidence regions
    mt = filter_to_autosomes(mt)
    mt = filter_low_conf_regions(mt)
    mt = mt.filter_rows(bi_allelic_expr(mt))

    logger.info("Starting sample QC...")
    mt = hl.sample_qc(mt)
    mt = mt.annotate_rows(variant_qc=hl.struct(af=hl.agg.mean(mt.GT.n_alt_alleles()) / 2))
    mt = mt.annotate_cols(sample_qc=mt.sample_qc.annotate(f_inbreeding=hl.agg.inbreeding(mt.GT, mt.variant_qc.af)))

    return mt.cols().select("sample_qc")


def main(args):

    hl.init(log="/outlier_detection.log", default_reference="GRCh38")

    data_source = args.data_source
    freeze = args.freeze
    pop_assignment_method = args.pop_assignment_method

    if args.run_mini_qc:
        # NOTE: we run outlier detection without adj filtration to get better separation between high and low quality samples
        # this is per Julia"s discussion with Konrad in #ukbb_qc
        # Need all workers for the mini qc
        mt = get_ukbb_data(data_source, freeze, split=True, adj=False)
        logger.info(f"Split hardcalls count: {mt.count()}")
        mt = annotate_interval_qc_filter(data_source, freeze, mt, autosomes_only=True)
        mt = mt.filter_rows(mt.interval_qc_pass)
        logger.info(f"mt count post interval filtering: {mt.count()}")

        logger.info("Filtering samples that fail hard filters...")
        mt = remove_hard_filter_samples(data_source, freeze, mt, non_refs_only=False)
        logger.info(f"mt count post filtering: {mt.count()}")

        logger.info("Running mini sample QC for platform- and population-specific filtering...")
        sample_qc_ht = run_sample_qc(mt)
        sample_qc_ht.write(
            qc_temp_data_prefix(data_source, freeze) + "outlier_sample_qc_intervals.ht", 
            args.overwrite
        )

    sample_qc_ht = hl.read_table(qc_temp_data_prefix(data_source, freeze) + "outlier_sample_qc_intervals.ht")
    strata = {}
    
    logger.info("Annotating platform assignments...")
    # NOTE using tranche as a proxy for platform in freeze 5
    #platform_ht = hl.read_table(platform_pca_results_ht_path(data_source, freeze))
    #sample_qc_ht = sample_qc_ht.annotate(qc_platform=platform_ht[sample_qc_ht.key].qc_platform)
    sample_map_ht = hl.read_table(array_sample_map_ht_path(data_source, freeze))
    
    logger.info("Annotating population assignments...")
    pop_ht = hl.read_table(ancestry_hybrid_ht_path(data_source, freeze))
    sample_qc_ht = sample_qc_ht.annotate(qc_pop=pop_ht[sample_qc_ht.key][pop_assignment_method])

    if not args.skip_platform_filter:
        strata["qc_pop"] = sample_qc_ht.qc_pop
    if not args.skip_population_filter:
        #strata["qc_platform"] = sample_qc_ht.qc_platform
        strata["qc_platform"] = sample_qc_ht.batch
    else:
        pop_assignment_method = None

    # Make qc_metrics a dict (needs to be dict for compute_stratified_metrics_filter)
    metrics = args.filtering_qc_metrics.split(",")
    qc_metrics = {}
    for m in metrics:
        qc_metrics[m] = sample_qc_ht.sample_qc[f"{m}"]

    # For each platform and population, aggregate sample QC metrics and calculate the MAD/mean/stdev
    logger.info("Flagging samples failing pop/platform-specific sample qc thresholds...")
    pop_platform_filter_ht = compute_stratified_metrics_filter(
        sample_qc_ht,
        qc_metrics,
        strata
    )
    pop_platform_filter_ht.write(
        platform_pop_outlier_ht_path(data_source, freeze, pop_assignment_method), 
        overwrite=args.overwrite
    )
    num_pass = pop_platform_filter_ht.aggregate(
        hl.agg.count_where(hl.len(pop_platform_filter_ht.qc_metrics_filters) == 0)
    )
    logger.info(f"{num_pass} exome samples found passing pop/platform-specific filtering")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--overwrite", help="Overwrite all data from this subset (default: False)", action="store_true")
    parser.add_argument("--slack_channel", help="Slack channel to post results and notifications to.")
    parser.add_argument("-s", "--data_source", help="Data source", choices=["regeneron", "broad"], default="broad")
    parser.add_argument("-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int)

    parser.add_argument("--run_mini_qc",
                                      help="Run mini sample qc needed for outlier filtering",
                                      action="store_true")
    parser.add_argument("--skip_platform_filter",
                                      help="Skip including platforms in outlier filtering",
                                      action="store_true")
    parser.add_argument("--skip_population_filter",
                                      help="Skip including population in outlier filtering",
                                      action="store_true")
    parser.add_argument("--pop_assignment_method", help="Population assignment method to use for stratification",
                        default="hybrid_pop", choices=["gnomad_pc_project_pop","HDBSCAN_pop_cluster","hybrid_pop"])
    parser.add_argument("--filtering_qc_metrics", help="List of QC metrics for filtering.", default=",".join(
    ["n_snp", "r_ti_tv", "r_insertion_deletion", "n_insertion", "n_deletion", "r_het_hom_var"]))

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
