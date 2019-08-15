from gnomad_hail import *
import hail as hl
from gnomad_hail.utils.sample_qc import *
from ukbb_qc.resources import *


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("relatedness")
logger.setLevel(logging.INFO)


def filter_low_conf_regions(mt: Union[hl.MatrixTable, hl.Table], filter_lcr: bool = True) -> Union[hl.MatrixTable, hl.Table]:
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
    from ukbb_qc.resources import lcr_intervals_path

    criteria = []
    if filter_lcr:
        lcr = hl.import_locus_intervals(lcr_intervals_path, reference_genome="GRCh38")
        criteria.append(hl.is_missing(lcr[mt.locus]))

    if criteria:
        filter_criteria = functools.reduce(operator.iand, criteria)
        if isinstance(mt, hl.MatrixTable):
            mt = mt.filter_rows(filter_criteria)
        else:
            mt = mt.filter(filter_criteria)

    return mt


def run_sample_qc(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Filter MTs to bi-allelic sites and remove problematic intervals, and performs sample QC
    :param MatrixTable mt: MT on which sample QC metrics need to be computed
    :return: MT filtered to autosomes and high-confidence regions, with computed sample QC column annotations
    :rtype: MatrixTable
    """
    mt = filter_to_autosomes(mt)
    mt = filter_low_conf_regions(mt)
    mt = mt.filter_rows((~mt.was_split) & (hl.len(mt.alleles) == 2))
    logger.info('starting sample QC...')
    mt = hl.sample_qc(mt)
    mt = mt.annotate_rows(variant_qc=hl.struct(af=hl.agg.mean(mt.GT.n_alt_alleles()) / 2))
    mt = mt.annotate_cols(sample_qc=mt.sample_qc.annotate(f_inbreeding=hl.agg.inbreeding(mt.GT, mt.variant_qc.af)))

    return mt.cols().select('sample_qc')


def main(args):
    data_source = args.data_source
    freeze = args.freeze
    pop_assignment_method = args.pop_assignment_method

    if args.run_mini_qc:
        # NOTE: we run outlier detection without adj filtration to get better separation between high and low quality samples
        # this is per Julia's discussion with Konrad in #ukbb_qc
        # Need all workers for the mini qc
        mt = get_ukbb_data(data_source, freeze, split=False, adj=False)
        logger.info('Filtering samples that fail hard filters...')
        qc_ht = hl.read_table(hard_filters_ht_path(data_source, freeze)).key_by('s')
        mt = mt.filter_cols(hl.len(qc_ht[mt.col_key].hard_filters) == 0).select_cols()
        logger.info('Running mini sample QC for platform- and population-specific filtering...')
        sample_qc_ht = run_sample_qc(mt)
        sample_qc_ht.write(qc_temp_data_prefix(data_source, freeze) + 'outlier_sample_qc.ht', args.overwrite)

    sample_qc_ht = hl.read_table(qc_temp_data_prefix(data_source, freeze) + 'outlier_sample_qc.ht')
    strata = []
    if not args.skip_platform_filter:
        logger.info('Annotating platform assignments...')
        platform_ht = hl.read_table(platform_pca_results_ht_path(data_source, freeze))
        sample_qc_ht = sample_qc_ht.annotate(qc_platform=platform_ht[sample_qc_ht.key].qc_platform)
        strata.append('qc_platform')

    if not args.skip_population_filter:
        logger.info('Annotating population assignments...')
        pop_ht = hl.read_table(ancestry_hybrid_ht_path(data_source, freeze))
        sample_qc_ht = sample_qc_ht.annotate(qc_pop=pop_ht[sample_qc_ht.key][pop_assignment_method])
        strata.append('qc_pop')
    else:
        pop_assignment_method = None

    # For each platform and population, aggregate sample QC metrics and calculate the MAD/mean/stdev
    logger.info('Flagging samples failing pop/platform-specific sample qc thresholds...')
    pop_platform_filter_ht = compute_stratified_metrics_filter(
        sample_qc_ht,
        args.filtering_qc_metrics.split(","),
        strata
    )
    pop_platform_filter_ht.write(platform_pop_outlier_ht_path(data_source, freeze, pop_assignment_method), overwrite=args.overwrite)

    num_pass = pop_platform_filter_ht.aggregate(hl.agg.count_where(hl.len(pop_platform_filter_ht.pop_platform_filters) == 0))
    logger.info(f'{num_pass} exome samples found passing pop/platform-specific filtering')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)

    parser.add_argument('--run_mini_qc',
                                      help="Run mini sample qc needed for outlier filtering",
                                      action='store_true')
    parser.add_argument('--skip_platform_filter',
                                      help="Skip including platforms in outlier filtering",
                                      action='store_true')
    parser.add_argument('--skip_population_filter',
                                      help="Skip including population in outlier filtering",
                                      action='store_true')
    parser.add_argument('--pop_assignment_method', help="Population assignment method to use for stratification",
                        default='hybrid_pop', choices=['gnomad_pc_project_pop','HDBSCAN_pop_cluster','hybrid_pop'])
    parser.add_argument('--filtering_qc_metrics', help="List of QC metrics for filtering.", default=",".join(
    ['n_snp', 'r_ti_tv', 'r_insertion_deletion', 'n_insertion', 'n_deletion', 'r_het_hom_var']))

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
