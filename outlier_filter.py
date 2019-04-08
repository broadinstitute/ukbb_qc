from gnomad_hail import *
import hail as hl
from gnomad_hail.utils.sample_qc import *
from ukbb_qc.resources import *

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("relatedness")
logger.setLevel(logging.INFO)


def run_sample_qc(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Filter MTs to bi-allelic sites and remove problematic intervals, and performs sample QC
    :param MatrixTable mt: MT on which sample QC metrics need to be computed
    :return: MT filtered to autosomes and high-confidence regions, with computed sample QC column annotations
    :rtype: MatrixTable
    """
    mt = filter_to_autosomes(mt) # filter_low_conf_regions(mt))
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)  # NOTE: this does not work on a split MT!
    mt = hl.sample_qc(mt)
    mt = mt.annotate_rows(variant_qc=hl.struct(af=hl.agg.mean(mt.GT.n_alt_alleles()) / 2))
    mt = mt.annotate_cols(sample_qc=mt.sample_qc.annotate(f_inbreeding=hl.agg.inbreeding(mt.GT, mt.variant_qc.af)))
    return mt


def main(args):
    data_source = args.data_source
    freeze = args.freeze

    mt = hl.read_matrix_table(hard_filters_mt_path(data_source, freeze))
    logger.info('Filtering samples that fail hard filters...')
    mt = mt.filter_cols(hl.len(mt.hard_filters) == 0).select_cols()

    logger.info('Running mini sample QC for platform- and population-specific filtering...')
    run_sample_qc(mt).cols().select('sample_qc').write(qc_temp_data_prefix(data_source, freeze) + 'outlier_sample_qc.ht', args.overwrite)
    sample_qc_ht = hl.read_table(qc_temp_data_prefix(data_source, freeze) + 'outlier_sample_qc.ht')

    strata = []
    if not args.skip_platform_filter:
        logger.info('Annotating platform assignments...')
        platform_ht = hl.read_table(platform_pca_results_ht_path(data_source, freeze))
        sample_qc_ht = sample_qc_ht.annotate(qc_platform=platform_ht[sample_qc_ht.key].qc_platform)
        strata.append('qc_platform')

    if not args.skip_population_filter:
        logger.info('Annotating population assignments...')
        pop_ht = hl.read_table(ancestry_pc_project_scores_ht_path(data_source, freeze))
        sample_qc_ht = sample_qc_ht.annotate(qc_pop=pop_ht[sample_qc_ht.key].pop.pop)
        strata.append('qc_pop')

    # For each platform and population, aggregate sample QC metrics and calculate the MAD/mean/stdev
    logger.info('Flagging samples failing pop/platform-specific sample qc thresholds...')
    pop_platform_filter_ht = compute_stratified_metrics_filter(
        sample_qc_ht,
        args.filtering_qc_metrics.split(","),
        strata
    )
    pop_platform_filter_ht.write(platform_pop_outlier_ht_path(data_source, freeze), overwrite=args.overwrite)

    num_pass = pop_platform_filter_ht.aggregate(hl.agg.count_where(hl.len(pop_platform_filter_ht.pop_platform_filters) == 0))
    logger.info(f'{num_pass} exome samples found passing pop/platform-specific filtering')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE)

    parser.add_argument('--skip_platform_filter',
                                      help="Skip including platforms in outlier filtering",
                                      action='store_true')
    parser.add_argument('--skip_population_filter',
                                      help="Skip including population in outlier filtering",
                                      action='store_true')
    parser.add_argument('--filtering_qc_metrics', help="List of QC metrics for filtering.", default=",".join(
    ['n_snp', 'r_ti_tv', 'r_insertion_deletion', 'n_insertion', 'n_deletion', 'r_het_hom_var']))

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)