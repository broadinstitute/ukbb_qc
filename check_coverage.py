from gnomad_hail import *
from ukbb_qc.resources import *


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("check_coverage")
logger.setLevel(logging.INFO)


def compute_cov_by_sample(mt: hl.MatrixTable):
    '''
    Annotate sample-specific coverage metrics by grouping an MT by the exome calling region type
    (i.e., gnomAD interval, UKBB island, or UKBB adjacent region) and GC content bin (either absolute or relative)
    and then computing the following annotations: mean DP and % genotypes >= 20x coverage.

    NOTE: Grouping by region_type and gc_bin and then running entries() shuffles too much to be successful;
    the more convoluted alternative method used here avoids long runtimes and shuffle errors

    :param mt: MatrixTable containing region_type and gc_bin annotations for each locus
    :return: Table with coverage metrics computed for each combination of region types and GC bin content, for each sample
    :rtype: Table
    '''
    region_ht = (mt.annotate_cols(data=hl.agg.group_by(hl.struct(region_type=mt.region_type, gc_bin=mt[gc_bin_name]),
                                                       hl.struct(mean_dp=hl.agg.mean(mt.DP),
                                                                 pct_gt_20x=hl.agg.filter(hl.is_defined(mt.DP),
                                                                                          hl.agg.fraction(
                                                                                              mt.DP >= 20)))))).cols()
    region_ht = region_ht.persist()
    region_ht = region_ht.annotate(data_items=hl.array(region_ht.data)).drop('data')
    region_ht = region_ht.explode('data_items')
    region_ht = region_ht.annotate(region_type=region_ht.data_items[0]['region_type'],
                                   gc_bin=region_ht.data_items[0]['gc_bin'],
                                   mean_dp=region_ht.data_items[1]['mean_dp'],
                                   pct_gt_20x=region_ht.data_items[1]['pct_gt_20x']).drop('data_items')
    return region_ht


def main(args):
    hl.init(log='/check_regeneron_vcf.log')

    mt = hl.read_matrix_table(get_ukbb_data('regeneron', split=False, raw=True))
    if args.test_run:
        mt = hl.filter_intervals(mt, [hl.parse_locus_interval('chr20', reference_genome='GRCh38')])

    if not args.skip_compute_interval_annotations:
        logger.info('Adding calling interval annotations')
        ht = hl.import_table(ukbb_calling_intervals_summary, impute=True)
        ht = ht.key_by(locus=hl.locus_interval(ht.seqnames, ht.start, ht.end, includes_end=True, reference_genome='GRCh38'))
        ht = ht.select('target_type', 'region_type', 'target_id')

        logger.info('Adding GC bin annotations')
        gc_ht = hl.import_bed(gc_content_bed, reference_genome='GRCh38', skip_invalid_intervals=True)
        gc_ht = gc_ht.annotate(gc_bin=hl.float32(gc_ht.target)).drop('target')
        # Merge GC bins by interval
        gc_ht = gc_ht.annotate(start_locus=hl.locus(gc_ht.interval.start.contig, gc_ht.interval.start.position, reference_genome='GRCh38'))
        gc_ht = gc_ht.annotate(**ht[gc_ht.start_locus])
        gc_ht = gc_ht.filter(hl.is_defined(gc_ht.target_id))
        gc_ht = gc_ht.group_by(gc_ht.target_id).aggregate(mean_gc=hl.agg.mean(gc_ht.gc_bin))
        gc_ht = gc_ht.key_by('mean_gc').add_index()
        bin_size = hl.floor(gc_ht.count() / 10)  # Hard-coding 10 relative GC bins
        gc_ht = gc_ht.annotate(rel_gc_bin=hl.floor(hl.float(gc_ht.idx + 1) / bin_size))
        gc_ht = gc_ht.key_by('target_id')

        ht = ht.annotate(**gc_ht[ht.target_id]).drop('idx')
        ht = ht.annotate(abs_gc_bin=(hl.case()
                                     .when((ht.mean_gc < 10), "[0, 0.1)")
                                     .when((ht.mean_gc < 20), "[0.1, 0.2)")
                                     .when((ht.mean_gc < 30), "[0.2, 0.3)")
                                     .when((ht.mean_gc < 40), "[0.3, 0.4)")
                                     .when((ht.mean_gc < 50), "[0.4, 0.5)")
                                     .when((ht.mean_gc < 60), "[0.5, 0.6)")
                                     .when((ht.mean_gc < 70), "[0.6, 0.7)")
                                     .when((ht.mean_gc < 80), "[0.7, 0.8)")
                                     .when((ht.mean_gc < 90), "[0.8, 0.9)")
                                     .when((ht.mean_gc <= 100), "[0.9, 1]")
                                     .default(hl.null(hl.tstr))))
        ht.write(intervals_ht(args.data_source, args.freeze))

    ht = hl.read_table(intervals_ht(args.data_source, args.freeze))
    mt = mt.annotate_rows(**ht[mt.locus])
    ht = ht.annotate(chr=ht.locus.start.contig).key_by('target_id').distinct()

    # NOTE: no median aggregator yet in Hail, so no median_dp=hl.agg.median(mt.DP)
    if not args.skip_interval_coverage:
        logger.info('Computing coverage metrics per target')
        target_mt = (mt.group_rows_by('target_id')
                     .aggregate_entries(mean_dp=hl.agg.mean(mt.DP),
                                        pct_gt_20x=hl.agg.filter(hl.is_defined(mt.DP),
                                                                 hl.agg.fraction(mt.DP >= 20))).result())
        target_mt = target_mt.annotate_rows(target_mean_dp=hl.agg.mean(target_mt.mean_dp),
                                            target_pct_gt_20x=hl.agg.mean(target_mt.pct_gt_20x))
        target_ht = target_mt.rows()
        target_ht = target_ht.annotate(**ht[target_ht.target_id])
        target_ht.write(coverage_by_target(args.data_source, args.freeze), overwrite=args.overwrite)
        target_ht = hl.read_table(coverage_by_target(args.data_source, args.freeze))
        target_ht.export(coverage_by_target(args.data_source, args.freeze, ht=False))

    if not args.skip_sample_coverage:
        logger.info('Computing sample coverage metrics per region type')
        rel_ht = compute_cov_by_sample(mt, 'rel_gc_bin')
        rel_ht.write(coverage_by_region_type(args.data_source, args.freeze, 'rel_gc_bin'))
        rel_ht = hl.read_table(coverage_by_region_type(args.data_source, args.freeze, 'rel_gc_bin'))
        rel_ht.export(coverage_by_region_type(args.data_source, args.freeze, 'rel_gc_bin', ht=False))

        abs_ht = compute_cov_by_sample(mt, 'abs_gc_bin')
        abs_ht.write(coverage_by_region_type(args.data_source, args.freeze, 'abs_gc_bin'))
        abs_ht = hl.read_table(coverage_by_region_type(args.data_source, args.freeze, 'abs_gc_bin'))
        abs_ht.export(coverage_by_region_type(args.data_source, args.freeze, 'abs_gc_bin', ht=False))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite pre-existing data', action='store_true')
    parser.add_argument('--data_source', help='Regeneron or Broad', choices=['regeneron', 'broad'])
    parser.add_argument('--freeze', help='UKBB tranche', default=CURRENT_FREEZE)
    parser.add_argument('--skip_compute_interval_annotations', help='Skip computing interval annotations', action='store_true')
    parser.add_argument('--skip_interval_coverage', help='Skip calculating interval coverage statistics', action='store_true')
    parser.add_argument('--skip_sample_coverage', help='Skip calculating sample coverage statistics', action='store_true')
    parser.add_argument('--test_run', help='Filter to chr20 for testing purposes', action='store_true')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)


# NOTES from Daniel:
# We could assess these intervals to see which actually have coverage (easy) and if they actually have some relevant
# functional annotation (harder, because there are lots of annotations that might have been used)
