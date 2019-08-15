from gnomad_hail import *
from ukbb_qc.resources import *
from ukbb_qc.utils import *
import hail as hl
import argparse

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("interval_qc")
logger.setLevel(logging.INFO)


def main(args):
    hl.init(log='/interval_qc.log')

    data_source = args.data_source
    freeze = args.freeze

    mt = get_ukbb_data(data_source, freeze, raw=True, split=False, adj=True)

    if args.chromosome is not None:
        mt = hl.filter_intervals(mt, [hl.parse_locus_interval(f'chr{args.chromosome}', reference_genome='GRCh38')])
    else:
        mt = filter_to_autosomes(mt)

    ht = hl.import_table(ukbb_calling_intervals_summary, impute=True)
    ht = ht.key_by(locus=hl.locus_interval(ht.seqnames, ht.start, ht.end, includes_end=True, reference_genome='GRCh38'))
    ht = ht.select('target_type', 'region_type', 'target_id')
    mt = mt.annotate_rows(**ht[mt.locus])
    ht = ht.annotate(chr=ht.locus.start.contig).key_by('target_id').distinct()

    target_mt = (mt.group_rows_by('target_id')
                 .partition_hint(1000)
                 .aggregate_entries(mean_dp=hl.agg.filter((hl.is_defined(mt.DP)) & (~hl.is_nan(mt.DP)),
                                                          hl.agg.mean(mt.DP)),
                                    pct_gt_20x=hl.agg.filter(hl.is_defined(mt.DP),
                                                             hl.agg.fraction(mt.DP >= 20))).result())
    target_mt = target_mt.naive_coalesce(1000)
    target_mt = target_mt.checkpoint(get_mt_checkpoint_path(data_source, freeze, "coverage_by_target_samples"),
                                     overwrite=args.overwrite)
    target_mt = target_mt.annotate_rows(target_num_mean_dp_defined=hl.agg.count_where(~hl.is_nan(target_mt.mean_dp)),
                                        target_mean_dp=hl.agg.filter(~hl.is_nan(target_mt.mean_dp),
                                                                     hl.agg.mean(target_mt.mean_dp)),
                                        target_pct_gt_20x=hl.agg.mean(target_mt.pct_gt_20x),
                                        pct_samples_15x=hl.agg.fraction(target_mt.mean_dp >= 15),
                                        pct_samples_20x=hl.agg.fraction(target_mt.mean_dp >= 20),
                                        pct_samples_25x=hl.agg.fraction(target_mt.mean_dp >= 25),
                                        pct_samples_30x=hl.agg.fraction(target_mt.mean_dp >= 30))

    target_ht = target_mt.rows()
    target_ht = target_ht.annotate(**ht[target_ht.target_id])
    target_ht = target_ht.naive_coalesce(100)

    target_ht.write(interval_qc_path(data_source, freeze, chrom=args.chromosome), overwrite=args.overwrite)
    target_ht = hl.read_table(interval_qc_path(data_source, freeze))
    target_ht.describe()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)

    parser.add_argument("-c", "--chromosome", type=int, choices=range(1,23),
                        help="If set it will only run the given chromosome")

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
