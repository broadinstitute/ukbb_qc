from gnomad_hail import *
from gnomad_hail.utils.sample_qc import *
from ukbb_qc.resources import *
import hail as hl
import argparse
from ukbb_qc.sanity_checks import *


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("process_raw")
logger.setLevel(logging.INFO)


def ld_prune_qc_mt(qc_mt: hl.MatrixTable, ld_r2: float = 0.1) -> hl.MatrixTable:
    """
    Returns MatrixTable for sample QC purposes
    :param MatrixTable mt: Raw MatrixTable to be filtered
    :param float ld_r2: r2 to LD prune the QC matrix to
    :return: LD pruned MatrixTable for sample QC purposes
    :rtype: MatrixTable
    """
    pruned_ht = hl.ld_prune(qc_mt.GT, r2=ld_r2)
    qc_mt = qc_mt.filter_rows(hl.is_defined(pruned_ht[qc_mt.row_key]))
    qc_mt = qc_mt.annotate_globals(qc_mt_params=qc_mt.qc_mt_params.annotate(ld_r2=ld_r2))
    qc_mt.annotate_cols(sample_callrate=hl.agg.fraction(hl.is_defined(qc_mt.GT)))
    return qc_mt


def main(args):
    hl.init(log='/process_raw.log')

    data_source = args.data_source
    freeze = args.freeze

    if not args.skip_compute_qc_mt:
        logger.info("Filtering to bi-allelic, high-callrate, common SNPs for sample QC...")
        mt = get_ukbb_data(data_source, freeze, raw=True, adj=True, split=False)

        lcr = get_lcr_intervals()
        mt = mt.filter_rows(hl.is_missing(lcr[mt.row_key]))
        # Todo: changed to use Laurent's code for the qc_mt, change is that it uses inbreeding_coeff and hardy_weinberg_threshold, is this OK? Do we want to apply_hard_filters?
        # Todo: add segdup and decoy filter
        qc_mt = get_qc_mt(mt, min_af=0.001, min_callrate=0.99, apply_hard_filters=False, ld_r2=None, filter_lcr=False,
            filter_decoy=False, filter_segdup=False)
        qc_mt = qc_mt.naive_coalesce(5000)
        qc_mt.write(qc_mt_path(data_source, freeze), overwrite=args.overwrite)

        qc_mt = hl.read_matrix_table(qc_mt_path(data_source, freeze))
        qc_mt = filter_to_autosomes(qc_mt)
        qc_ht = hl.sample_qc(qc_mt).cols().select('sample_qc')
        qc_ht = qc_ht.transmute(sample_qc=qc_ht.sample_qc.select('call_rate', 'gq_stats', 'dp_stats'))
        qc_ht.write(qc_ht_path(data_source, freeze), overwrite=args.overwrite)

    if not args.skip_ld_prune_qc_mt:
        logger.info("LD prune qc MatrixTable sample QC...")
        qc_mt = hl.read_matrix_table(qc_mt_path(data_source, freeze))
        qc_mt = qc_mt.unfilter_entries()
        pruned_qc_mt = ld_prune_qc_mt(qc_mt)
        variants, samples = pruned_qc_mt.count()
        logger.info(f'{variants} variants found in LD-pruned MT')
        pruned_qc_mt.write(qc_mt_path(data_source, freeze, ld_pruned=True), overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)

    parser.add_argument('--skip_compute_qc_mt', help='Skip Compute matrix to be used in sample qc', action='store_true')
    parser.add_argument('--skip_ld_prune_qc_mt', help='Skip LD prunes the qc matrix', action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
