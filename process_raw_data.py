from gnomad_hail import *
from gnomad_hail.utils.sample_qc import *
from resources import *
import hail as hl
import argparse


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("process_raw")
logger.setLevel(logging.INFO)


# Note: The following functions are slightly modified from Laurent's myoseq sample qc should probably move to a common location
def compute_qc_mt(mt: hl.MatrixTable, min_af: float = 0.001, min_callrate: float = 0.99) -> hl.MatrixTable:
    """
    Returns MatrixTable for sample QC purposes
    Default criteria: callrate > 0.99, AF > 0.001, SNPs only, bi-allelics only
    :param MatrixTable mt: Raw MatrixTable to be filtered
    :param float min_af: Minimum allele frequency for variant filtering
    :param float min_callrate: Minimum callrate for variant filtering
    :return: MatrixTable for sample QC purposes
    :rtype: MatrixTable
    """
    qc_mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]) &
                        (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > min_af) &
                        (hl.agg.fraction(hl.is_defined(mt.GT)) > min_callrate))
    qc_mt = qc_mt.annotate_globals(
        qc_mt_params=hl.struct(
            min_af=min_af,
            min_callrate=min_callrate
        )
    )
    qc_mt.annotate_cols(sample_callrate=hl.agg.fraction(hl.is_defined(qc_mt.GT)))
    return qc_mt


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
        qc_mt = compute_qc_mt(get_ukbb_data(data_source, freeze, adj=True, split=False))
        qc_mt = qc_mt.naive_coalesce(5000)
        qc_mt.write(qc_mt_path(data_source, freeze), overwrite=args.overwrite)

        qc_mt = hl.read_matrix_table(qc_mt_path(data_source, freeze))
        qc_mt = filter_to_autosomes(qc_mt)
        qc_ht = hl.sample_qc(qc_mt).cols().select('sample_qc')
        qc_ht.write(qc_ht_path(data_source, freeze), overwrite=args.overwrite)

    if not args.skip_ld_prune_qc_mt:
        logger.info("LD prune qc MatrixTable sample QC...")
        qc_mt = hl.read_matrix_table(qc_mt_path(data_source, freeze))
        pruned_qc_mt = ld_prune_qc_mt(qc_mt)
        variants, samples = pruned_qc_mt.count()
        logger.info(f'{variants} variants found in LD-pruned MT')
        pruned_qc_mt = filter_to_adj(pruned_qc_mt)
        pruned_qc_mt.write(qc_mt_path(data_source, freeze, ld_pruned=True), overwrite=args.overwrite)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE)

    parser.add_argument('--skip_compute_qc_mt', help='Skip Compute matrix to be used in sample qc', action='store_true')
    parser.add_argument('--skip_ld_prune_qc_mt', help='Skip LD prunes the qc matrix', action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)