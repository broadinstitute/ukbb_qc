from gnomad_hail import *
import hail as hl
from ukbb_qc.resources import *
from ukbb_qc.call_sex import *

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("hardcalls")
logger.setLevel(logging.INFO)


def main(args):
    hl.init()

    data_source = args.data_source
    freeze = args.freeze

    logger.info('IMPUTING SEX...')
    logger.info('Filtering to high-callrate, common, biallelic SNPs on relevant chromosomes...')
    mt = hl.read_matrix_table(raw_mt_path(data_source, freeze))
    mt = hl.variant_qc(mt)
    mt.describe()

    sex_check_intervals = [hl.parse_locus_interval(x, reference_genome='GRCh38') for x in ['chr20', 'chrX', 'chrY']]
    sex_mt = hl.filter_intervals(mt, sex_check_intervals)

    sex_mt = sex_mt.filter_rows((hl.len(sex_mt.alleles) == 2) & hl.is_snp(sex_mt.alleles[0], sex_mt.alleles[1]))
    sex_mt = sex_mt.filter_rows((hl.agg.fraction(hl.is_defined(sex_mt.GT)) > 0.99) & (sex_mt.variant_qc.AF[0] > 0.05))
    # NOTE: Deleted 'PASS' requirement for now -- let's see how this goes

    # NOTE: coverage estimation is now based on high-callrate, common, biallelic SNPs
    sex_ht = run_impute_sex(sex_mt, data_source, freeze)
    sex_colnames = ['f_stat', 'is_female', 'sex', 'normalized_Y_cov']
    sex_ht = sex_ht.select(*sex_colnames)
    mt = mt.annotate_cols(**sex_ht[mt.col_key])

    # TODO: check distributions here

    if args.write_hardcalls:
        logger.info("Generating hardcalls...")
        mt = get_ukbb_data(data_source, freeze, split=False, raw=True)
        ht = hl.read_table(hard_filters_ht_path(data_source, freeze))
        mt = annotate_adj(mt.select_cols(sex=ht[mt.row_key].sex))
        mt = mt.select_entries(GT=mt.GT, adj=mt.adj)  # Note: this is different from gnomAD hardcalls file because no PGT or PID
        mt = adjust_sex_ploidy(mt, mt.sex)
        mt = mt.select_cols().naive_coalesce(10000)
        mt.write(get_ukbb_data_path(data_source, freeze, hardcalls=True, split=False), args.overwrite)

    if args.split_hardcalls:
        logger.info("Running split_multi on the hardcalls...")
        mt = get_ukbb_data(data_source, freeze, split=False)
        mt = hl.split_multi_hts(mt)
        mt.write(get_ukbb_data_path(data_source, freeze, hardcalls=True, split=True), args.overwrite)

    if args.write_nonrefs:
        logger.info("Creating sparse MT with only non-ref genotypes...")
        mt = get_ukbb_data(data_source, freeze, split=False, raw=True).select_cols()
        mt = mt.annotate_entries(is_missing=hl.is_missing(mt.GT))
        mt = mt.filter_entries(mt.is_missing | mt.GT.is_non_ref())
        mt = annotate_adj(mt)
        mt = mt.naive_coalesce(10000)
        mt.write(get_ukbb_data_path(data_source, freeze, split=False, non_refs_only=True), args.overwrite)

    if args.split_nonrefs:
        logger.info("Running split_multi on non-ref MT...")
        mt = get_ukbb_data(data_source, freeze, split=False, non_refs_only=True)
        mt = mt.drop('PL')  # Note: I guess we need to check if PL has issues rather than just always dropping
        mt = hl.split_multi_hts(mt)
        mt = mt.filter_entries(mt.is_missing | mt.GT.is_non_ref())
        mt.write(get_ukbb_data_path(data_source, freeze, split=True, non_refs_only=True), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-o', '--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE)

    parser.add_argument('--write_hardcalls', help='Creates a hardcalls mt', action='store_true')
    parser.add_argument('--split_hardcalls', help='Creates a split hardcalls mt from the hardcalls mt', action='store_true')
    parser.add_argument('--write_nonrefs', help='Creates a sparse mt with only non-ref genotypes', action='store_true')
    parser.add_argument('--split_nonrefs', help='Creates a split sparse mt with only non-ref genotypes', action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)