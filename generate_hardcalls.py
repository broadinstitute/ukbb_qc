from gnomad_hail import *
import hail as hl
from ukbb_qc.resources import *
from ukbb_qc.call_sex import *

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("hardcalls")
logger.setLevel(logging.INFO)


def main(args):
    hl.init(log='/generate_hardcalls.log')

    data_source = args.data_source
    freeze = args.freeze

    if args.impute_sex:
        logger.info('Imputing sex...')
        mt = get_ukbb_data(data_source, freeze, split=False, raw=True)

        # NOTE: correct densified sparse mt for hardcalls
        mt = mt.key_rows_by('locus', 'alleles')
        mt = mt.annotate_entries(GT=hl.experimental.lgt_to_gt(mt.LGT, mt.LA))
        sex_check_intervals = [hl.parse_locus_interval(x, reference_genome='GRCh38') for x in ['chr20', 'chrX', 'chrY']]
        mt = hl.filter_intervals(mt, sex_check_intervals)
        mt = hl.variant_qc(mt)
        run_impute_sex(mt, data_source, freeze)
    # NOTE: check distributions here before continuing with hardcalls

    if args.write_hardcalls:
        logger.info("Generating hardcalls...")
        mt = get_ukbb_data(data_source, freeze, split=False, raw=True)
        ht = hl.read_table(sex_ht_path(data_source, freeze))
        mt = annotate_adj(mt.select_cols(sex=ht[mt.col_key].sex))
        mt = mt.select_entries(GT=mt.GT, adj=mt.adj)  # Note: this is different from gnomAD hardcalls file because no PGT or PID
        mt = adjust_sex_ploidy(mt, mt.sex)
        mt = mt.select_cols().naive_coalesce(10000)
        mt.write(get_ukbb_data_path(data_source, freeze, hardcalls=True, split=False), args.overwrite)

    if args.split_hardcalls:
        logger.info("Running split_multi on the hardcalls...")
        mt = get_ukbb_data(data_source, freeze, split=False, raw=False)
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
        if data_source == 'regeneron':
            mt = mt.drop('PL')  # Note: I guess we need to check if PL has issues rather than just always dropping
        mt = hl.split_multi_hts(mt)
        mt = mt.filter_entries(mt.is_missing | mt.GT.is_non_ref())
        mt.write(get_ukbb_data_path(data_source, freeze, split=True, non_refs_only=True), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-o', '--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)

    parser.add_argument('--impute_sex', help='Impute sex on raw MT (prerequisite for creating hardcalls)', action='store_true')
    parser.add_argument('--write_hardcalls', help='Creates a hardcalls mt', action='store_true')
    parser.add_argument('--split_hardcalls', help='Creates a split hardcalls mt from the hardcalls mt', action='store_true')
    parser.add_argument('--write_nonrefs', help='Creates a sparse mt with only non-ref genotypes', action='store_true')
    parser.add_argument('--split_nonrefs', help='Creates a split sparse mt with only non-ref genotypes', action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
