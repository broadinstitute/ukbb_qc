from gnomad_hail import *
from resources import *


def main(args):
    hl.init()

    data_source = args.data_source
    freeze = args.freeze

    if args.write_hardcalls:
        mt = hl.read_matrix_table(hard_filters_mt_path(data_source, freeze))
        mt = mt.select_entries(GT=mt.GT, adj=mt.adj)  # Note: this is different from gnomAD hardcalls file because no PGT or PID
        mt = adjust_sex_ploidy(mt, mt.sex)
        mt = mt.select_cols().naive_coalesce(10000)
        mt.write(get_ukbb_data_path(data_source, freeze, hardcalls=True, split=False), args.overwrite)

    if args.split_hardcalls:
        mt = get_ukbb_data(data_source, freeze, split=False)
        mt = hl.split_multi_hts(mt)
        mt.write(get_ukbb_data_path(data_source, freeze, hardcalls=True, split=True), args.overwrite)

    if args.write_nonrefs:  # CPU-hours: 600 (E)
        mt = get_ukbb_data(data_source, freeze, split=False, raw=True).select_cols()
        mt = mt.annotate_entries(is_missing=hl.is_missing(mt.GT))
        mt = mt.filter_entries(mt.is_missing | mt.GT.is_non_ref())
        mt = annotate_adj(mt)
        mt = mt.naive_coalesce(10000)
        mt.write(get_ukbb_data_path(data_source, freeze, split=False, non_refs_only=True), args.overwrite)

    if args.split_nonrefs:  # CPU-hours: 300 (E)
        mt = get_ukbb_data(data_source, freeze, split=False, non_refs_only=True)
        mt = mt.drop('PL')  # Note: I guess we need to check if PL has issues rather than just always dropping
        mt = hl.split_multi_hts(mt)
        mt = mt.filter_entries(mt.is_missing | mt.GT.is_non_ref())
        mt.write(get_ukbb_data_path(data_source, freeze, split=True, non_refs_only=True), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--data_source', help='Source of the data, either broad or regeneron')
    parser.add_argument('--freeze', help='Data freeze to use', default=CURRENT_FREEZE)

    parser.add_argument('--write_hardcalls', help='Creates a hardcalls mt', action='store_true')
    parser.add_argument('--split_hardcalls', help='Creates a split hardcalls mt from the hardcalls mt', action='store_true')
    parser.add_argument('--write_nonrefs', help='Creates a sparse mt with only non-ref genotypes', action='store_true')
    parser.add_argument('--split_nonrefs', help='Creates a split sparse mt with only non-ref genotypes', action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)