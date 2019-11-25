from gnomad_hail import *
import hail as hl
from ukbb_qc.resources import *
from ukbb_qc.call_sex import *
from gnomad_qc.annotations.generate_qc_annotations import *
from ukbb_qc.utils import *

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("hardcalls")
logger.setLevel(logging.INFO)


def annotate_adj(
        mt: hl.MatrixTable,
        adj_gq: int = 20,
        adj_dp: int = 10,
        adj_ab: float = 0.2,
        haploid_adj_dp: int = 10
) -> hl.MatrixTable:
    """
    Annotate genotypes with adj criteria (assumes diploid)
    Defaults correspond to gnomAD values.
    """
    return mt.annotate_entries(adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, adj_gq, adj_dp, haploid_adj_dp))


def main(args):
    hl.init(log='/generate_hardcalls.log')

    data_source = args.data_source
    freeze = args.freeze

    if args.impute_sex:
        logger.info('Imputing sex...')
        mt = get_ukbb_data(data_source, freeze, split=False, raw=True)

        # NOTE: correct densified sparse mt for hardcalls
        mt = mt.key_rows_by('locus', 'alleles')
        # NOTE: as of 11/25 densify no longer converts LGT to GT
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
        mt = mt.annotate_entries(GT=hl.experimental.lgt_to_gt(mt.LGT, mt.LA))
        mt = annotate_adj(mt.select_cols(sex=ht[mt.col_key].sex))
        mt = mt.select_entries(GT=mt.GT, adj=mt.adj)  # Note: this is different from gnomAD hardcalls file because no PGT or PID
        mt = adjust_sex_ploidy(mt, mt.sex)
        mt = mt.select_cols().naive_coalesce(15000)
        mt.write(get_ukbb_data_path(data_source, freeze, hardcalls=True, split=False), args.overwrite)

    if args.split_hardcalls:
        logger.info("Running split_multi on the hardcalls...")
        #mt = get_ukbb_data(data_source, freeze, split=False, raw=False)
        #mt = hl.split_multi_hts(mt)
        mt = get_ukbb_data(data_source, freeze, split=False, raw=True)

        # Add allele data to mt
        allele_data = hl.struct(nonsplit_alleles=mt.alleles,
                            has_star=hl.any(lambda a: a == '*', mt.alleles))
        mt = mt.annotate_rows(allele_data=allele_data.annotate(**add_variant_type(mt.alleles)))

        sex_ht = hl.read_table(sex_ht_path(data_source, freeze))
        mt = hl.experimental.sparse_split_multi(mt)
        mt = annotate_adj(mt.select_cols(sex=sex_ht[mt.col_key].sex))
        mt = mt.select_entries(GT=mt.GT, adj=mt.adj)
        mt = adjust_sex_ploidy(mt, mt.sex)
        mt = mt.select_cols().naive_coalesce(15000)
        mt = mt.checkpoint(get_ukbb_data_path(data_source, freeze, hardcalls=True, split=True), args.overwrite)

        # Finish generating allele data
        ht = mt.rows().select('allele_data')
        allele_type = (hl.case()
                   .when(hl.is_snp(ht.alleles[0], ht.alleles[1]), 'snv')
                   .when(hl.is_insertion(ht.alleles[0], ht.alleles[1]), 'ins')
                   .when(hl.is_deletion(ht.alleles[0], ht.alleles[1]), 'del')
                   .default('complex')
                   )
        ht = ht.annotate(allele_data=ht.allele_data.annotate(allele_type=allele_type,
                                                              was_mixed=ht.allele_data.variant_type == 'mixed'))
        ht.write(var_annotations_ht_path(data_source, freeze, 'allele_data'), args.overwrite)

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
