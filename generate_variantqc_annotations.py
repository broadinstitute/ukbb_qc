from gnomad_hail import *
from gnomad_qc.annotations.generate_qc_annotations import *
from ukbb_qc.resources import *
import argparse

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variantqc_annotations")
logger.setLevel(logging.INFO)


# NOTE: removed {"release_samples_raw": mt.meta.release} from sample_group_filters
def generate_call_stats(mt: hl.MatrixTable) -> hl.Table:
    """
    Add qc_callstats for 3 categories (high-quality samples, release, and all in VDS)
    """
    sample_group_filters = {
        "qc_samples_raw": mt.meta.high_quality,
        "all_samples_raw": True
    }
    mt = mt.select_cols(**sample_group_filters)
    mt = unphase_mt(mt.select_rows())
    call_stats_expression = []
    for group in sample_group_filters.keys():
        callstats = hl.agg.filter(mt[group], hl.agg.call_stats(mt.GT, mt.alleles))
        call_stats_expression.append(callstats.annotate(meta={'group': group}))

    return mt.annotate_rows(qc_callstats=call_stats_expression).rows()

# NOTE: altered from gnomAD version to create HTs for annotation
def annotate_truth_data(mt: hl.MatrixTable) -> hl.Table:
    """
    Writes bi-allelic sites MT with the following annotations:
     - truth_data (presence in Omni, HapMap, 1KG high conf SNVs, Mills)

    :param MatrixTable mt: Full MT
    :return: Table with qc annotations
    :rtype: Table
    """
    mt = mt.select_rows()

    truth_mtes = {
        'hapmap': hapmap_mt_path(),
        'omni': omni_mt_path(),
        'mills': mills_mt_path(),
        'kgp_high_conf_snvs': kgp_high_conf_snvs_mt_path()
    }
    truth_htes = {key: hl.split_multi_hts(hl.read_matrix_table(path).repartition(1000).rows(), left_aligned=False)
                  for key, path in truth_mtes.items()}
    truth_htes = truth_htes.update({'ukbb_array': var_annotations_ht_path(data_source, freeze, 'array_concordance')})  # TODO: formalize code to create this resource

    return mt.annotate_rows(truth_data=hl.struct(**{root: hl.is_defined(truth_ht[mt.row_key])
                                                    for root, truth_ht in truth_htes.items()})).rows()


def generate_array_concordant_variants(mt: hl.MatrixTable) -> hl.Table:
    """

    :param mt:
    :return: Table containing concordant variants found in both the UKBB array data and the UKBB exome data
    :rtype: Table
    """

    mt = hl.read_matrix_table(array_mt_path(liftover=True))
    return

# def generate_de_novos(mt: hl.MatrixTable, fam_file: str, freq_data: hl.Table) -> hl.Table:
#     mt = mt.select_cols()
#     fam_ht = read_fam(fam_file).key_by()
#     fam_ht = fam_ht.select(s=[fam_ht.s, fam_ht.pat_id, fam_ht.mat_id]).explode('s').key_by('s')
#     mt = mt.filter_cols(hl.is_defined(fam_ht[mt.s]))
#     mt = mt.select_rows()
#     mt = hl.split_multi_hts(mt)
#     mt = mt.annotate_rows(family_stats=freq_data[mt.row_key].family_stats)
#     ped = hl.Pedigree.read(fam_file, delimiter='\\t')
#
#     de_novo_table = hl.de_novo(mt, ped, mt.family_stats[0].unrelated_qc_callstats.AF[1])
#     de_novo_table = de_novo_table.key_by('locus', 'alleles').collect_by_key('de_novo_data')
#
#     return de_novo_table



def main(args):
    hl.init(log='/generate_variantqc_annotations.log')

    data_source = args.data_source
    freeze = args.freeze

    # if args.vep:  # CPU-hours: 250 (E), 600 (G)
    #     mt = get_gnomad_data(data_type).rows().select()
    #     hl.vep(mt, vep_config).write(annotations_ht_path(data_type, 'vep'), args.overwrite)
    #
    #     mt = get_gnomad_data(data_type).rows().select()
    #     hl.vep(mt, vep_config, csq=True).write(annotations_ht_path(data_type, 'vep_csq'), args.overwrite)

    if args.generate_allele_data:
        mt = get_ukbb_data(data_source, freeze, split=False)
        generate_allele_data(mt).write(var_annotations_ht_path(data_source, freeze, 'allele_data'), overwrite=args.overwrite)

    if args.generate_qc_annotations:
        # Turn on spark speculation: --properties 'spark:spark.speculation=true,spark:spark.speculation.quantile=0.9,spark:spark.speculation.multiplier=3'
        # TODO: replace this and following 2 lines with get_ukbb_data(data_source, freeze, non_refs_only=True, meta_root='meta') once 'PL' is retained for Broad split nonrefs
        mt = get_ukbb_data(data_source, freeze, split=False, non_refs_only=True, meta_root='meta')
        mt = hl.split_multi_hts(mt)
        mt = mt.filter_entries(mt.is_missing | mt.GT.is_non_ref())
        mt = generate_qc_annotations(mt, all_annotations=args.calculate_all_annotations, medians=args.calculate_medians)
        mt.write(var_annotations_ht_path(data_source, freeze, 'qc_stats'), stage_locally=True, overwrite=args.overwrite)

    # if args.generate_qual_hists:  # CPU-hours: 4000 (E), 8000 (G)
    #     mt = get_gnomad_data(data_type, raw=True, split=False, release_samples=True)
    #     ht = generate_qual_hists(mt)
    #     write_temp_gcs(ht, annotations_ht_path(data_type, 'qual_hists'), args.overwrite)

    if args.generate_call_stats:
        mt = get_ukbb_data(data_source, freeze, meta_root='meta')
        generate_call_stats(mt).write(var_annotations_ht_path(data_source, freeze, 'call_stats'), overwrite=args.overwrite)

    if args.generate_family_stats:  # CPU-hours: 8K (E), 13K (G)
        mt = get_ukbb_data(data_source, freeze, meta_root='meta')
        ht, sample_table = generate_family_stats(mt, inferred_ped_path(data_source, freeze), args.include_adj_family_stats)
        ht.write(var_annotations_ht_path(data_source, freeze, 'family_stats'), stage_locally=True, overwrite=args.overwrite)
        sample_table.write(sample_annotations_table_path(data_source, freeze, 'family_stats'), stage_locally=True, overwrite=args.overwrite)

    # if args.generate_de_novos:  # (2.2 min/part @ 100K = 3K CPU-hours) + (7.4 m/p = 12K) + (34 m/p = ~44K) = 59K
    #     # Turn on spark speculation?
    #     mt = get_gnomad_data(data_type, raw=True, split=False)
    #     freq_data = hl.read_table(annotations_ht_path(data_type, 'family_stats'))
    #     mt = generate_de_novos(mt, fam_path(data_type), freq_data)
    #     mt.write(annotations_ht_path(data_type, 'de_novos'), args.overwrite)

    if args.annotate_truth_data:
        mt = get_ukbb_data(data_source, freeze, meta_root=None)
        annotate_truth_data(mt).write(var_annotations_ht_path(data_source, freeze, 'truth_data'), overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE)
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')

    # parser.add_argument('--vep', help='Runs VEP', action='store_true')
    parser.add_argument('--generate_allele_data', help='Calculates allele data', action='store_true')
    parser.add_argument('--generate_qc_annotations', help='Calculates QC annotations', action='store_true')
    # parser.add_argument('--generate_qual_hists', help='Calculates GQ, DP, AB histograms per variant', action='store_true')
    parser.add_argument('--generate_call_stats', help='Calculates call stats', action='store_true')
    parser.add_argument('--generate_family_stats', help='Calculates family stats', action='store_true')
    parser.add_argument('--include_adj_family_stats', help='Also calculate family stats for adj genotypes', action='store_true')
    # parser.add_argument('--generate_de_novos', help='Calculates de novo data', action='store_true')
    parser.add_argument('--annotate_truth_data', help='Annotates MT with truth data', action='store_true')

    parser.add_argument('--calculate_medians', help='Calculate metric medians (warning: slow)', action='store_true')
    parser.add_argument('--calculate_all_annotations', help='Calculation many more annotations (warning: slow)', action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)

