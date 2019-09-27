from gnomad_hail import *
from gnomad_qc.annotations.generate_qc_annotations import *
from ukbb_qc.resources import *
from ukbb_qc.utils import *
import argparse

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variantqc_annotations")
logger.setLevel(logging.INFO)

# TODO: Remove hardcoded paths

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
        # 'hapmap': hapmap_mt_path(),
        'omni': omni_mt_path(),
        'mills': mills_mt_path(),
        'kgp_high_conf_snvs': kgp_high_conf_snvs_mt_path(),
        'ukbb_array': 'gs://broad-ukbb/broad.freeze_4/temp/exome_subset_concordance.rekeyed.mt'
    }

    truth_htes = {key: hl.split_multi_hts(hl.read_matrix_table(path).repartition(1000).rows(), left_aligned=False)
                  for key, path in truth_mtes.items()}
    # TODO: formalize code to create this resource
    truth_htes.update({'hapmap': hl.read_table(hapmap_ht_path()),
                       'sib_singletons': hl.read_table('gs://broad-ukbb/broad.freeze_4/variant_qc/variant_annotations/sibling_singletons.train.ht'),
                       'ukbb_array_con':hl.read_table('gs://broad-ukbb/broad.freeze_4/temp/array_variant_concordance_callrate_0.95_non_ref_con_0.9.ht'),
                       'ukbb_array_con_common':hl.read_table('gs://broad-ukbb/broad.freeze_4/temp/array_variant_concordance_callrate_0.95_non_ref_con_0.9_AF_0.001.ht')})

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


def generate_sibling_singletons(mt, relatedness_ht, num_var_per_sibs_cutoff=None):
    mt = filter_to_autosomes(mt)
    relatedness_ht = annotate_relationship(relatedness_ht)
    relatedness_siblings_ht = relatedness_ht.filter(relatedness_ht.relationship_classification == 'Full-sibling')
    relatedness_siblings_i = relatedness_siblings_ht.i.s.collect()
    relatedness_siblings_j = relatedness_siblings_ht.j.s.collect()

    relatedness_siblings_dict = {}
    num_sibling = 0
    already_in = set([])
    sibling_to_filter = set([])
    for i,j in zip(relatedness_siblings_i,relatedness_siblings_j):
        if i not in already_in and j not in already_in:
            already_in.add(i)
            already_in.add(j)
            relatedness_siblings_dict[i] = num_sibling
            relatedness_siblings_dict[j] = num_sibling
            num_sibling += 1
        elif i in already_in and j not in already_in:
            sibling_to_filter.add(j)
        elif j in already_in and i not in already_in:
            sibling_to_filter.add(i)

    relatedness_siblings_dict = hl.literal(relatedness_siblings_dict)

    # This will filter the MT to only the calls with two non ref genotypes since if it is in both
    # siblings there will be 2 and we only want singletons so we can next check that both siblings have the
    # Variant
    mt = mt.filter_cols(~hl.literal(sibling_to_filter).contains(mt.s))
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.n_non_ref == 2)

    mt = mt.annotate_cols(sibling_index=relatedness_siblings_dict.get(mt.s))
    mt = mt.filter_cols(hl.is_defined(mt.sibling_index))

    dataset_result = (mt.group_cols_by(mt.sibling_index)
                      .aggregate_cols(sibling_pair=hl.agg.collect_as_set(mt.s))
                      .aggregate_entries(sibling_n_alt=hl.agg.sum(mt.GT.n_alt_alleles()),
                                         sibling_non_ref=hl.agg.count_where(mt.GT.is_non_ref()))
                      .result())
    dataset_result = dataset_result.filter_entries((dataset_result.sibling_non_ref == 2) &
                                                   (dataset_result.sibling_n_alt == 2))
    dataset_result = dataset_result.annotate_rows(
        sibling_singleton=hl.agg.count_where(hl.is_defined(dataset_result.sibling_non_ref)))
    dataset_result = dataset_result.filter_rows(dataset_result.sibling_singleton == 1)
    dataset_result = dataset_result.annotate_cols(num_variants=hl.agg.count_where(dataset_result.sibling_non_ref == 2))

    if num_var_per_sibs_cutoff:
        dataset_result = dataset_result.filter_cols(dataset_result.num_variants <= num_var_per_sibs_cutoff)
        dataset_result = dataset_result.annotate_rows(
            sibling_singleton=hl.agg.count_where(hl.is_defined(dataset_result.sibling_non_ref)))
        dataset_result = dataset_result.filter_rows(dataset_result.sibling_singleton == 1)

    return dataset_result


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
    hl._set_flags(newaggs=None)
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
        mt = get_ukbb_data(data_source, freeze, non_refs_only=True, meta_root='meta')
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

    if args.generate_sibling_singletons:
        relatedness_ht = hl.read_table(relatedness_ht_path(data_source, freeze))
        mt = get_ukbb_data(data_source, freeze, split=True)
        sib_mt = generate_sibling_singletons(mt, relatedness_ht, args.num_var_per_sibs_cutoff)
        sib_mt.rows().naive_coalesce(500).write(var_annotations_ht_path(data_source, freeze, 'sib_singletons'), overwrite=args.overwrite)

    # if args.generate_de_novos:  # (2.2 min/part @ 100K = 3K CPU-hours) + (7.4 m/p = 12K) + (34 m/p = ~44K) = 59K
    #     # Turn on spark speculation?
    #     mt = get_gnomad_data(data_type, raw=True, split=False)
    #     freq_data = hl.read_table(annotations_ht_path(data_type, 'family_stats'))
    #     mt = generate_de_novos(mt, fam_path(data_type), freq_data)
    #     mt.write(annotations_ht_path(data_type, 'de_novos'), args.overwrite)

    # TODO: fix resources and hardcoded paths
    if args.create_truth_data:
        # Split truth data set randomly
        sib_ht = hl.read_table(var_annotations_ht_path(data_source, freeze, 'sib_singletons'))
        sib_ht_test = sib_ht.sample(p=args.test_train_split)
        sib_ht_train = sib_ht.anti_join(sib_ht_test)
        logger.info(f'Keeping {sib_ht_test.count()} variants for testing and {sib_ht_train.count()} variants for training out of {sib_ht.count()} total variants')
        sib_ht_test.write(var_annotations_ht_path(data_source, freeze, 'sib_singletons.test'), overwrite=args.overwrite)
        sib_ht_train.write(var_annotations_ht_path(data_source, freeze, 'sib_singletons.train'), overwrite=args.overwrite)

    if args.annotate_truth_data:
        mt = get_ukbb_data(data_source, freeze, meta_root=None)
        mt = annotate_truth_data(mt).checkpoint(var_annotations_ht_path(data_source, freeze, 'truth_data'), overwrite=args.overwrite)
        mt.summarize()


# TODO: add groupings
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('-o','--overwrite', help='Overwrite data', action='store_true')

    # parser.add_argument('--vep', help='Runs VEP', action='store_true')
    parser.add_argument('--generate_allele_data', help='Calculates allele data', action='store_true')
    parser.add_argument('--generate_qc_annotations', help='Calculates QC annotations', action='store_true')
    # parser.add_argument('--generate_qual_hists', help='Calculates GQ, DP, AB histograms per variant', action='store_true')
    parser.add_argument('--generate_call_stats', help='Calculates call stats', action='store_true')
    parser.add_argument('--generate_family_stats', help='Calculates family stats', action='store_true')
    parser.add_argument('--generate_sibling_singletons', help='Creates a hail Table of variants that are sibling singletons', action='store_true')
    parser.add_argument('--num_var_per_sibs_cutoff', help='Percentage of truth data to hold back for testing', default=40)
    parser.add_argument('--include_adj_family_stats', help='Also calculate family stats for adj genotypes', action='store_true')
    # parser.add_argument('--generate_de_novos', help='Calculates de novo data', action='store_true')
    parser.add_argument('--create_truth_data', help='Create additional UKBB truth data by selecting sibling singletons and array-concordant variants', action='store_true')
    parser.add_argument('--test_train_split', help='Percentage of truth data to hold back for testing', default=0.2)

    parser.add_argument('--annotate_truth_data', help='Annotates MT with truth data', action='store_true')

    parser.add_argument('--calculate_medians', help='Calculate metric medians (warning: slow)', action='store_true')
    parser.add_argument('--calculate_all_annotations', help='Calculation many more annotations (warning: slow)', action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)

