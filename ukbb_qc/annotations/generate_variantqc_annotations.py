from gnomad_hail import *
from gnomad_qc.annotations.generate_qc_annotations import *
import gnomad_hail.resources.basics as gres
from gnomad_hail.utils.generic import vep_or_lookup_vep, vep_struct_to_csq
from ukbb_qc.resources.resources import *
from ukbb_qc.utils.utils import *
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


def annotate_truth_data(ht: hl.Table, truth_tables: Dict[str, hl.Table]) -> hl.Table:
    """
    Writes bi-allelic sites MT with the following annotations by default:
     - Omni
     - HapMap
     - 1KG high conf SNVs
     - Mills

    Additional truthset tables can be passed with `truth_tables`

    :param ht: Table to annotate
    :param truth_tables: A dictionary containing additional truth set hail tables keyed by the annotation to use in HT
    :return: Table with qc annotations
    """
    truth_ht =  hl.read_table(truth_ht_path)
    ht = ht.join(truth_ht, how="left")

    ht = ht.annotate(
        **{
            root: hl.is_defined(truth[ht.key]) for root, truth in truth_tables.items()
        }
    )

    return ht


def generate_sibling_singletons(mt, relatedness_ht, num_var_per_sibs_cutoff=None):
    """
    Gets a MatrixTable of all the variants that are found in sibling pairs and not in the rest of the input MatrixTable
    Columns of the MatrixTable are sibling pairs not individuals

    :param MatrixTable mt: Full MT
    :param Table relatedness_ht: output Table from relatedness inference, contains sample pairs with kin and ibd
    :param int num_var_per_sibs_cutoff: Max number of sibling singletons in a pair for the pair to be included in truth set
    :return: MatrixTable of sibling pairs filtered to variant unique to the sibling pair (sibling singletons)
    :rtype: MatrixTable
    """
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


def generate_qual_hists(mt: hl.MatrixTable) -> hl.Table:
    #mt = hl.split_multi_hts(mt)
    return mt.annotate_rows(
        gq_hist_alt=hl.agg.filter(mt.GT.is_non_ref(), hl.agg.hist(mt.GQ, 0, 100, 20)),
        gq_hist_all=hl.agg.hist(mt.GQ, 0, 100, 20),
        dp_hist_alt=hl.agg.filter(mt.GT.is_non_ref(), hl.agg.hist(mt.DP, 0, 100, 20)),
        dp_hist_all=hl.agg.hist(mt.DP, 0, 100, 20),
        ab_hist_alt=hl.agg.filter(mt.GT.is_het(), hl.agg.hist(mt.AD[1] / hl.sum(mt.AD), 0, 1, 20))
    ).rows()


def main(args):
    hl.init(log='/generate_variantqc_annotations.log')
    data_source = args.data_source
    freeze = args.freeze

    if args.vep:  # CPU-hours: 250 (E), 600 (G)
        logger.info(f'Running VEP on the MT...')
        ht = get_ukbb_data(data_source, freeze).select_rows()
        ht = vep_or_lookup_vep(ht, reference='GRCh38').write(var_annotations_ht_path(data_source, freeze, 'vep'), args.overwrite)
        ht = ht.annotate(vep_csq=vep_struct_to_csq(ht.vep))
        ht.write(var_annotations_ht_path(data_source, freeze, 'vep'), args.overwrite)

    if args.generate_allele_data:
        mt = get_ukbb_data(data_source, freeze, split=False)
        generate_allele_data(mt).write(var_annotations_ht_path(data_source, freeze, 'allele_data'), overwrite=args.overwrite)

    if args.generate_qc_annotations:
        # Turn on spark speculation: --properties 'spark:spark.speculation=true,spark:spark.speculation.quantile=0.9,spark:spark.speculation.multiplier=3'
        mt = get_ukbb_data(data_source, freeze, non_refs_only=True, meta_root='meta')
        mt = generate_qc_annotations(mt, all_annotations=args.calculate_all_annotations, medians=args.calculate_medians)
        mt.write(var_annotations_ht_path(data_source, freeze, 'qc_stats'), stage_locally=True, overwrite=args.overwrite)

    if args.generate_qual_hists:  # CPU-hours: 4000 (E), 8000 (G)
        mt = get_ukbb_data(data_source, freeze, raw=True, split=False)
        mt = mt.drop('gvcf_info')
        mt = mt.key_rows_by('locus', 'alleles')
        mt = hl.experimental.sparse_split_multi(mt)
        ht = generate_qual_hists(mt)
        write_temp_gcs(ht, var_annotations_ht_path(data_source, freeze, 'qual_hists'), args.overwrite)

    if args.generate_call_stats:
        if data_source == "regeneron":
            mt = get_ukbb_data(data_source, freeze)
            meta_ht = hl.read_table(meta_ht_path(data_source, freeze))
            meta_ht = meta_ht.annotate(high_quality=~meta_ht.is_filtered)
            mt = mt.annotate_cols(**{'meta': meta_ht[mt.s]})
        else:
            mt = get_ukbb_data(data_source, freeze, meta_root='meta')
        generate_call_stats(mt).write(var_annotations_ht_path(data_source, freeze, 'call_stats'), overwrite=args.overwrite)

    if args.generate_family_stats:  # CPU-hours: 8K (E), 13K (G)
        if data_source == "regeneron":
            mt = get_ukbb_data(data_source, freeze)
            meta_ht = hl.read_table(meta_ht_path(data_source, freeze))
            meta_ht = meta_ht.annotate(high_quality=~meta_ht.is_filtered)
            mt = mt.annotate_cols(**{'meta': meta_ht[mt.s]})
        else:
            mt = get_ukbb_data(data_source, freeze, meta_root='meta')
        ht, sample_table = generate_family_stats(mt, inferred_ped_path(data_source, freeze), args.include_adj_family_stats)
        ht.write(var_annotations_ht_path(data_source, freeze, 'family_stats'), stage_locally=True, overwrite=args.overwrite)
        sample_table.write(sample_annotations_table_path(data_source, freeze, 'family_stats'), stage_locally=True, overwrite=args.overwrite)

    if args.generate_sibling_singletons:
        relatedness_ht = hl.read_table(relatedness_ht_path(data_source, freeze))
        mt = get_ukbb_data(data_source, freeze)
        sib_mt = generate_sibling_singletons(mt, relatedness_ht, args.num_var_per_sibs_cutoff)
        sib_mt = sib_mt.checkpoint(get_mt_checkpoint_path(data_source, freeze, 'sib_singletons'), overwrite=args.overwrite)
        sib_mt.rows().naive_coalesce(500).write(var_annotations_ht_path(data_source, freeze, 'sib_singletons'), overwrite=args.overwrite)

        # Split sibling singleton data set randomly
        sib_ht = hl.read_table(var_annotations_ht_path(data_source, freeze, 'sib_singletons'))
        sib_ht_test = sib_ht.sample(p=args.test_train_split)
        sib_ht_train = sib_ht.anti_join(sib_ht_test)
        logger.info(f'Keeping {sib_ht_test.count()} variants for testing and {sib_ht_train.count()} variants for training out of {sib_ht.count()} total variants')
        sib_ht_test.write(var_annotations_ht_path(data_source, freeze, 'sib_singletons.test'), overwrite=args.overwrite)
        sib_ht_train.write(var_annotations_ht_path(data_source, freeze, 'sib_singletons.train'),
                           overwrite=args.overwrite)

    if args.generate_array_concordant_ht:
        variants = hl.read_table(array_variant_concordance_path(data_source, freeze))
        callrate_cutoff = variants.callrate_cutoff.take(1)[0]
        af_cutoff = variants.af_cutoff.take(1)[0]

        exome_mt = hl.read_matrix_table(
            get_mt_checkpoint_path(
                data_source,
                freeze,
                name=f"exome_subset_concordance_callrate_{callrate_cutoff}_af_{af_cutoff}"
            )
        )

        variants = variants.annotate(AF=exome_mt.rows()[variants.key].variant_qc.AF[1])
        variants = variants.filter(
            (variants.prop_gt_con_non_ref > args.concordance_cutoff)
            & (variants.AF > args.variant_qc_af_cutoff)
        )
        variants = variants.repartition(1000)
        variants.write(
            var_annotations_ht_path(
                data_source,
                freeze,
                f'array_con_con_{args.concordance_cutoff}_AF_{args.variant_qc_af_cutoff}'
            ),
            overwrite=args.overwrite
        )

    if args.annotate_truth_data:
        ht = get_ukbb_data(data_source, freeze, meta_root=None).select_rows()
        truth_ht = annotate_truth_data(
            ht,
            {
                'sib_singletons': hl.read_table(
                    var_annotations_ht_path(data_source, freeze, 'sib_singletons.train')
                ),
                'ukbb_array_con_common': hl.read_table(
                    var_annotations_ht_path(
                        data_source,
                        freeze,
                        'array_variant_concordance_callrate_0.95_non_ref_con_0.9_AF_0.001'
                    )
                )
            }
        )

        truth_ht = truth_ht.checkpoint(
            var_annotations_ht_path(data_source, freeze, 'truth_data'),
            overwrite=args.overwrite
        )
        truth_ht.summarize()


# TODO: add groupings
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('-o','--overwrite', help='Overwrite data', action='store_true')

    parser.add_argument('--vep', help='Runs VEP', action='store_true')
    parser.add_argument('--generate_allele_data', help='Calculates allele data', action='store_true')
    parser.add_argument('--generate_qc_annotations', help='Calculates QC annotations', action='store_true')
    parser.add_argument('--generate_qual_hists', help='Calculates GQ, DP, AB histograms per variant', action='store_true')
    parser.add_argument('--generate_call_stats', help='Calculates call stats', action='store_true')
    parser.add_argument('--generate_family_stats', help='Calculates family stats', action='store_true')
    parser.add_argument('--generate_sibling_singletons', help='Creates a hail Table of variants that are sibling singletons', action='store_true')
    parser.add_argument('--generate_array_concordant_ht', help='Creates a hail Table of array concordant variants', action='store_true')
    parser.add_argument(
        "--concordance_cutoff",
        help="Array exome concordance cutoff for variant QC HT.",
        type=float,
        default=0.9
    )
    parser.add_argument(
        "--variant_qc_af_cutoff",
        help="Allele frequency cutoff used for variant QC HT creation, must be equal to or greater than af_cutoff.",
        type=float,
        default=0.001
    )
    parser.add_argument('--num_var_per_sibs_cutoff', help='Max number of sibling singletons in a pair for the pair to be included in truth set', default=40)
    parser.add_argument('--include_adj_family_stats', help='Also calculate family stats for adj genotypes', action='store_true')
    parser.add_argument('--test_train_split', help='Percentage of truth data to hold back for testing', default=0.2)
    parser.add_argument('--annotate_truth_data', help='Creates a HT of UKBB variants annotated with truth sites', action='store_true')
    parser.add_argument('--calculate_medians', help='Calculate metric medians (warning: slow)', action='store_true')
    parser.add_argument('--calculate_all_annotations', help='Calculation many more annotations (warning: slow)', action='store_true')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)

