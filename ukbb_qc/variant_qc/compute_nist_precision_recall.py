from variant_qc.variant_qc_functions import *

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("make_fake_pr_curves.log")
logger.setLevel(logging.INFO)

# TODO: fix partitioning for all input MTs/HTs
# TODO: compute sites-level concordance between Broad NA12878 and NIST NA12878 once NA12878 is joint-called in Broad:
# write_truth_concordance()

truth_samples = {
    'HG001': 'HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_all.vcf.gz',
    'HG002': 'HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_all.vcf.gz',
    'HG003': 'HG003_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_all.vcf.gz',
    'HG004': 'HG004_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC_CHROM1-22_v.3.3.2_all.vcf.gz',
    'HG005': 'HG005_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_all.vcf.gz',
    'HG006': 'HG006_GIAB_GRCh38_highconf_CG-IllFB-IllSNT-10X_testing_CHROM1-22_v.3.3.2_all.vcf.gz',
    'HG007': 'HG007_GIAB_GRCh38_highconf_CG-IllFB-IllSNT-10X_CHROM1-22_v.3.3.2_all.vcf.gz'
}


def pre_process_mt(mt_path: str) -> hl.Table:
    '''
    Splits MT, min reps variants, and creates HT

    :param mt_path:
    :return:
    :rtype: Table
    '''
    mt = hl.read_matrix_table(mt_path)
    mt = hl.split_multi_hts(mt)
    mt = mt.key_rows_by(**hl.min_rep(mt.locus, mt.alleles))
    return mt.rows()


def create_regeneron_ht(data_source, freeze):
    nf_mt = get_ukbb_data(data_source, freeze, adj=True, split=True)
    nf_ht = nf_mt.annotate_rows(singleton=hl.agg.sum(nf_mt.GT.n_alt_alleles()) == 1).rows()
    gl_ht = pre_process_mt(f'gs://broad-ukbb/regeneron.freeze_{freeze}/data/regeneron.freeze_{freeze}.gl.mt')
    ht = nf_ht.annotate(pass_regeneron=hl.is_defined(gl_ht.index(nf_ht.key)), score=nf_ht.qual)

    return ht.checkpoint(f'{variant_qc_prefix(data_source, freeze)}/assessment/regeneron.freeze_{freeze}.nf.gl_annotated.ht', overwrite=args.overwrite)


def create_binned_concordance(data_source: str, freeze: int, metric: str, nbins: int, overwrite: bool, truth_sample: str, filter_difficult_regions: bool) -> None:
    """
    Creates and writes a concordance table binned by rank (both absolute and relative) for a given data type, truth sample and metric.

    :param str data_source: 'broad' or 'regeneron'
    :param int freeze: UKBB tranche version
    :param str truth_sample: Which truth sample concordance to load
    :param str metric: One of the evaluation metrics (or a RF hash)
    :param int nbins: Number of bins for the rank
    :param bool filter_difficult_regions: Remove variants annotated in NIST dataset as occurring at least 1 "difficult region" track
    :param bool overwrite: Whether to overwrite existing table
    :return: Nothing -- just writes the table
    :rtype: None
    """

    if hl.hadoop_exists(binned_concordance_path(data_source, freeze, truth_sample, metric)+'/_SUCCESS') and not overwrite:
        logger.warn(f"Skipping binned concordance creation as {binned_concordance_path(data_source, freeze, truth_sample, metric)} exists and overwrite=False")
    else:
        ht = hl.read_table(f'{variant_qc_prefix(data_source, freeze)}/assessment/joint_nist_{data_source}.{truth_sample}.{args.run_hash}.ht')
        if filter_difficult_regions:
            ht = ht.filter(hl.is_defined(ht.info.difficultregion), keep=False)

        if metric == 'vqsr':
            metric_ht = hl.read_table(score_ranking_path(data_source, freeze, metric))
        elif metric == 'gnomad_exomes':
            metric_ht = hl.read_table('gs://gnomad-public/release/2.1.1/liftover_grch38/ht/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht')
        elif metric == 'regeneron':
            metric_ht = hl.read_table(f'{variant_qc_prefix(data_source, freeze)}/assessment/regeneron.freeze_{freeze}.nf.gl_annotated.ht')
            metric_ht = add_rank(metric_ht, -1.0 * metric_ht.score)
        else:
            metric_ht = hl.read_table(rf_path(data_source, freeze, 'rf_result', run_hash=metric))

        metric_snvs, metrics_indels = metric_ht.aggregate([hl.agg.count_where(hl.is_snp(metric_ht.alleles[0], metric_ht.alleles[1])),
                                                           hl.agg.count_where(~hl.is_snp(metric_ht.alleles[0], metric_ht.alleles[1]))])

        snvs, indels = ht.aggregate([hl.agg.count_where(hl.is_snp(ht.alleles[0], ht.alleles[1])),
                                     hl.agg.count_where(~hl.is_snp(ht.alleles[0], ht.alleles[1]))])

        ht = ht.annotate_globals(global_counts=hl.struct(snvs=metric_snvs, indels=metrics_indels),
                                 counts=hl.struct(snvs=snvs, indels=indels))

        ht = ht.annotate(snv=hl.is_snp(ht.alleles[0], ht.alleles[1]),
                         score=metric_ht[ht.key].score,
                         global_rank=metric_ht[ht.key].rank,
                         # TP => allele is found in both data sets
                         tp=ht[f'pass_{data_source}'] & ht.pass_nist,
                         # FP => allele is found only in test data set
                         fp=ht[f'pass_{data_source}'] & ~ht.pass_nist,
                         # FN => allele is found only in truth data set
                         fn=~ht[f'pass_{data_source}'] & ht.pass_nist
                         )

        # TODO: replace above block with the following when sites-concordance is available
        # ht = ht.annotate(snv=hl.is_snp(ht.alleles[0], ht.alleles[1]),
        #                  score=metric_ht[ht.key].score,
        #                  global_rank=metric_ht[ht.key].rank,
        #                  # TP => allele is found in both data sets
        #                  n_tp=ht.concordance[3][3] + ht.concordance[3][4] + ht.concordance[4][3] + ht.concordance[4][4],
        #                  # FP => allele is found only in test data set
        #                  n_fp=hl.sum(ht.concordance[3][:2]) + hl.sum(ht.concordance[4][:2]),
        #                  # FN => allele is found only in truth data set
        #                  n_fn=hl.sum(ht.concordance[:2].map(lambda x: x[3] + x[4]))
        #                  )

        ht = add_rank(ht, -1.0*ht.score)

        ht = ht.annotate(rank=[
            hl.tuple(['global_rank', (ht.global_rank + 1) / hl.cond(
                             ht.snv,
                             ht.globals.global_counts.snvs,
                             ht.globals.global_counts.indels
                         )]),
            hl.tuple(['truth_sample_rank', (ht.rank + 1) / hl.cond(
                             ht.snv,
                             ht.globals.counts.snvs,
                             ht.globals.counts.indels
                         )])
        ])

        ht = ht.explode(ht.rank)
        ht = ht.annotate(
            rank_name=ht.rank[0],
            bin=hl.int(ht.rank[1] * nbins)
        )

        ht = ht.group_by(
            'rank_name',
            'snv',
            'bin'
        ).aggregate(
            # Look at site-level metrics -> tp > fp > fn -- only important for multi-sample comparisons
            tp=hl.agg.count_where(ht.tp),
            fp=hl.agg.count_where(ht.fp),
            fn=hl.agg.count_where(ht.fn),
            min_score=hl.agg.min(ht.score),
            max_score=hl.agg.max(ht.score),
            n_alleles=hl.agg.count()
        ).repartition(5)

        metric = metric + '.hc_regions' if filter_difficult_regions else metric
        ht.write(binned_concordance_path(data_source, freeze, truth_sample, metric), overwrite=overwrite)


def main(args):
    hl._set_flags(newaggs=None)
    data_source = args.data_source
    freeze = args.freeze

    if args.import_truth_vcfs:
        for sample_id, file_handle in truth_samples.items():
            mt = hl.import_vcf(f'gs://broad-ukbb/resources/nist/{file_handle}', force_bgz=True,
                              reference_genome='GRCh38', header_file='gs://broad-ukbb/resources/nist/nist_header_for_import.txt')
            mt = mt.checkpoint(f'gs://broad-ukbb/resources/nist/{sample_id}.mt', overwrite=args.overwrite)
            mt.summarize()
            # NOTE: splitting multi, min-repping directly on VCF import, and writing doesn't work -- errors out;
            # requires a MT write without splitting and then splitting multi on HT
            # This is likely because of invalid PLs in the NIST genotypes

            mt = hl.read_matrix_table(f'gs://broad-ukbb/resources/nist/{sample_id}.mt')
            mt = hl.split_multi_hts(mt)
            mt = mt.key_rows_by(**hl.min_rep(mt.locus, mt.alleles))

            ht = mt.rows()
            ht = ht.annotate(pass_nist=hl.len(ht.filters) == 0)
            ht.summarize()
            ht.write(f'gs://broad-ukbb/resources/nist/{sample_id}.ht', overwrite=args.overwrite)


    if args.match_raw_sites:
        for sample_id in truth_samples.keys():
            nist_ht = hl.read_table(f'gs://broad-ukbb/resources/nist/{sample_id}.ht')
            if args.run_hash == 'vqsr':
                ht = hl.read_table(score_ranking_path(data_source, freeze, args.run_hash))
                ht = ht.annotate(pass_broad=hl.len(ht.filters) == 0)
            elif args.run_hash == 'gnomad_exomes':
                ht = hl.read_table('gs://gnomad-public/release/2.1.1/liftover_grch38/ht/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht')
                ht = ht.select('was_split', 'singleton', 'score', 'filters')
                ht = ht.annotate(pass_broad=hl.len(ht.filters) == 0)
                # NOTE: some sites will be missing annotations (filters, was_split, etc.) because they were not successfully lifted over (ReverseComplementedAlleles)
            elif args.run_hash == 'regeneron':
                if not hl.hadoop_exists(f'{variant_qc_prefix(data_source, freeze)}/assessment/regeneron.freeze_{freeze}.nf.gl_annotated.ht/_SUCCESS'):
                    ht = create_regeneron_ht('regeneron', freeze)
                else:
                    ht = hl.read_table(f'{variant_qc_prefix(data_source, freeze)}/assessment/regeneron.freeze_{freeze}.nf.gl_annotated.ht')
            else:
                ht = hl.read_table(rf_path(data_source, freeze, 'rf_result', run_hash=args.run_hash))
                ht = ht.annotate(pass_broad=(ht.rf_probability.get('TP') > args.tp_threshold))

            # Match sites on locus and allele
            joint_ht = (nist_ht.join(ht, how='inner').naive_coalesce(500)
                        .checkpoint(f'{variant_qc_prefix(data_source, freeze)}/assessment/joint_nist_{data_source}.{sample_id}.{args.run_hash}.ht', overwrite=args.overwrite))
            filter_summary = joint_ht.group_by(joint_ht.pass_nist, joint_ht[f'pass_{data_source}']).aggregate(n=hl.agg.count())
            filter_summary.show()

            joint_ht = joint_ht.annotate(is_snp=hl.is_snp(joint_ht.alleles[0], joint_ht.alleles[1]))
            filter_summary = joint_ht.group_by(joint_ht.is_snp, joint_ht.was_split_1,
                                               joint_ht.pass_nist, joint_ht[f'pass_{data_source}']).aggregate(n=hl.agg.count())
            filter_summary.show(20)

    if args.filter_difficult_regions:
        for sample_id in truth_samples.keys():
            joint_ht = hl.read_table(f'{variant_qc_prefix(data_source, freeze)}/assessment/joint_nist_{data_source}.{sample_id}.{args.run_hash}.ht')
            joint_ht = joint_ht.filter(hl.is_defined(joint_ht.info.difficultregion), keep=False)

            filter_summary = joint_ht.group_by(joint_ht.pass_nist, joint_ht[f'pass_{data_source}']).aggregate(
                n=hl.agg.count())
            filter_summary.show()

            joint_ht = joint_ht.annotate(is_snp=hl.is_snp(joint_ht.alleles[0], joint_ht.alleles[1]))
            filter_summary = joint_ht.group_by(joint_ht.is_snp, joint_ht.was_split_1,
                                               joint_ht.pass_nist, joint_ht[f'pass_{data_source}']).aggregate(
                n=hl.agg.count())
            filter_summary.show(20)

    if args.compute_binned_concordance:
        for sample_id in truth_samples.keys():
            create_binned_concordance(data_source, freeze, args.run_hash, args.n_bins, overwrite=args.overwrite,
                                      truth_sample=sample_id, filter_difficult_regions=args.filter_difficult_regions)

# TODO: make the filtering of "difficult" region variants required for future tranche runs


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Current freeze #', default=CURRENT_FREEZE, type=int)
    parser.add_argument('--run_hash', help='Run hash for RF results to be ranked.')
    parser.add_argument('--n_bins', help='Number of bins for the binned file (default: 100)', default=100, type=int)
    parser.add_argument('--tp_threshold', help='Probability threshold for random forest filtering', type=float)

    parser.add_argument('--import_truth_vcfs', help='Import NIST VCFs', action='store_true')
    parser.add_argument('--match_raw_sites', help='Compare filter status in Broad and NIST VCFs', action='store_true')
    parser.add_argument('--compute_binned_concordance', help='Compute site concordance between Broad and NIST calls with binned ranking annotations', action='store_true')
    parser.add_argument('--filter_difficult_regions', help='Remove NIST variants that occur in at least one "difficult region" track', action='store_true')

    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()


    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
