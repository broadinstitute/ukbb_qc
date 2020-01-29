from variant_qc.variant_qc_functions import *

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("make_fake_pr_curves.log")
logger.setLevel(logging.INFO)


# TODO: compute sites-level concordance between Broad NA12878 and Regeneron NA12878
def write_truth_concordance(data_source: str, freeze: int, truth_sample: str, overwrite: bool) -> None:
    sample_mapping = {
        'NA12878': 'Coriell_NA12878_NA12878'
    }

    truth_mt = hl.read_matrix_table('gs://broad-ukbb/broad.freeze_4/data/NovaSeq_Coriell_Freeze.NF.truth_mt')
    truth_mt = truth_mt.filter_cols(truth_mt.s == sample_mapping[truth_sample])
    truth_mt = unphase_truth_mt(hl.split_multi_hts(truth_mt, left_aligned=False))
    truth_mt = truth_mt.key_cols_by(s=hl.str(truth_sample))
    truth_mt = truth_mt.repartition(1000 if data_type == 'genomes' else 100, shuffle=False)

    broad_ht = hl.read_table('gs://broad-ukbb/broad.freeze_4/temp/joint_broad_regeneron.NA12878.ht')
    broad_mt = hl.read_matrix_table(hardcalls_mt_path(data_source, freeze, split=True))
    # TODO: once Broad callset has native NA12878, filter instead to NA12878 Broad calls
    # broad_mt = broad_mt.key_cols_by(s=hl.str(truth_sample))

    # truth_mt = hl.read_matrix_table(NA12878_mt_path() if truth_sample == 'NA12878' else syndip_mt_path())
    # truth_mt = truth_mt.key_cols_by(s=hl.str(truth_sample))
    # if data_type == 'exomes':
    #     exome_calling_intervals = hl.import_locus_intervals(exome_calling_intervals_path, skip_invalid_intervals=True)
    #     truth_mt = truth_mt.filter_rows(hl.is_defined(exome_calling_intervals[truth_mt.locus]))
    # truth_mt = unphase_mt(hl.split_multi_hts(truth_mt, left_aligned=False))
    #
    # sample_concordance_ht, sites_concordance_ht = compute_concordance(mt, truth_mt, name=truth_sample)
    # sites_concordance_ht.write(annotations_ht_path(data_type, f'{truth_sample}_concordance'), overwrite=overwrite)
    # sample_concordance_ht.write(sample_annotations_table_path(data_type, f'{truth_sample}_concordance'), overwrite=overwrite)



def create_binned_concordance(data_source: str, freeze: int, metric: str, nbins: int, overwrite: bool, truth_sample: str = 'NA12878') -> None:
    """
    Creates and writes a concordance table binned by rank (both absolute and relative) for a given data type, truth sample and metric.

    :param str data_source: 'broad' or 'regeneron'
    :param int freeze: UKBB tranche version
    :param str truth_sample: Which truth sample concordance to load
    :param str metric: One of the evaluation metrics (or a RF hash)
    :param int nbins: Number of bins for the rank
    :param bool overwrite: Whether to overwrite existing table
    :return: Nothing -- just writes the table
    :rtype: None
    """

    if hl.hadoop_exists(binned_concordance_path(data_source, freeze, truth_sample, metric)+'/_SUCCESS') and not overwrite:
        logger.warn(f"Skipping binned concordance creation as {binned_concordance_path(data_source, freeze, truth_sample, metric)} exists and overwrite=False")
    else:
        ht = hl.read_table(f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/joint_broad_regeneron.NA12878.ht')
        lcr = hl.import_locus_intervals(lcr_intervals_path, reference_genome='GRCh38')
        count1 = ht.count()
        # high_conf_intervals = hl.import_locus_intervals(NA12878_high_conf_regions_bed_path)  # Check 38 paths
        # segdup = hl.import_locus_intervals(segdup_intervals_path)  # Check 38 paths
        # ht = ht.filter(
        #     hl.is_missing(lcr[ht.locus])  # TODO: add in 38 paths as they become available
        # )
        logger.info(f'Found {ht.count()} variants out of {count1} after removing LCR regions')

        if metric in ['vqsr', 'rf_2.0.2', 'rf_2.0.2_beta', 'cnn']:
            metric_ht = hl.read_table(score_ranking_path(data_type, metric))
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
                         tp=ht.pass_broad & ht.pass_reg,
                         # FP => allele is found only in test data set
                         fp=ht.pass_broad & ~hl.is_defined(ht.pass_reg),  # pass_reg is missing for samples that do not pass in Regeneron
                         # FN => allele is found only in truth data set
                         fn=~ht.pass_broad & ht.pass_reg
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

        ht.write(binned_concordance_path(data_source, freeze, truth_sample, metric), overwrite=overwrite)


def get_binned_concordance_pd(data_source: str, freeze: int, truth_samples: List[str], models: Union[Dict[str, str], List[str]]) -> pd.DataFrame:
    """
    Creates a pandas DF containing the binned concordance results for all given truth samples / models.

    :param str data_source: 'broad' or 'regeneron'
    :param int freeze: UKBB tranche version
    :param list of str truth_samples: List of truth samples to include
    :param list of str or dict of str -> str models: Models to include. Either a list of the model ids, or a dict with model id -> model name for display
    :return: Pandas dataframe with binned concordance results
    :rtype: DataFrame
    """

    def get_binned_concordance_ht(data_source: str, freeze:int, truth_samples: List[str], models: Dict[str, str]) -> hl.Table:
        """
        Combines binned concordance results for multiple truth samples and/or models into a single Table.
        """
        hts = []
        for truth_sample in truth_samples:
            for model_id, model_name in models.items():
                ht = hl.read_table(binned_concordance_path(data_source, freeze, truth_sample, model_id))
                ht = ht.annotate(truth_sample=truth_sample, model=model_name)
                hts.append(ht)

        return hts[0].union(*hts[1:])

    def compute_cumul_metrics(df: pd.DataFrame) -> pd.DataFrame:
        """
        Computes cumulative metrics on a pandas DF.
        """
        df = df.sort_values(by=['bin'])
        df['cum_tp'] = df['tp'].cumsum()
        df['cum_fp'] = df['fp'].cumsum()
        total_pos = df['tp'].sum() + df['fn'].sum()
        total_neg = df['fp'].sum()
        df['cum_tn'] = total_neg - df['cum_fp']
        df['cum_fn'] = total_pos - df['cum_tp']
        df['precision'] = df['cum_tp'] / (df['cum_tp'] + df['cum_fp'])
        df['recall'] = df['cum_tp'] / (df['cum_tp'] + df['cum_fn'])
        df['cum_alleles'] = df['n_alleles'].cumsum()
        return df[['bin', 'min_score', 'max_score', 'n_alleles', 'tp', 'fp', 'fn', 'cum_alleles', 'cum_tp', 'cum_fp', 'cum_fn', 'cum_tn', 'precision', 'recall']]

    if not isinstance(models, dict):
        models = {m: m for m in models}

    df = get_binned_concordance_ht(data_source, freeze, truth_samples, models).to_pandas()
    df = df.groupby(['rank_name', 'truth_sample', 'model', 'snv']).apply(compute_cumul_metrics)
    return df.fillna(-1).groupby(['rank_name', 'truth_sample', 'model', 'snv'])


def main(args):
    hl._set_flags(newaggs=None)
    data_source = args.data_source
    freeze = args.freeze

    if args.import_vcfs:
        rg = hl.get_reference('GRCh37')
        grch37_contigs = [x for x in rg.contigs if not x.startswith('GL') and not x.startswith('M')]
        contig_dict = dict(zip(grch37_contigs, ['chr' + x for x in grch37_contigs]))

        gl_mt = hl.import_vcf('gs://broad-ukbb/broad.freeze_4/data/vcfs/NovaSeq_Coriell_Freeze.GL.vcf.gz', force_bgz=True, reference_genome='GRCh38', array_elements_required=False, contig_recoding=contig_dict)
        gl_mt = gl_mt.checkpoint('gs://broad-ukbb/broad.freeze_4/data/NovaSeq_Coriell_Freeze.GL.mt', overwrite=args.overwrite)
        gl_mt.describe()
        gl_mt.summarize(entries=False)
        nf_mt = hl.import_vcf('gs://broad-ukbb/broad.freeze_4/data/vcfs/NovaSeq_Coriell_Freeze.NF.vcf.gz', force_bgz=True, reference_genome='GRCh38', array_elements_required=False, contig_recoding=contig_dict)
        nf_mt = nf_mt.checkpoint('gs://broad-ukbb/broad.freeze_4/data/NovaSeq_Coriell_Freeze.NF.mt', overwrite=args.overwrite)
        nf_mt.describe()
        nf_mt.summarize(entries=False)

    if args.match_raw_sites:
        broad_ht = hl.read_table(rf_path(data_source, freeze, 'rf_result', run_hash=args.run_hash))
        broad_ht = broad_ht.annotate(pass_broad=(broad_ht.rf_probability.get('TP') > 0.03))
        nf_mt = hl.read_matrix_table('gs://broad-ukbb/broad.freeze_4/data/NovaSeq_Coriell_Freeze.NF.mt')
        gl_ht = hl.read_matrix_table('gs://broad-ukbb/broad.freeze_4/data/NovaSeq_Coriell_Freeze.GL.mt').rows()
        gl_ht = gl_ht.annotate(pass_reg=True)

        regeneron_mt = nf_mt.filter_cols(nf_mt.s == 'Coriell_NA12878_NA12878')
        regeneron_ht = regeneron_mt.filter_rows(hl.agg.sum(regeneron_mt.GT.n_alt_alleles()) > 0, keep=True).rows()
        logger.info(f'Found {regeneron_ht.count()} variants in NA12878 out of {regeneron_mt.count_rows()} in callset')
        regeneron_ht = regeneron_ht.join(gl_ht, how='left')

        joint_ht = regeneron_ht.join(broad_ht, how='inner').checkpoint(f'gs://broad-ukbb/{data_source}.freeze_{freeze}/temp/joint_broad_regeneron.NA12878.ht', overwrite=args.overwrite)
        filter_summary = joint_ht.group_by(joint_ht.pass_broad, joint_ht.pass_reg).aggregate(n=hl.agg.count())
        filter_summary.show()

    if args.compute_binned_concordance:
        create_binned_concordance(data_source, freeze, args.run_hash, args.n_bins, overwrite=args.overwrite)


    # Make tables for plotting:
    # +1. Modify create_binned_concordance
    # +2. Run get_binned_concordance_pd out of the box (variant_qc_plots.py)
    # +3. Run plot_concordance_pr (as in example)
    # 4. Add plots to compare against gnomAD, other RF model

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Current freeze #', default=CURRENT_FREEZE, type=int)
    parser.add_argument('--run_hash', help='Run hash for RF results to be ranked.')
    parser.add_argument('--n_bins', help='Number of bins for the binned file (default: 100)', default=100, type=int)

    parser.add_argument('--import_vcfs', help='Import NA12878 VCFs', action='store_true')
    parser.add_argument('--match_raw_sites', help='Compare filter status in Broad and Regeneron NA12878 VCFs', action='store_true')
    parser.add_argument('--compute_binned_concordance', help='Compute site concordance between Broad and Regeneron NA12878 calls with binned ranking annotations', action='store_true')

    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()


    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
