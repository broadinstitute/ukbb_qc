from gnomad_hail import *
from gnomad_hail.utils import rf
from gnomad_qc.variant_qc.variantqc import sample_rf_training_examples
from ukbb_qc.resources import *
from ukbb_qc.variant_qc_functions import *
from pprint import pformat
import json

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("assess_variant_qc_results.log")
logger.setLevel(logging.INFO)


def main(args):
    hl._set_flags(newaggs=None)
    data_source = 'broad'
    freeze = args.freeze
    run_hash = args.run_hash


    if args.join_tables:
        gnomad_ht = hl.read_table('gs://gnomad-public/release/2.1.1/liftover_grch38/ht/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht')
        # print(hl.eval(gnomad_ht.freq_meta))
        gnomad_ht = gnomad_ht.annotate(pass_gnomad=(hl.len(gnomad_ht.filters) == 0))
        ukbb_ht = hl.read_table(rf_path(data_source, freeze, 'rf_result', run_hash=run_hash))
        ukbb_ht = ukbb_ht.annotate(pass_broad=(ukbb_ht.rf_probability.get('TP') > args.tp_threshold))
        regeneron_ht = hl.read_matrix_table(f'gs://broad-ukbb/regeneron.freeze_{freeze}/data/regeneron.freeze_{freeze}.nf.mt').rows()
        regeneron_filt_ht = hl.read_matrix_table(f'gs://broad-ukbb/regeneron.freeze_{freeze}/data/regeneron.freeze_{freeze}.gl.mt').rows()
        regeneron_ht = regeneron_ht.annotate(pass_regeneron=hl.is_defined(regeneron_filt_ht.index(regeneron_ht.key)))

        joint_ht = gnomad_ht.join(ukbb_ht, how='outer')
        joint_ht = joint_ht.join(regeneron_ht, how='outer')

        # Filter to properly covered intervals
        coverage_ht = hl.read_table(f'{sample_qc_prefix(data_source, freeze)}/interval_qc/coverage_by_target.ht')
        bad_intervals_ht = coverage_ht.filter(coverage_ht.pct_samples_20x < 0.85, keep=True).key_by('locus')
        count1 = joint_ht.count()
        logger.info(f'Found {bad_intervals_ht.count()} low-coverage intervals out of {coverage_ht.count()} total intervals')
        joint_ht = joint_ht.filter(hl.is_defined(bad_intervals_ht[joint_ht.locus]), keep=False)
        logger.info(f'Kept {joint_ht.count()} variants out of {count1} variants after removing poorly covered intervals in Broad callset')

        # nas = joint_ht.filter(hl.is_missing(joint_ht.pass_ukbb) & hl.is_missing(joint_ht.pass_gnomad))
        # nas.summarize()
        # NOTE: ~60k variants found without pass status in either callset -- these are variants lost in liftover from exome hg19 data

        joint_ht = joint_ht.checkpoint(f'{variant_qc_prefix(data_source, freeze)}/assessment/joint_gnomad_broad_regeneron.{run_hash}.ht', overwrite=args.overwrite)
        filter_summary = joint_ht.group_by(joint_ht.pass_gnomad, joint_ht.pass_broad, joint_ht.pass_regeneron).aggregate(n=hl.agg.count())
        filter_summary.show(50)


    if args.bin_rf_probs_overall:
        joint_ht = hl.read_table(f'{variant_qc_prefix(data_source, freeze)}/assessment/joint_gnomad_broad_regeneron.{run_hash}.ht')
        joint_ht = (joint_ht.group_by(tp_bin = joint_ht.rf_probability_1.get('TP') // 0.01)
                    .aggregate(n_pass_gnomad=hl.agg.count_where(hl.is_defined(joint_ht.pass_gnomad) & hl.is_defined(joint_ht.pass_broad) & joint_ht.pass_gnomad),
                               n_fail_gnomad=hl.agg.count_where(hl.is_defined(joint_ht.pass_gnomad) & hl.is_defined(joint_ht.pass_broad) & ~joint_ht.pass_gnomad),
                               n_unique_broad_vs_gnomad=hl.agg.count_where(~hl.is_defined(joint_ht.pass_gnomad) & hl.is_defined(joint_ht.pass_broad)),
                               n_pass_regeneron=hl.agg.count_where(hl.is_defined(joint_ht.pass_regeneron) & hl.is_defined(joint_ht.pass_broad) & joint_ht.pass_regeneron),
                               n_fail_regeneron=hl.agg.count_where(hl.is_defined(joint_ht.pass_regeneron) & hl.is_defined(joint_ht.pass_broad) & ~joint_ht.pass_regeneron),
                               n_unique_broad_vs_regeneron=hl.agg.count_where(~hl.is_defined(joint_ht.pass_regeneron) & hl.is_defined(joint_ht.pass_broad)),
                               n=hl.agg.count()))
        joint_ht = joint_ht.checkpoint(f'{variant_qc_prefix(data_source, freeze)}/assessment/joint_gnomad_broad_regeneron.{run_hash}.probs.ht', overwrite=args.overwrite)
        joint_ht.show(50)
        joint_ht.export(f'{variant_qc_prefix(data_source, freeze)}/assessment/joint_gnomad_broad_regeneron.{run_hash}.probs.tsv.gz')


    if args.aggregate_overall:
        joint_ht = hl.read_table(f'{variant_qc_prefix(data_source, freeze)}/assessment/joint_gnomad_broad_regeneron.{run_hash}.ht')
        # Annotate by gnomAD adj AF:
        joint_ht = joint_ht.annotate(gnomad_af_bin=hl.case()
                                     .when((joint_ht.freq[0].AF <= 0.00001), "(0, 0.00001]")
                                     .when((joint_ht.freq[0].AF > 0.00001) & (joint_ht.freq[0].AF <= 0.0001), "(0.00001, 0.0001]")
                                     .when((joint_ht.freq[0].AF > 0.0001) & (joint_ht.freq[0].AF <= 0.001), "(0.0001, 0.001]")
                                     .when((joint_ht.freq[0].AF > 0.001) & (joint_ht.freq[0].AF <= 0.01), "(0.001, 0.01]")
                                     .when((joint_ht.freq[0].AF > 0.01) & (joint_ht.freq[0].AF <= 0.02), "(0.01, 0.02]")
                                     .when((joint_ht.freq[0].AF > 0.02) & (joint_ht.freq[0].AF <= 0.03), "(0.02, 0.03]")
                                     .when((joint_ht.freq[0].AF > 0.03) & (joint_ht.freq[0].AF <= 0.04), "(0.03, 0.04]")
                                     .when((joint_ht.freq[0].AF > 0.04) & (joint_ht.freq[0].AF <= 0.05), "(0.03, 0.05]")
                                     .default(hl.null(hl.tstr)))
        joint_ht = (joint_ht.group_by(joint_ht.gnomad_af_bin, joint_ht.variant_type, joint_ht.pass_gnomad, joint_ht.pass_broad, joint_ht.pass_regeneron)
                    .aggregate(n_shared_gnomad=hl.agg.count_where(hl.is_defined(joint_ht.pass_gnomad) & hl.is_defined(joint_ht.pass_broad)),
                               n_gnomad_unique=hl.agg.count_where(hl.is_defined(joint_ht.pass_gnomad) & ~hl.is_defined(joint_ht.pass_broad)),
                               n_shared_regeneron=hl.agg.count_where(hl.is_defined(joint_ht.pass_regeneron) & hl.is_defined(joint_ht.pass_broad)),
                               n_regeneron_unique=hl.agg.count_where(hl.is_defined(joint_ht.pass_regeneron) & ~hl.is_defined(joint_ht.pass_broad)),
                               n=hl.agg.count()))
        joint_ht = joint_ht.annotate(frac_shared_gnomad=joint_ht.n_shared_gnomad/joint_ht.n,
                                     frac_gnomad_unique=joint_ht.n_gnomad_unique/joint_ht.n,
                                     frac_shared_regeneron=joint_ht.n_shared_regeneron/joint_ht.n,
                                     frac_regeneron_unique=joint_ht.n_regeneron_unique/joint_ht.n)
        joint_ht = joint_ht.checkpoint(f'{variant_qc_prefix(data_source, freeze)}/assessment/joint_gnomad_broad_regeneron.{run_hash}.agg.ht', overwrite=args.overwrite)
        joint_ht.show(50)
        joint_ht.export(f'{variant_qc_prefix(data_source, freeze)}/assessment/joint_gnomad_broad_regeneron.{run_hash}.agg.txt.bgz')


    ### PICK UP HERE ###
    if args.aggregate_nfe:

        if args.filter_nfe:
            # Read in gnomAD and filter to NFE, identify those variants
            gnomad_mt = get_gnomad_data('exomes')
            rows, cols = gnomad_mt.count()
            gnomad_mt = gnomad_mt.filter_cols(gnomad_mt.meta.pop == 'nfe', keep=True)
            logger.info(f'Found {gnomad_mt.count_cols()} NFE samples out of {cols} original samples in gnomAD exomes')
            gnomad_mt = gnomad_mt.filter_rows(hl.agg.sum(gnomad_mt.GT.n_alt_alleles()) > 0, keep=True)
            logger.info(f'Found {gnomad_mt.count_rows()} variants in NFE samples out of {rows} original variants in all gnomAD exomes')
            gnomad_mt.rows().write('gs://broad-ukbb/broad.freeze_4/temp/gnomad_nfe.ht', overwrite=args.overwrite)

            # Read in UKBB and filter to NFE, identify those variants
            # this reads in split hardcalls (non-adj)
            ukbb_mt = get_ukbb_data(data_source, freeze)
            rows, cols = ukbb_mt.count()

            # Keep only NFE found with gnomAD PC project
            # hard code path until we replace meta ht with pop adj meta ht
            meta_ht = hl.read_table(f'gs://broad-ukbb/regeneron.freeze_{freeze}/sample_qc/meta_w_pop_adj.ht')
            meta_ht = meta_ht.annotate(is_nfe=(meta_ht.gnomad_pc_project_pop == 'nfe'))
            logger.info(f'NFE samples in Regeneron meta table: {meta_ht.aggregate(hl.agg.counter(meta_ht.is_nfe))}')

            ukbb_mt = ukbb_mt.annotate_cols(**{'meta': meta_ht[ukbb_mt.s]})
            ukbb_mt = ukbb_mt.filter_cols(ukbb_mt.meta.is_nfe)    
            logger.info(f'Found {ukbb_mt.count_cols()} NFE samples out of {cols} original samples in UKBB exomes')
            ukbb_mt = ukbb_mt.annotate_rows(nfe_ac=hl.agg.sum(ukbb_mt.GT.n_alt_alleles()))
            ukbb_mt = ukbb_mt.filter_rows(ukbb_mt.nfe_ac > 0, keep=True)
            logger.info(f'Found {ukbb_mt.count_rows()} variants in NFE samples out of {rows} original variants in all UKBB exomes')
            ukbb_mt.rows().write(f'gs://broad-ukbb/broad.freeze_{freeze}/temp/ukbb_nfe.ht', overwrite=args.overwrite)

        # Filter joint ht to those variants
        gnomad_ht = hl.read_table(f'gs://broad-ukbb/broad.freeze_{freeze}/temp/gnomad_nfe.ht').select()
        rows = gnomad_ht.count()
        liftover_ht = hl.read_table(
            'gs://gnomad-public/release/2.1.1/liftover_grch38/ht/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht').key_by('original_locus', 'original_alleles').rename({'locus' : 'locus_grch38', 'alleles' : 'alleles_grch38'})
        gnomad_ht = gnomad_ht.annotate(**liftover_ht[gnomad_ht.key])
        gnomad_ht = gnomad_ht.filter(hl.is_defined(gnomad_ht.locus_grch38))
        logger.info(
            f'Found {gnomad_ht.count()} variants in lifted-over NFE samples out of {rows} original variants in NFE gnomAD exomes')
        gnomad_ht = gnomad_ht.key_by('locus_grch38', 'alleles_grch38')

        ukbb_ht = hl.read_table(f'gs://broad-ukbb/broad.freeze_{freeze}/temp/ukbb_nfe.ht').select()
        nfe_ht = ukbb_ht.join(gnomad_ht, how='outer').checkpoint(f'gs://broad-ukbb/broad.freeze_{freeze}/temp/joint_gnomad_ukbb.nfe.ht',
                                                                overwrite=args.overwrite)

        #joint_ht = hl.read_table(f'gs://broad-ukbb/broad.freeze_{freeze}/temp/joint_gnomad_ukbb.ht')
        joint_ht = hl.read_table(f'{variant_qc_prefix(data_source, freeze)}/assessment/joint_gnomad_broad_regeneron.{run_hash}.ht')
        nfe_ht = hl.read_table(f'gs://broad-ukbb/broad.freeze_{freeze}/temp/joint_gnomad_ukbb.nfe.ht')

        joint_ht = joint_ht.filter(hl.is_defined(nfe_ht[joint_ht.key]))
        logger.info(f'Found {joint_ht.count()} variants in NFE samples in both gnomAD and UKBB')

        # # Note: freq[2] contains the NFE variant allele frequencies in gnomAD
        joint_ht = joint_ht.filter((joint_ht.freq[2].AF <= 0.05) & (joint_ht.freq[2].AF > 0.00001), keep=True)
        joint_ht = joint_ht.annotate(gnomad_af_bin=hl.case()
                                    .when((joint_ht.freq[0].AF <= 0.00001), "(0, 0.00001]")
                                    .when((joint_ht.freq[0].AF > 0.00001) & (joint_ht.freq[0].AF <= 0.0001), "(0.00001, 0.0001]")
                                    .when((joint_ht.freq[0].AF > 0.0001) & (joint_ht.freq[0].AF <= 0.001), "(0.0001, 0.001]")
                                    .when((joint_ht.freq[0].AF > 0.001) & (joint_ht.freq[0].AF <= 0.01), "(0.001, 0.01]")
                                    .when((joint_ht.freq[0].AF > 0.01) & (joint_ht.freq[0].AF <= 0.02), "(0.01, 0.02]")
                                    .when((joint_ht.freq[0].AF > 0.02) & (joint_ht.freq[0].AF <= 0.03), "(0.02, 0.03]")
                                    .when((joint_ht.freq[0].AF > 0.03) & (joint_ht.freq[0].AF <= 0.04), "(0.03, 0.04]")
                                    .when((joint_ht.freq[0].AF > 0.04) & (joint_ht.freq[0].AF <= 0.05), "(0.03, 0.05]")
                                    .default(hl.null(hl.tstr)))
        joint_ht = (joint_ht.group_by(joint_ht.gnomad_af_bin, joint_ht.variant_type, joint_ht.pass_gnomad, joint_ht.pass_broad)
                    .aggregate(n_shared=hl.agg.count_where(hl.is_defined(joint_ht.pass_broad) & hl.is_defined(joint_ht.pass_gnomad)),
                            n_gnomad_unique=hl.agg.count_where(
                                    hl.is_defined(joint_ht.pass_gnomad) & ~hl.is_defined(joint_ht.pass_broad)),
                            n=hl.agg.count()))
        joint_ht = joint_ht.annotate(frac_shared=joint_ht.n_shared/joint_ht.n,
                                        frac_gnomad_unique=joint_ht.n_gnomad_unique/joint_ht.n)
        joint_ht = joint_ht.checkpoint(f'gs://broad-ukbb/broad.freeze_{freeze}/temp/joint_gnomad_ukbb.nfe.agg.ht', overwrite=args.overwrite)
        joint_ht.show(50)
        joint_ht.export(f'gs://broad-ukbb/broad.freeze_{freeze}/temp/joint_gnomad_ukbb.nfe.agg.tsv.gz')


    if args.bin_rf_probs_nfe:
        joint_ht = hl.read_table('gs://broad-ukbb/broad.freeze_4/temp/joint_gnomad_ukbb.nfe.ht')
        joint_ht = joint_ht.group_by(tp_bin = joint_ht.rf_probability // 0.01).aggregate(n_pass_gnomad=hl.agg.count_where(hl.is_defined(joint_ht.pass_gnomad) & joint_ht.pass_gnomad),
                                                                                                  n_fail_gnomad=hl.agg.count_where(hl.is_defined(joint_ht.pass_gnomad) & ~joint_ht.pass_gnomad),
                                                                                                  n_unique_ukbb=hl.agg.count_where(~hl.is_defined(joint_ht.pass_gnomad)),
                                                                                                  n=hl.agg.count())
        joint_ht = joint_ht.checkpoint('gs://broad-ukbb/broad.freeze_4/temp/joint_gnomad_ukbb.nfe.probs.ht', overwrite=args.overwrite)
        joint_ht.show(50)
        joint_ht.export('gs://broad-ukbb/broad.freeze_4/temp/joint_gnomad_ukbb.nfe.probs.tsv.gz')


    if args.eval_sib_singletons:
        sib_ht = hl.read_table(f'{variant_qc_prefix(data_source, freeze)}/variant_annotations/sibling_singletons.test.ht')
        joint_ht = hl.read_table(f'{variant_qc_prefix(data_source, freeze)}/assessment/joint_gnomad_broad_regeneron.{run_hash}.ht')
        sib_ht = sib_ht.annotate(**joint_ht[sib_ht.key])

        filter_summary = sib_ht.group_by(sib_ht.allele_type_1, sib_ht.was_split_1, sib_ht.pass_gnomad, sib_ht.pass_broad, sib_ht.pass_regeneron).aggregate(n=hl.agg.count())
        filter_summary.show(50)
        filter_summary.export(f'{variant_qc_prefix(data_source, freeze)}/assessment/sib_singletons_summmary.{run_hash}.tsv.gz')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--freeze', help='Current freeze #', default=CURRENT_FREEZE, type=int)
    parser.add_argument('--run_hash', help='Run hash for RF results to be ranked.')

    parser.add_argument('--join_tables', help='Annotate filtering status and join variant HTs for gnomAD and UKBB exomes.', action='store_true')
    parser.add_argument('--tp_threshold', help='Probability threshold for random forest filtering', type=float)

    parser.add_argument('--bin_rf_probs_overall', help='Aggregate gnomAD pass/fail statistics by RF probabilities for all variants.', action='store_true')
    parser.add_argument('--bin_rf_probs_nfe', help='Aggregate gnomAD pass/fail statistics by RF probabilities for variants found in NFE samples.', action='store_true')

    parser.add_argument('--aggregate_overall', help='Aggregate gnomAD variant recovery statistics for overall callset.', action='store_true')
    parser.add_argument('--aggregate_nfe', help='Aggregate gnomAD variant recovery statistics for non-Finnish Europeans in callset.', action='store_true')
    parser.add_argument('--filter_nfe', help='Identify variants found in non-Finnish Europeans in gnomAD and UKBB callsets.', action='store_true')
    parser.add_argument('--eval_sib_singletons', help='Evaluate performance of RF and Regeneron filtering on sibling singleton variants.', action='store_true')


    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()


    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
