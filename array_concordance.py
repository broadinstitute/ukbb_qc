from gnomad_hail import *
from gnomad_hail.utils.liftover import *
from ukbb_qc.resources import *
import hail as hl
import argparse


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("load_array_data")
logger.setLevel(logging.INFO)


def main(args):
    hl.init(log='/array_concordance.log')

    if args.import_arrays:
        logger.info('Importing array data (autosomes only) into hail Matrix Table...')

        array_files = [(get_array_data_path('bed', f'chr{chrom}'),
                        get_array_data_path('bim', f'chr{chrom}')) for chrom in range(1, 23)]

        mts = [hl.import_plink(bed=f[0],
                               bim=f[1],
                               fam=get_array_data_path('fam','chr22'),
                               reference_genome="GRCh37",
                               a2_reference=False,
                               skip_invalid_loci=True) for f in array_files]

        array_mt = mts[0].union_rows(*mts[1:])
        array_mt.write(array_mt_path(), overwrite=args.overwrite)
        array_variants, array_samples = array_mt.count()
        logger.info(f'{array_variants} variants and {array_samples} samples found in array data')


    if args.liftover_arrays:
        logger.info('Lifting over array data from GRCh37 to GRCh38...')
        array_mt = hl.read_matrix_table(array_mt_path())
        rg37 = hl.get_reference('GRCh37')
        rg38 = hl.get_reference('GRCh38')
        rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)

        array_mt = array_mt.annotate_rows(new_locus=hl.liftover(array_mt.locus, 'GRCh38', include_strand=True),
                                          old_locus=array_mt.locus)

        # Filter out any liftover that has changed strands because they were kept in the liftover file
        array_mt = array_mt.filter_rows(hl.is_defined(array_mt.new_locus) & ~array_mt.new_locus.is_negative_strand)
        array_mt = array_mt.key_rows_by(locus=array_mt.new_locus.result, alleles=array_mt.alleles)

        array_mt = array_mt.checkpoint(array_mt_path(liftover=True), overwrite=args.overwrite)
        array_variants, array_samples = array_mt.count()
        logger.info(f'{array_variants} variants and {array_samples} samples found in liftover array data')


    if args.array_concordance:
        logger.info('Checking concordance between exome and array data...')
        array_mt = hl.read_matrix_table(array_mt_path(liftover=True))
        logger.info('Mapping array sample names to exome sample names...')
        sample_map = hl.import_table(array_sample_map, delimiter=',')
        sample_map = sample_map.key_by(s=sample_map.eid_26041)
        array_mt = array_mt.key_cols_by(s=sample_map[array_mt.s].eid_sample)
        array_mt = array_mt.filter_cols(hl.is_defined(array_mt.s))
        array_variants, array_samples = array_mt.count()
        logger.info(f'{array_samples} samples found in array to exome name map')
        array_mt = hl.variant_qc(array_mt)

        rg38 = hl.get_reference('GRCh38')
        rg38.add_sequence(
            'gs://hail-common/references/Homo_sapiens_assembly38.fasta.gz',
            'gs://hail-common/references/Homo_sapiens_assembly38.fasta.fai'
            )

        # NOTE: This step doesn't currently work
        # see https://github.com/hail-is/hail/issues/5371
        logger.info('Checking SNPs for reference mismatches')
        array_missmatch_ht = annotate_snp_mismatch(array_mt.rows(), rg38)
        mismatch = check_mismatch(array_missmatch_ht)
        logger.info('{} total SNPs'.format(mismatch['total_variants']))
        logger.info('{} SNPs on minus strand'.format(mismatch['negative_strand']))
        logger.info('{} reference mismatches in SNPs'.format(mismatch['total_mismatch']))
        logger.info('{} mismatches on minus strand'.format(mismatch['negative_strand_mismatch']))

        array_mt = array_mt.filter_rows((array_mt.variant_qc.call_rate > args.call_rate_cutoff) &
                                        ~array_missmatch_ht[array_mt.row_key].reference_mismatch)
        exome_mt = get_ukbb_data(args.data_source, args.freeze, adj=True, split=False)
        exome_mt = hl.variant_qc(exome_mt)

        # NOTE: Filter to autosomes because of adjusted sex ploidies in hardcalls mt (hail throws a ploidy 0 error)
        exome_mt = exome_mt.filter_rows(exome_mt.locus.in_autosome())

        # Renaming exome samples to match array data
        exome_mt = exome_mt.key_cols_by(s=exome_mt.s.split("_")[1])
        exome_mt = exome_mt.filter_rows((hl.len(exome_mt.alleles) == 2) &
                                        hl.is_snp(exome_mt.alleles[0], exome_mt.alleles[1]) &
                                        hl.is_defined(array_mt.index_rows(exome_mt.row_key)) &
                                        (exome_mt.variant_qc.call_rate > args.call_rate_cutoff))
        exome_mt = exome_mt.checkpoint(get_mt_checkpoint_path(args.data_source, args.freeze, name="exome_subset_concordance"), overwrite=args.overwrite)
        array_mt = array_mt.filter_rows(hl.is_defined(exome_mt.index_rows(array_mt.row_key)))

        array_mt = array_mt.checkpoint(get_mt_checkpoint_path(args.data_source, args.freeze, name="array_subset_concordance"), overwrite=args.overwrite)

        summary, samples, variants = hl.concordance(array_mt, exome_mt)
        samples = samples.annotate(
            num_gt_con=samples.concordance[2][2] + samples.concordance[3][3] + samples.concordance[4][4],
            num_gt=hl.sum(samples.concordance[2][2:]) + hl.sum(samples.concordance[3][2:]) + hl.sum(samples.concordance[4][2:]))
        samples = samples.annotate(prop_gt_con=samples.num_gt_con / samples.num_gt)
        min_prop_gt_con = samples.aggregate(hl.agg.min(samples.prop_gt_con))
        logger.info(f'Minimum proportion concordance: {min_prop_gt_con}')
        samples.write(array_sample_concordance_path(args.data_source, args.freeze), overwrite=args.overwrite)
        variants.write(array_variant_concordance_path(args.data_source, args.freeze), overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')

    array_import = parser.add_argument_group('Import array plink files and write Matrix Table')
    array_import.add_argument('-i', '--import_arrays', help='Import array data', action='store_true')
    array_import.add_argument('-l', '--liftover_arrays', help='Liftover array Matrix Table to GRCh38', action='store_true')

    concordance = parser.add_argument_group("Compute array concordance with exomes")
    concordance.add_argument('-s', '--data_source', help='Data source', choices=['regeneron', 'broad'], default='broad')
    concordance.add_argument('-f', '--freeze', help='Data freeze to use', default=CURRENT_FREEZE, type=int)
    concordance.add_argument('-c', '--array_concordance', help='Compute array concordance', action='store_true')
    concordance.add_argument('--call_rate_cutoff', help='Call rate cutoff.', type=float,
                        default=0.95)

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
