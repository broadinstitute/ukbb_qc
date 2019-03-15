from gnomad_hail import *

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("check_regeneron_vcf")
logger.setLevel(logging.INFO)

def main(args):
    hl.init(log='/check_regeneron_vcf.log')

    rg = hl.get_reference('GRCh37')
    grch37_contigs = [x for x in rg.contigs if not x.startswith('GL') and not x.startswith('M')]

    contig_dict = dict(zip(grch37_contigs, ['chr'+x for x in grch37_contigs]))
    print(contig_dict)

    nf_mt = hl.import_vcf('gs://gnomad-ukbb/data/vcfs/nf/NF_by_chrom/UKB_Freeze_Four.*.NF.vcf.gz',
                       force_bgz=True, array_elements_required=False, reference_genome='GRCh38', contig_recoding=contig_dict)
    # NOTE: These VCFs have been deleted
    nf_mt.write('gs://broad-ukbb/data/regeneron.freeze_4.nf.mt', overwrite=args.overwrite)
    nf_mt = hl.read_matrix_table('gs://broad-ukbb/data/regeneron.freeze_4.nf.mt')
    nf_mt.describe()
    nf_variants, nf_samples = nf_mt.count()
    logger.info(f'{nf_variants} variants and {nf_samples} samples found in NF VCF')

    gl_mt = hl.import_vcf('gs://broad-ukbb/data/vcfs/gl/GL_by_chrom/UKB_Freeze_Four.*.GL.vcf.gz',
                       force_bgz=True, array_elements_required=False, reference_genome='GRCh38', contig_recoding=contig_dict)
    gl_mt.write('gs://broad-ukbb/data/regeneron.freeze_4.gl.mt', overwrite=args.overwrite)
    gl_mt = hl.read_matrix_table('gs://broad-ukbb/data/regeneron.freeze_4.gl.mt')
    gl_mt.describe()
    gl_variants, gl_samples = gl_mt.count()
    logger.info(f'{gl_variants} variants and {gl_samples} samples found in GL VCF')


    # Are all the variants in GL found in NF?
    nf_ht = nf_mt.rows()
    gl_ht = gl_mt.rows()
    joint_ht = nf_ht.join(gl_ht)
    logger.info(f'{joint_ht.count()} variants found in joint VCF')

# TODO: update filepaths after Julia's updated resources file gets merged

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite pre-existing data', action='store_true')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
