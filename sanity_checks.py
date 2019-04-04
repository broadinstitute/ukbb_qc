import argparse
import logging
from resources import *


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sanity_checks")
logger.setLevel(logging.INFO)


def summary(mt: hl.MatrixTable) -> hl.Struct:
    """
    Gets a summary of variants in mt, prints number of variants to stdout, and checks that each chromosome has variant calls
    :param MatrixTable mt: Raw MatrixTable to be checked
    :return: Struct of MatrixTable variant summary
    :rtype: Struct
    """

    var_summary = hl.summarize_variants(mt, show=False)
    logger.info('Dataset has {} variants'.format(var_summary.n_variants))
   
    # check that all contigs have variant calls 
    for contig in var_summary.contigs:
        if var_summary.contigs[contig] == 0:
            logger.info('{} has no variants called'.format(contig))

    return hl.summarize_variants(mt, show=False)


def alt_alleles(mt: hl.MatrixTable) -> hl.expr.StructExpression:
    """
    Gets allele count distributions with n_alt_alleles() function
    :param MatrixTable mt: MatrixTable to be checked
    :return: StructExpression of allele count distribution
    :rtype: StructExpression
    """

    return mt.aggregate_entries(hl.agg.counter(mt.GT.n_alt_alleles()))
    
 
def adj_check(mt: hl.MatrixTable, mt_adj: hl.MatrixTable) -> bool:
    """
    Checks if MatrixTable has been filtered using adj criteria by checking allele counts pre and post adj filtration
    :param MatrixTable mt: MatrixTable to be checked
    :param MatrixTable mt_adj: MatrixTable filtered using adj criteria
    :return: Bool of whether MatrixTable has been adj filtered
    :rtype: bool
    """

    pre = alt_alleles(mt)
    logger.info('\nAllele distribution pre adj filtration: {}'.format(pre))
    post = alt_alleles(mt_adj)
    logger.info('\nAllele distribution post adj filtration: {}'.format(post))

    adj = True
    for field in pre:
        if field in post:
            if post[field] < pre[field]:
                logger.info('Allele {} had fewer allele counts post adj filtration'.format(field))
                adj = False
            elif post[field] > pre[field]:
                logger.info('Allele {} had more alleles post adj filtration'.format(field))
                adj = False
            else:
                logger.info('Allele {} has the same allele count pre/post adj filtration'.format(field))
        else:
            logger.info('Allele {} was present pre-adj filtration but not post'.format(field))

    return adj


def write_samples(sample_qc_path: str, sample_set: set, name: str) -> None:
    """
    Writes sample mismatches to output txt
    :param str sample_qc_path: Path to output sample_qc bucket
    :param set sample_set: Set of sample mismatches
    :param str name: Name of output txt
    :return: None
    :rtype: None
    """

    out = f'{sample_qc_path}/{name}.txt'
    with hl.hadoop_open(out, 'w') as o:
        o.write('\n'.join(sample_set) + '\n')


def sample_check(mt: hl.MatrixTable, ht: hl.Table, sample_qc_path: str) -> bool:
    """
    Checks for sample mismatch between MatrixTable and provided list of samples; if mismatch, writes to output txt
    :param MatrixTable mt: MatrixTable containing samples to be checked
    :param Table ht: Table with one column containing expected samples
    :param str sample_qc_path: Path to output sample_qc bucket (for output txts)
    :return: Bool with whether there was a sample mismatch
    :rtype: bool
    """
   
    # create bool to return 
    sample_mismatch = False

    # get samples from mt
    mt_samples = mt.aggregate_cols(hl.agg.collect_as_set(mt.s))
   
    # get samples from samptable
    ht_samples = ht.aggregate(hl.agg.collect_as_set(ht.f0)) 

    # check for samples in mt that aren't in sample list
    mt_extras = mt_samples.difference(ht_samples)
    if len(mt_extras) > 0:
        sample_mismatch = True
        write_samples(sample_qc_path, mt_extras, 'mt_unique_samples')
        
    # check for samples that are in sample list that aren't in mt
    ht_extras = ht_samples.difference(mt_samples)
    if len(ht_extras) > 0:
        sample_mismatch = True
        write_samples(sample_qc_path, ht_extras, 'sample_list_unique_samples')

    return sample_mismatch
    

def hail_qc(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Annotates MatrixTable with variant and sample QC
    :param MatrixTable mt: MatrixTable to be annotated
    :return: MatrixTable with variant/sample QC annotations
    :rtype: MatrixTable
    """

    mt = hl.variant_qc(mt)
    mt = hl.sample_qc(mt)
    return mt


def main(args):

    hl.init()

    datasource = args.datasource
    if args.freeze:
        freeze = args.freeze
    else:
        freeze = CURRENT_FREEZE

    # get output path for sample qc files
    sample_qc_path = sample_qc_prefix(datasource, freeze)
    
    logger.info('Loading raw matrixtable of dataset and reading in list of provided samples')
    mt = hl.read_matrix_table(raw_mt_path(datasource, freeze))
    ht = hl.import_table(sample_list_path(datasource, freeze), no_header = True)

    logger.info('Getting variant summary for mt')
    var_summary = summary(mt)
    logger.info('\nVariant summary: {}'.format(var_summary))

    logger.info('Annotating mt with sample and variant qc')
    mt = hail_qc(mt)
   
    logger.info('Annotating with adj')
    mt = annotate_adj(mt)

    # NOTE writing out mt here as no other annotations are added to mt after this point
    logger.info('Writing out adj, sample/variant qc annotated mt')
    mt = mt.checkpoint(adj_mt_path(datasource, freeze), overwrite = args.overwrite)   
 
    logger.info('Checking for adj filtration')
    mt_adj = filter_to_adj(mt)
    adj = adj_check(mt, mt_adj)
   
    if adj:
        logger.info('Dataset was adj filtered')
    else:
        logger.info('Dataset was not adj filtered')

    logger.info('Checking if all samples provided are in dataset')
    if sample_check(mt, ht, sample_qc_path):
        logger.info('\nSample mismatch between dataset and provided list of samples')
    else:
        logger.info('\nAll samples in dataset match provided list of samples')
    
    logger.info('Sanity checks complete')
    

if __name__ == '__main__':

    # Create argument parser
    parser = argparse.ArgumentParser(description='This script does basic sanity checks on a MatrixTable')
    parser.add_argument('-d', '--datasource', help='Data source', choices=['regeneron', 'broad'], default='broad')
    parser.add_argument('-f', '--freeze', help='Current freeze #', type=int)
    parser.add_argument('-o', '--overwrite', help='Overwrite pre-existing data', action='store_true', default=True)
    args = parser.parse_args()
    main(args)
