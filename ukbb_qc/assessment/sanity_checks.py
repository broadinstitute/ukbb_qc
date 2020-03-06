from gnomad.utils.gnomad_functions import *
from ukbb_qc.resources import *


def summarize_mt(mt: hl.MatrixTable) -> hl.Struct:
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

    return var_summary


def check_adj(mt: hl.MatrixTable, mt_adj: hl.MatrixTable) -> bool:
    """
    Checks if MatrixTable has been filtered using adj criteria by checking allele counts pre and post adj filtration
    :param MatrixTable mt: MatrixTable to be checked
    :param MatrixTable mt_adj: MatrixTable filtered using adj criteria
    :return: Bool of whether MatrixTable has been adj filtered
    :rtype: bool
    """

    pre = mt.aggregate_entries(hl.agg.counter(mt.GT.n_alt_alleles()))
    logger.info('\nAllele distribution pre adj filtration: {}'.format(pre))
    post = mt_adj.aggregate_entries(hl.agg.counter(mt_adj.GT.n_alt_alleles()))
    logger.info('\nAllele distribution post adj filtration: {}'.format(post))

    adj = False
    if sum(pre.values()) != sum(post.values()):
        adj = True

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
   
    sample_mismatch = False
    mt_samples = mt.aggregate_cols(hl.agg.collect_as_set(mt.s))
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
