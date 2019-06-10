from ukbb_qc.resources import *


def remove_hard_filter_samples(data_source: str, freeze: int, mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Removes samples that failed hard filters from mt

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param MatrixTable mt: Input MatrixTable
    :return: MatrixTable with samples removed
    :rtype: MatrixTable
    """

    # read in hard filters ht
    ht = hl.read_table(hard_filters_ht_path(data_source, freeze))

    # remove unnecessary fields from hard filters ht
    ht = ht.select('hard_filters')

    # remove hard filtered samples
    mt = mt.filter_cols(hl.len(ht[mt.col_key].hard_filters) == 0)
    return mt
