from ukbb_qc.resources import *


def remove_hard_filter_samples(data_source: str, freeze: int, t: Union[hl.MatrixTable, hl.Table]) -> Union[hl.MatrixTable, hl.Table]:
    """
    Removes samples that failed hard filters from MatrixTable or Table. 
    Assumes MatrixTable's col key/Table's row key is sample.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param MatrixTable/able t: Input MatrixTable or Table
    :return: MatrixTable or Table with samples removed
    :rtype: MatrixTable or Table
    """

    # read in hard filters ht
    ht = hl.read_table(hard_filters_ht_path(data_source, freeze))

    # remove unnecessary fields from hard filters ht
    ht = ht.select('hard_filters')

    # remove hard filtered samples
    if isinstance(t, hl.MatrixTable):
        t = t.filter_cols(hl.len(ht[t.col_key].hard_filters) == 0) 
    else:
        t.filter((hl.len(ht[t.row_key].hard_filters) == 0)
    return t
