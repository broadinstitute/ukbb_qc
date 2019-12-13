import logging
from ukbb_qc.resources import *


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("utils")
logger.setLevel(logging.INFO)


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

    # get number of hard filtered samples
    hf_samples = ht.aggregate(hl.agg.count_where(hl.len(ht.hard_filters) > 0))
    logger.info(f'Removing {hf_samples} samples that failed hard filters')

    # remove hard filtered samples
    ht = ht.filter(hl.len(ht.hard_filters) == 0)
    if isinstance(t, hl.MatrixTable):
        t = t.filter_cols(hl.is_defined(ht[t.col_key]))
        t = t.annotate_rows(non_refs=hl.agg.count_where(t.GT.is_non_ref()))
        t = t.filter_rows(t.non_refs > 0).drop('non_refs')
    else:
        t.filter((hl.len(ht[t.key].hard_filters) == 0))
    return t


def interval_qc_filter(data_source: str,
                       freeze: int,
                       t: Union[hl.MatrixTable, hl.Table],
                       pct_samples_20x: float = 0.85) -> Union[hl.MatrixTable, hl.Table]:
    """
    Removes poorly covered intervals from a MatrixTable or Table.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param MatrixTable/able t: Input MatrixTable or Table
    :param float pct_samples_20x: Percent of samples with coverage greater than 20X over the interval for filtering
    :return: MatrixTable or Table with samples removed
    :rtype: MatrixTable or Table
    """
    interval_qc_ht = hl.read_table(interval_qc_path(data_source, freeze))
    good_intervals_ht = interval_qc_ht.filter(interval_qc_ht.pct_samples_20x > pct_samples_20x).key_by('interval')

    if isinstance(t, hl.MatrixTable):
        t = t.filter_rows(hl.is_defined(good_intervals_ht[t.locus]))
    else:
        t = t.filter(hl.is_defined(good_intervals_ht[t.locus]))
    return t


def annotate_relationship(relatedness_ht,
                          first_degree_threshold=[0.1767767, 0.4],
                          second_degree_threshold=0.08838835,
                          ibd2_parent_offspring_threshold=0.14):
    relatedness_ht = relatedness_ht.annotate(relationship_classification=hl.case()
                                             .when((relatedness_ht.kin > second_degree_threshold) &
                                                   (relatedness_ht.kin < first_degree_threshold[0]), 'Second-degree')
                                             .when((relatedness_ht.kin > first_degree_threshold[0]) &
                                                   (relatedness_ht.kin < first_degree_threshold[1]) &
                                                   (relatedness_ht.ibd2 >= ibd2_parent_offspring_threshold),
                                                   'Full-sibling')
                                             .when((relatedness_ht.kin > first_degree_threshold[0]) &
                                                   (relatedness_ht.kin < first_degree_threshold[1]) &
                                                   (relatedness_ht.ibd2 < ibd2_parent_offspring_threshold),
                                                   'Parent-child')
                                             .default("None"))

    return relatedness_ht
