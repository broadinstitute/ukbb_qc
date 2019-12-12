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


def import_capture_regions(interval_path: str, output_path: str, ow: bool) -> None:
    '''
    Imports capture region text file into Table and writes Table at specified path

    :param str interval_path: Path to input file
    :param str output_path: Path to output file
    :param bool ow: Whether to overwrite data
    :return: None
    :rtype: None
    '''

    logger.info('Importing capture table')
    capture_ht = hl.import_table(
                    interval_path,
                    no_header=True,
                    impute=True,
                    min_partitions=10)

    # seqnames	start	end	width	strand	target_type	region_type	target_id
    #chr1	11719	12377	659	*	processed_transcript|transcribed_unprocessed_pseudogene|ice_target_1	gnomad	target_1
    capture_ht = capture_ht.transmute(interval=hl.parse_locus_interval(
                                    hl.format('[%s:%s-%s]',
                                               capture_ht.f0,
                                               capture_ht.f1,
                                               capture_ht.f2),
                                    reference_genome='GRCh38'))
    capture_ht = capture_ht.select('interval').key_by('interval')

    capture_ht.describe()
    logger.info('Writing capture ht')
    capture_ht.write(output_path, overwrite=ow)
