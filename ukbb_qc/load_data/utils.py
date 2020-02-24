import hail as hl
import logging
from gnomad_hail.utils.generic import file_exists
from ukbb_qc.resources.basics import (
                            array_sample_map_path,
                            excluded_samples_path, ukbb_phenotype_path
                            )
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("load_data")
logger.setLevel(logging.INFO)


# Sample resources
def import_array_exome_id_map_ht(freeze: int = CURRENT_FREEZE) -> hl.Table:
    """
    Imports file linking array IDs to exome IDs into Table

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes. Default is CURRENT_FREEZE.
    :return: Table with array IDs mapped to exome IDs
    :rtype: hl.Table
    """
    sample_map_ht = hl.import_table(array_sample_map_path(freeze), delimiter=",", quote='"')
    sample_map_ht = sample_map_ht.key_by(s=sample_map_ht.eid_sample)
    sample_map_ht = sample_map_ht.transmute(
        batch_num=sample_map_ht.batch,
        batch=sample_map_ht['batch.c']
        )
    logger.info(
            f"Total number of IDs in the array to exome sample map: {sample_map_ht.count()}..."
        )

    if file_exists(excluded_samples_path(freeze)):
        excluded_samples_ht = hl.import_table(excluded_samples_path(freeze), no_header=True)
        excluded_samples = hl.literal(ht.aggregate(hl.agg.collect_as_set(ht.f0)))
        logger.info(
            f"Total number of samples to exclude: {hl.eval(hl.len(excluded_samples))}"
            )
        sample_map_ht = sample_map_ht.annotate(withdrawn_consent=excluded_samples.contains(sample_map_ht.eid_26041))
        withdrawn_ids = sample_map_ht.aggregate(hl.agg.count_where(sample_map_ht.withdrawn_consent))
        logger.info(
            f"Total number of IDs with withdrawn consents in sample map ht: {withdrawn_ids}"
            )
    return sample_map_ht


def import_phenotype_ht() -> hl.Table:
    """
    Imports UKBB phenotype file as a Table.

    :return: None
    :rtype: None
    """
    phenotype_ht = hl.import_table(ukbb_phenotype_path, impute=True)
    phenotype_ht = phenotype_ht.key_by(s_old=hl.str(phenotype_ht['f.eid']))
    phenotype_ht.write(phenotype_ht_path(), overwrite=True)


# Interval resources
def import_capture_intervals(interval_path: str, output_path: str, header: bool, overwrite: bool) -> None:
    """
    Imports capture intervals text file into Table and writes Table at specified path

    :param str interval_path: Path to input file
    :param str output_path: Path to output file
    :param bool header: Whether capture intervals file has a header line
    :param bool overwrite: Whether to overwrite data
    :return: None
    :rtype: None
    """
    logger.info("Importing capture table")
    capture_ht = hl.import_table(
        interval_path, no_header=not header, impute=True, min_partitions=10
    )

    # seqnames	start	end	width	strand	target_type	region_type	target_id
    capture_ht = capture_ht.transmute(
        interval=hl.parse_locus_interval(
            hl.format("[%s:%s-%s]", capture_ht.f0, capture_ht.f1, capture_ht.f2),
            reference_genome="GRCh38",
        )
    )
    capture_ht = capture_ht.select("interval").key_by("interval")

    capture_ht.describe()
    logger.info("Writing capture ht")
    capture_ht.write(output_path, overwrite=overwrite)
