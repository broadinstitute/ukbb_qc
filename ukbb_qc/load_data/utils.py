import hail as hl
import logging
from ukbb_qc.resources.basics import array_sample_map
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("load_data")
logger.setLevel(logging.INFO)


# Sample resources
def import_array_exome_id_map_ht(freeze: int = CURRENT_FREEZE, array_sample_map: str): -> hl.Table
    """
    Imports file linking array IDs to exome IDs into Table

    :param int freeze: One of the data freezes. Default is CURRENT_FREEZE.
    :param str array_sample_map: Path to linking csv for data freeze.
    :return: Table with array IDs mapped to exome IDs
    :rtype: hl.Table
    """
    sample_map_ht = hl.import_table(array_sample_map(freeze), delimiter=",", quote='"')
    sample_map_ht = sample_map_ht.key_by(s=sample_map_ht.eid_sample)
    logger.info(
            f"Total number of IDs in the array to exome sample map: {sample_map_ht.count()}..."
        )
    return sample_map_ht


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
