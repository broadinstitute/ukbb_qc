import hail as hl
import logging
from typing import Optional
from gnomad.utils.generic import file_exists
from ukbb_qc.resources.basics import (
    array_sample_map_path,
    capture_ht_path,
    excluded_samples_path,
    phenotype_ht_path,
    ukbb_phenotype_path,
)
from ukbb_qc.resources.variant_qc import var_annotations_ht_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("load_data")
logger.setLevel(logging.INFO)


# Sample resources
def import_array_exome_id_map_ht(freeze: int = CURRENT_FREEZE) -> hl.Table:
    """
    Imports file linking array IDs to exome IDs into Table.
    Writes Table with array IDs mapped to exome IDs.

    :param int freeze: One of the data freezes. Default is CURRENT_FREEZE.
    :return: Table with array IDs mapped to exome IDs
    :rtype: hl.Table
    """
    data_source = "broad"

    sample_map_ht = hl.import_table(
        array_sample_map_path(freeze), delimiter=",", quote='"'
    )
    sample_map_ht = sample_map_ht.key_by("eid_sample")
    sample_map_ht = sample_map_ht.transmute(
        batch_num=sample_map_ht.batch,
        batch=sample_map_ht["batch.c"],
        ukbb_app_26041_id=sample_map_ht.eid_26041,
    )
    logger.info(
        f"Total number of IDs in the array to exome sample map: {sample_map_ht.count()}..."
    )

    if file_exists(excluded_samples_path()):
        excluded_samples_ht = hl.import_table(excluded_samples_path(), no_header=True)
        excluded_samples = hl.literal(ht.aggregate(hl.agg.collect_as_set(ht.f0)))
        logger.info(
            f"Total number of samples to exclude: {hl.eval(hl.len(excluded_samples))}"
        )
        sample_map_ht = sample_map_ht.annotate(
            withdrawn_consent=excluded_samples.contains(sample_map_ht.ukbb_app_26041_id)
        )
        withdrawn_ids = sample_map_ht.aggregate(
            hl.agg.count_where(sample_map_ht.withdrawn_consent)
        )
        logger.info(
            f"Total number of IDs with withdrawn consents in sample map ht: {withdrawn_ids}"
        )
    return sample_map_ht


def import_phenotype_ht() -> None:
    """
    Imports UKBB phenotype file as a Table.

    :return: None
    :rtype: None
    """
    phenotype_ht = hl.import_table(ukbb_phenotype_path, impute=True)
    phenotype_ht = phenotype_ht.key_by(s_old=hl.str(phenotype_ht["f.eid"]))
    phenotype_ht.write(phenotype_ht_path(), overwrite=True)


# Interval resources
def import_capture_intervals(interval_path: str, header: bool, overwrite: bool) -> None:
    """
    Imports capture intervals text file into Table and writes Table at specified path

    :param str interval_path: Path to input file
    :param str output_path: Path to output file
    :param bool header: Whether capture intervals file has a header line
    :param bool overwrite: Whether to overwrite data
    :return: None
    :rtype: None
    """
    data_source = "broad"

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
    capture_ht.write(capture_ht_path(data_source), overwrite=overwrite)


def import_vqsr(
    freeze: int,
    vqsr_path: str,
    vqsr_type: str = "AS",
    num_partitions: int = 5000,
    overwrite: bool = False,
    import_header_path: Optional[str] = None,
) -> None:
    """
    Imports vqsr site vcf into a HT

    :param freeze: One of the data freezes. Default is CURRENT_FREEZE.
    :param vqsr_path: Path to input vqsr site vcf. This can be specified as Hadoop glob patterns
    :param vqsr_type: One of `AS` (allele specific) or `AS_TS` (allele specific with transmitted singletons)
    :param num_partitions: Number of partitions to use for the VQSR HT
    :param overwrite: Whether to overwrite imported VQSR HT
    :param import_header_path: Optional path to a header file to use for import
    :return: None
    """
    data_source = "broad"

    logger.info(f"Importing VQSR annotations for {vqsr_type} VQSR...")
    mt = hl.import_vcf(
        vqsr_path,
        force_bgz=True,
        reference_genome="GRCh38",
        header_file=import_header_path,
    ).naive_coalesce(num_partitions)

    ht = mt.rows()
    row_count1 = ht.count()
    ht = hl.split_multi_hts(ht).checkpoint(
        var_annotations_ht_path(
            data_source, freeze, "vqsr" if vqsr_type == "AS" else "AS_TS_vqsr"
        ),
        overwrite=overwrite,
    )
    row_count2 = ht.count()
    logger.info(
        f"Found {row_count1} unsplit and {row_count2} split variants with VQSR annotations"
    )
