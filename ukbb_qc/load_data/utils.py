import logging
from typing import Optional

import hail as hl

from gnomad.resources.grch38.reference_data import clinvar
from gnomad.utils.file_utils import file_exists
from gnomad.utils.sparse_mt import split_info_annotation
from ukbb_qc.resources.basics import (
    array_sample_map_path,
    capture_ht_path,
    excluded_samples_path,
    phenotype_ht_path,
    ukbb_phenotype_path,
)
from ukbb_qc.resources.sample_qc import (
    get_ukbb_array_pcs_path,
    get_ukbb_array_pcs_ht_path,
)
from ukbb_qc.resources.variant_qc import clinvar_pathogenic_ht_path, var_annotations_ht_path
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
        excluded_samples = hl.literal(
            excluded_samples_ht.aggregate(hl.agg.collect_as_set(excluded_samples_ht.f0))
        )
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


def load_ukbb_array_pcs() -> hl.Table:
    """
    Imports UKBB genotype array sample QC information into a Table and extracts the ancestry PCs.

    :return: array ancestry PCs Table
    :rtype: Table
    """
    array_pcs = hl.import_table(
        get_ukbb_array_pcs_path(), impute=True, delimiter="\t", no_header=True,
    )
    array_pcs = array_pcs.select(
        s=hl.str(array_pcs.f1),
        scores=hl.array([array_pcs[f"f{i+26}"] for i in range(0, 40)]),
    )
    array_pcs = array_pcs.key_by("s")
    array_pcs.write(get_ukbb_array_pcs_ht_path())


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


# Variant resources
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
    ).repartition(num_partitions)

    ht = mt.rows()

    ht = ht.annotate(
        info=ht.info.annotate(
            AS_VQSLOD=ht.info.AS_VQSLOD.map(lambda x: hl.float(x)),
            AS_QUALapprox=ht.info.AS_QUALapprox.split("\|")[1:].map(
                lambda x: hl.int(x)
            ),
            AS_VarDP=ht.info.AS_VarDP.split("\|")[1:].map(lambda x: hl.int(x)),
            AS_SB_TABLE=ht.info.AS_SB_TABLE.split("\|").map(
                lambda x: x.split(",").map(lambda y: hl.int(y))
            ),
        )
    )

    ht = ht.checkpoint(
        var_annotations_ht_path(
            f'{"vqsr" if vqsr_type == "AS" else "AS_TS_vqsr"}.unsplit',
            data_source,
            freeze,
        ),
        overwrite=overwrite,
    )

    unsplit_count = ht.count()
    ht = hl.split_multi_hts(ht)

    ht = ht.annotate(
        info=ht.info.annotate(**split_info_annotation(ht.info, ht.a_index)),
    )

    ht = ht.checkpoint(
        var_annotations_ht_path(
            "vqsr" if vqsr_type == "AS" else "AS_TS_vqsr", data_source, freeze,
        ),
        overwrite=overwrite,
    )
    split_count = ht.count()
    logger.info(
        f"Found {unsplit_count} unsplit and {split_count} split variants with VQSR annotations"
    )


def load_clinvar_path() -> hl.Table:
    """
    Filters the most recent GRCh38 ClinVar dataset in gnomad_methods resources to only pathogenic variants.

    Removes variants with no star assertions and filtering to only (likely) pathogenic variants and no
    conflicting clinical interpretations
reference_data.clinvar
    :return: Table of ClinVar variants filtered to only pathogenic
    """
    logger.info(
        f"Filtering ClinVar version {} to only pathogenic variants"
    )

    clinvar_ht = clinvar.ht()
    logger.info(f"Found {clinvar_ht.count()} variants before filtering")

    no_star_assertions = hl.literal(
        {
            "no_assertion_provided",
            "no_assertion_criteria_provided",
            "no_interpretation_for_the_single_variant",
        }
    )clinvar.default_version
    clinvar_ht = clinvar_ht.filter(
        hl.set(clinvar_ht.info.CLNREVSTAT).intersection(no_star_assertions).length()
        > 0,
        keep=False,
    )
    logger.info(
        f"Found {clinvar_ht.count()} variants after removing variants without assertions"
    )
    clinvar_ht = clinvar_ht.filter(
        clinvar_ht.info.CLNSIG.map(lambda x: x.lower())
        .map(lambda x: x.contains("pathogenic"))
        .any(lambda x: x),
        keep=True,
    )
    logger.info(
        f"Found {clinvar_ht.count()} variants after filtering to (likely) pathogenic variants"
    )
    clinvar_ht = clinvar_ht.filter(
        hl.is_defined(clinvar_ht.info.CLNSIGCONF), keep=False
    )
    clinvar_ht.write(clinvar_pathogenic_ht_path(clinvar.default_version), overwrite=True)
    logger.info(
        f"Found {clinvar_ht.count()} variants after filtering to variants without CLNSIGCONF (conflicting clinical \
        interpretations)"
    )


