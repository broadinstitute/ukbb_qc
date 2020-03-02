from gnomad_hail.resources.resource_utils import DataException
import gnomad_hail.resources.grch38 as grch38
from .resource_utils import CURRENT_FREEZE, CURRENT_HAIL_VERSION, DATA_SOURCES, FREEZES
from typing import Dict


def get_truth_sample_info(
    data_source: str, freeze: int = CURRENT_FREEZE
) -> Dict[str, str]:
    truth_samples = {
        "syndip": {
            "s": "CHMI_CHMI3_Nex1",
            "truth_mt": grch38.syndip,
            "mt": truth_sample_mt_path(data_source, freeze, "syndip"),
            "bed": grch38.syndip_hc_intervals,
        },
        "na12878": {
            "s": "Coriell_NA12878_NA12878",
            "truth_mt": grch38.na12878_giab,
            "mt": truth_sample_mt_path(data_source, freeze, "na12878"),
            "bed": grch38.na12878_giab_hc_intervals,
        },
    }
    return truth_samples


def variant_qc_prefix(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    if data_source not in DATA_SOURCES:
        raise DataException("This data_source is currently not present")
    if freeze not in FREEZES:
        raise DataException("This freeze is currently not present")

    return f"gs://broad-ukbb/{data_source}.freeze_{freeze}/variant_qc"


def var_annotations_ht_path(data_source: str, freeze: int, annotation_type: str) -> str:
    """
    Get variant-level annotations

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param str annotation_type: Tag describing variant-level annotations in HT, e.g., "vep", "qc_stats", "call_stats", "vqsr"
    :return: Path to annotations Table
    :rtype: str
    """
    return f"{variant_qc_prefix(data_source, freeze)}/variant_annotations/{annotation_type}.ht"


def info_ht_path(data_source: str, freeze: int, split: bool = True) -> str:
    """
    Get variant-level annotations

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param bool split: Whether the dataset should be split (only applies to raw=False)
    :return: Path to annotations Table
    :rtype: str
    """
    return f'{variant_qc_prefix(data_source, freeze)}/variant_annotations/info{"_split" if split else ""}.ht'


def truth_sample_mt_path(data_source: str, freeze: int, truth_sample: str) -> str:
    """
    Get truth sample path

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param str truth_sample: Name of the truth sample
    :return: Path to annotations Table
    :rtype: str
    """
    return f"{variant_qc_prefix(data_source, freeze)}/truth_samples/{truth_sample}.ht"


def sample_annotations_table_path(
    data_source: str, freeze: int, annotation_type: str
) -> str:
    """
    Get sample-level annotations

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param str annotation_type: Tag describing variant-level annotations in HT, e.g., "family_stats"
    :return: Path to annotations HT
    :rtype: str
    """
    return f"{variant_qc_prefix(data_source, freeze)}/sample_annotations/{annotation_type}.ht"


def rf_run_hash_path(data_source: str, freeze: int = CURRENT_FREEZE):
    """
    Returns the path to the json file containing the RF runs list.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Path to json file
    :rtype: str
    """

    return f"{variant_qc_prefix(data_source, freeze)}/rf/runs.json"


def rf_annotated_path(
    data_source: str, freeze: int = CURRENT_FREEZE, adj: bool = False
) -> str:
    """
    Returns the path to the RF-ready annotated HT

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param bool adj: Whether to load 'adj' or 'raw'
    :return:
    """

    return f'{variant_qc_prefix(data_source, freeze)}/rf/rf_annotated.{"adj" if adj else "raw"}.ht'


def rf_path(
    data_source: str,
    freeze: int = CURRENT_FREEZE,
    data: str = "rf_result",
    run_hash: str = None,
) -> str:
    """

    Gets the path to the desired RF data.
    Data can take the following values:
        - 'training': path to the training data for a given run
        - 'model': path to pyspark pipeline RF model
        - 'rf_result' (default): path to HT containing result of RF filtering

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param str data: One of 'pre_rf', 'training', 'model' or 'rf_result' (default)
    :param str run_hash: Hash of RF run to load
    :return:
    """

    extension = "model" if data == "model" else "ht"
    return f"{variant_qc_prefix(data_source, freeze)}/rf/{run_hash}/{data}.{extension}"


def score_ranking_path(
    data_source: str, freeze: int, data: str, binned: bool = False,
) -> str:
    """
    Returns the path to non-RF metrics score rankings Tables, e.g.:
    * vqsr

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param data: The score data to return
    :param binned: Whether to get the binned data
    :return: Path to desired hail Table
    :rtype: str
    """

    return f'{variant_qc_prefix(data_source, freeze)}/score_rankings/{data}{"_binned" if binned else ""}.ht'


def binned_concordance_path(
    data_source: str, freeze: int, truth_sample: str, metric: str
):
    """
    :param str data_source: 'broad' or 'regeneron'
    :param int freeze: UKBB tranche version
    :param str truth_sample: Which truth sample concordance to analyze (e.g., "NA12878" or "syndip")
    :param str metric: One of the evaluation metrics (or a RF hash)
    :return: Path to Hail Table
    :rtype: str
    """

    return f"{variant_qc_prefix(data_source, freeze)}/rf/{truth_sample}.{metric}.binned_concordance.ht"


def omni_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return "gs://gnomad-public/truth-sets/hail-{0}/1000G_omni2.5.hg38.mt".format(
        hail_version
    )


def mills_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return "gs://gnomad-public/truth-sets/hail-{0}/Mills_and_1000G_gold_standard.indels.hg38.mt".format(
        hail_version
    )


def hapmap_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return "gs://gnomad-public/truth-sets/hail-{0}/hapmap_3.3.hg38.mt".format(
        hail_version
    )


def hapmap_ht_path():
    return "gs://broad-ukbb/resources//hapmap_3.3.hg38.ht"


def dbsnp_ht_path():
    return "gs://gnomad-public/resources/grch38/dbsnp_b151_grch38_all_20180418.ht"


def kgp_high_conf_snvs_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return "gs://gnomad-public/truth-sets/hail-{0}/1000G_phase1.snps.high_confidence.hg38.mt".format(
        hail_version
    )
