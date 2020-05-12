import hail as hl
from typing import Any, Dict, Optional, Union

from gnomad.resources.resource_utils import DataException
import gnomad.resources.grch38 as grch38
from .resource_utils import CURRENT_FREEZE, CURRENT_HAIL_VERSION, DATA_SOURCES, FREEZES


SYNDIP = "CHMI_CHMI3_Nex1"
"""
String representation for syndip truth sample
"""

NA12878 = "Coriell_NA12878_NA12878"
"""
String representation for NA12878 truth sample
"""

TRUTH_SAMPLES = [SYNDIP, NA12878]
"""
List containing string representation for truth samples (syndip, NA12878)
"""


def clinvar_pathogenic_ht_path(version: str) -> str:
    """
    Returns path to a ClinVar Table that has been filtered to only pathogenic variants

    :param str version: ClinVar version
    :return: Path to ClinVar pathogenic variant Table
    """
    return f"gs://broad-ukbb/resources/clinvar_{version}.pathogenic.ht"


def get_truth_sample_data(
    data_source: str,
    freeze: int = CURRENT_FREEZE,
    truth_sample: str = None,
    data_type: str = None,
) -> Union[str, hl.Table, hl.MatrixTable]:
    """
    Returns relevant data for truth samples. 

    Current truth samples available are syndip and na12878.

    The following information is available for each truth sample:
    - s: sample name in the callset
    - truth_mt: truth sample MatrixTable
    - hc_intervals: high confidence interval Table in truth sample
    - callset_truth_mt: truth sample MatrixTable (subset from callset)

    :param str data_source: One of 'regeneron' or 'broad'
    :param str freeze: One of the data freezes
    :param str truth_sample: Name of the truth sample. One of 'syndip' or 'na12878'
    :param str data_type: Truth sample data type. One of 's', 'truth_mt', 'hc_intervals', or 'callset_truth_mt'. 
    :return: Sample name, Table, or MatrixTable of requested truth sample data
    :rtype: Union[str, hl.Table, hl.MatrixTable]
    """
    truth_samples = {
        "syndip": {
            "s": SYNDIP,
            "truth_mt": grch38.syndip.mt(),
            "hc_intervals": grch38.syndip_hc_intervals.ht(),
        },
        "na12878": {
            "s": NA12878,
            "truth_mt": grch38.na12878_giab.mt(),
            "hc_intervals": grch38.na12878_giab_hc_intervals.ht(),
        },
    }

    if not truth_sample or not data_type:
        raise DataException("Must specify both truth_sample and desired data_type")

    if truth_sample not in truth_samples:
        raise DataException("This truth sample is not present")

    truth_samples_info = truth_samples[truth_sample]
    if (data_type not in truth_samples_info) or (
        data_type == "callset_truth_mt"
        and not file_exists(
            f"{truth_sample_mt_path(data_source, freeze, truth_sample)}_SUCCESS"
        )
    ):
        raise DataException(
            f"This data type is not present for truth sample: {truth_sample}"
        )

    if data_type == "callset_truth_mt":
        return hl.read_matrix_table(
            truth_sample_mt_path(data_source, freeze, truth_sample)
        )
    else:
        return truth_samples_info[data_type]


def variant_qc_prefix(data_source: str, freeze: int = CURRENT_FREEZE) -> str:
    """
    Returns path to variant QC bucket

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Path to variant QC bucket
    :rtype: str
    """
    if data_source not in DATA_SOURCES:
        raise DataException("This data_source is currently not present")
    if freeze not in FREEZES:
        raise DataException("This freeze is currently not present")

    return f"gs://broad-ukbb/{data_source}.freeze_{freeze}/variant_qc"


def var_annotations_ht_path(
    annotation_type: str, data_source: str, freeze: int = CURRENT_FREEZE
) -> str:
    """
    Get variant-level annotations

    :param str annotation_type: Tag describing variant-level annotations in HT, e.g., "vep", "qc_stats", "call_stats", "vqsr"
    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Path to annotations Table
    :rtype: str
    """
    return f"{variant_qc_prefix(data_source, freeze)}/variant_annotations/{annotation_type}.ht"


def info_ht_path(
    data_source: str, freeze: int = CURRENT_FREEZE, split: bool = True
) -> str:
    """
    Get variant-level annotations

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param bool split: Whether the dataset should be split
    :return: Path to variant info Table
    :rtype: str
    """
    return f'{variant_qc_prefix(data_source, freeze)}/variant_annotations/info{"_split" if split else ""}.ht'


def truth_sample_mt_path(
    truth_sample: str, data_source: str, freeze: int = CURRENT_FREEZE
) -> str:
    """
    Get path to the truth sample MT that is subset from the full callset

    :param str truth_sample: Name of the truth sample
    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Path to callset truth sample MT
    :rtype: str
    """
    return f"{variant_qc_prefix(data_source, freeze)}/truth_samples/{truth_sample}.mt"


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

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param bool adj: Whether to load 'adj' or 'raw'
    :return: Table with RF annotations
    :rtype: Table
    """
    return f'{variant_qc_prefix(data_source, freeze)}/rf/rf_annotated.{"adj" if adj else "raw"}.ht'


def rf_path(
    data_source: str,
    freeze: int = CURRENT_FREEZE,
    data: str = "rf_result",
    run_hash: Optional[str] = None,
) -> str:
    """
    Gets the path to the desired RF data.

    Data can take the following values:
        - 'training': path to the training data for a given run
        - 'model': path to pyspark pipeline RF model
        - 'rf_result' (default): path to HT containing result of RF filtering

    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param str data: One of 'training', 'model' or 'rf_result' (default)
    :param str run_hash: Hash of RF run to load
    :return: Path to desired RF data
    :rtype: str
    """

    extension = "model" if data == "model" else "ht"
    return f"{variant_qc_prefix(data_source, freeze)}/rf/{run_hash}/{data}.{extension}"


def score_quantile_bin_path(
    metric: str,
    data_source: str,
    freeze: int = CURRENT_FREEZE,
    aggregated: bool = False,
) -> str:
    """
    Returns the path to a Table containing RF or VQSR scores and annotated with a bin based on quantiles of the metric scores.

    :param str metric: The score data (RF hash or VQSR) to return
    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :param bool aggregated: Whether to get the aggregated data.
         If True, will return the path to Table grouped by quantile bin that contains aggregated variant counts per bin.
    :return: Path to desired hail Table
    :rtype: str
    """

    return f'{variant_qc_prefix(data_source, freeze)}/score_rankings/{metric}{"_aggregated" if aggregated else ""}.ht'


def binned_concordance_path(
    truth_sample: str, metric: str, data_source: str, freeze: int = CURRENT_FREEZE
):
    """
    Returns the path to a truth sample concordance Table (containing TP, FP, FN) between a truth sample within the callset and the samples truth data grouped by quantile bins of a metric (RF or VQSR scores)
    
    :param str truth_sample: Which truth sample concordance to analyze (e.g., "NA12878" or "syndip")
    :param str metric: One of the evaluation metrics (RF hash or vqsr)
    :param str data_source: One of 'regeneron' or 'broad' 
    :param int freeze: One of the data freezes
    :return: Path to binned truth data concordance Hail Table
    :rtype: str
    """
    return f"{variant_qc_prefix(data_source, freeze)}/rf/{truth_sample}.{metric}.binned_concordance.ht"
