from gnomad_hail import *
from ukbb_qc.resources import *

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variant_qc_functions")
logger.setLevel(logging.INFO)

def get_rf_runs(data_source: str, freeze: int = CURRENT_FREEZE) -> Dict:
    """

    Loads RF run data from JSON file.

    :param str data_source: 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Dictionary containing the content of the JSON file, or an empty dictionary if the file wasn't found.
    :rtype: dict
    """

    from ukbb_qc.resources import rf_run_hash_path

    json_file = rf_run_hash_path(data_source, freeze)
    if hl.utils.hadoop_exists(json_file):
        with hl.hadoop_open(rf_run_hash_path(data_source, freeze)) as f:
            return json.load(f)
    else:
        logger.warning(f"File {json_file} could not be found. Returning empty RF run hash dict.")
        return {}


