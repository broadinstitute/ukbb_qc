INTERVAL_QC_PARAMETERS = (85, 20)
"""
Parameters used for interval QC (sample percentage cutoff, coverage cutoff).
"""

SEXES_UKBB = ["XX", "XY"]
"""
Sample sexes used in VCF export.

Used to stratify frequency annotations (AC, AN, AF) for each sex in UKBB.
"""

UKBB_POPS = {
    "EUR": "European",
    "CSA": "Central/South Asian",
    "AFR": "African",
    "EAS": "East Asian",
    "MID": "Middle Eastern",
    "AMR": "Admixed American",
}
"""
Dictionary of pan-ancestry population labels and their descriptions.

Information taken from:
 https://pan.ukbb.broadinstitute.org/docs/technical-overview
"""
