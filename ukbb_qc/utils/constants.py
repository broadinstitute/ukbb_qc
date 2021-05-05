INTERVAL_QC_PARAMETERS = (85, 20)
"""
Parameters used for interval QC (sample percentage cutoff, coverage cutoff).
"""

UKBB_SEXES = ["XX", "XY"]
"""
Sample sexes used in VCF export.

Used to stratify frequency annotations (AC, AN, AF) for each sex in UKBB.
No longer necessary -- this constant also exists in gnomAD repo now:
https://github.com/broadinstitute/gnomad_methods/blob/master/gnomad/resources/grch38/gnomad.py#L29
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
