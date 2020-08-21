INTERVAL_QC_PARAMETERS = (85, 20)
"""
Parameters used for interval QC (sample percentage cutoff, coverage cutoff).
"""

SEXES_UKBB = ["XX", "XY"]
"""
Sample sexes used in VCF export.

Used to stratify frequency annotations (AC, AN, AF) for each sex in UKBB.
"""

ANNOTATIONS_HISTS = {
    "AS_FS": (0, 50, 50),
    "InbreedingCoeff": (-0.25, 0.25, 50),
    "AS_MQ": (0, 80, 40),
    "AS_MQRankSum": (-15, 15, 60),
    "AS_QD": (0, 40, 40),
    "AS_ReadPosRankSum": (-15, 15, 60),
    "AS_SOR": (0, 10, 50),
    "AS_BaseQRankSum": (-15, 15, 60),
    "AS_VarDP": (1, 9, 32),
    "AS_VQSLOD": (-30, 30, 60),
    "rf_tp_probability": (0, 1, 50),
    "AS_pab_max": (0, 1, 50),
}
"""
Dictionary of metric names and their histogram values (start, end, bins).
"""
