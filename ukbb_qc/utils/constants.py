# UKB constants
GROUPS = ["adj", "raw"]
SEXES = ["male", "female"]
POPS = [
    "0",
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10",
    "11",
    "12",
    "13",
    "afr",
    "amr",
    "asj",
    "eas",
    "fin",
    "nfe",
    "oth",
    "sas",
]
FAF_POPS = ["0", "5", "6", "13", "nfe"]

# gnomAD constants
GNOMAD_POPS = ["afr", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"]
GNOMAD_FAF_POPS = ["afr", "amr", "eas", "nfe", "sas"]
GNOMAD_NFE_SUBPOPS = ["onf", "bgr", "swe", "nwe", "seu", "est"]
GNOMAD_EAS_SUBPOPS = ["kor", "oea", "jpn"]

# constants for both UKB and gnomAD
HISTS = ["gq_hist_alt", "gq_hist_all", "dp_hist_alt", "dp_hist_all", "ab_hist_alt"]
SORT_ORDER = ["popmax", "group", "pop", "subpop", "sex"]

# pop names dict (used in VCF export)
pop_names = {
    "0": "hybrid population inference cluster 0",
    "1": "hybrid population inference cluster 1",
    "2": "hybrid population inference cluster 2",
    "3": "hybrid population inference cluster 3",
    "4": "hybrid population inference cluster 4",
    "5": "hybrid population inference cluster 5",
    "6": "hybrid population inference cluster 6",
    "7": "hybrid population inference cluster 7",
    "8": "hybrid population inference cluster 8",
    "9": "hybrid population inference cluster 9",
    "10": "hybrid population inference cluster 10",
    "11": "hybrid population inference cluster 11",
    "12": "hybrid population inference cluster 12",
    "13": "hybrid population inference cluster 13",
    "oth": "Other",
    "afr": "African-American/African",
    "ami": "Amish",
    "amr": "Latino",
    "eas": "East Asian",
    "fin": "Finnish",
    "eur": "European",
    "nfe": "Non-Finnish European",
    "sas": "South Asian",
    "mde": "Middle Eastern",
    "asj": "Ashkenazi Jewish",
    "uniform": "Uniform",
    "sas_non_consang": "South Asian (F < 0.05)",
    "consanguineous": "South Asian (F > 0.05)",
    "exac": "ExAC",
    "bgr": "Bulgarian (Eastern European)",
    "deu": "German",
    "est": "Estonian",
    "esp": "Spanish",
    "gbr": "British",
    "nwe": "North-Western European",
    "seu": "Southern European",
    "ita": "Italian",
    "swe": "Swedish",
    "chn": "Chinese",
    "kor": "Korean",
    "hkg": "Hong Kong",
    "sgp": "Singaporean",
    "twn": "Taiwanese",
    "jpn": "Japanese",
    "oea": "Other East Asian",
    "oeu": "Other European",
    "onf": "Other Non-Finnish European",
    "unk": "Unknown",
}

# INFO dict used in VCF export for variant annotations
INFO_DICT = {
    "FS": {
        "Description": "Phred-scaled p-value of Fisher's exact test for strand bias"
    },
    "InbreedingCoeff": {
        "Description": "Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation"
    },
    "MQ": {
        "Description": "Root mean square of the mapping quality of reads across all samples"
    },
    "MQRankSum": {
        "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities"
    },
    "QD": {
        "Description": "Variant call confidence normalized by depth of sample reads supporting a variant"
    },
    "ReadPosRankSum": {
        "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference read position bias"
    },
    "SOR": {"Description": "Strand bias estimated by the symmetric odds ratio test"},
    "VQSR_POSITIVE_TRAIN_SITE": {
        "Description": "Variant was used to build the positive training set of high-quality variants for VQSR"
    },
    "VQSR_NEGATIVE_TRAIN_SITE": {
        "Description": "Variant was used to build the negative training set of low-quality variants for VQSR"
    },
    "AS_BaseQRankSum": {
        "Description": "Allele specific Z-score from Wilcoxon rank sum test of alternate vs. reference base qualities"
    },
    "ClippingRankSum": {
        "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference number of hard clipped bases"
    },
    "VarDP": {
        "Description": "Depth over variant genotypes (does not include depth of reference samples)"
    },
    "AS_VarDP": {
        "Description": "Allele specific depth over variant genotypes (does not include depth of reference samples)"
    },
    "AS_VQSLOD": {
        "Description": "Allele specific log-odds ratio of being a true variant versus being a false positive under the trained VQSR Gaussian mixture model"
    },
    "AS_VQSR_culprit": {
        "Description": "Allele specific worst-performing annotation in the VQSR Gaussian mixture model"
    },
    "lcr": {"Description": "Variant falls within a low complexity region"},
    "fail_interval_qc": {
        "Description": "Variant falls within a region where less than 85% of samples had a mean coverage of 20X"
    },
    "nonpar": {
        "Description": "Variant (on sex chromosome) falls outside a pseudoautosomal region"
    },
    "rf_positive_label": {
        "Description": "Variant was labelled as a positive example for training of random forest model"
    },
    "rf_negative_label": {
        "Description": "Variant was labelled as a negative example for training of random forest model"
    },
    "rf_label": {"Description": "Random forest training label"},  # export?
    "rf_train": {
        "Description": "Variant was used in training random forest model"
    },  # export?
    "rf_tp_probability": {
        "Description": "Random forest prediction probability for a site being a true variant"
    },
    "transmitted_singleton": {
        "Description": "Variant was a callset-wide doubleton that was transmitted within a family (i.e., a singleton amongst unrelated sampes in cohort)"
    },
    "original_alleles": {"Description": "Alleles before splitting multiallelics"},
    "variant_type": {
        "Description": "Variant type (snv, indel, multi-snv, multi-indel, or mixed)"
    },
    "allele_type": {
        "Number": "A",
        "Description": "Allele type (snv, ins, del, or mixed)",
    },
    "n_alt_alleles": {
        "Number": "A",
        "Description": "Total number of alternate alleles observed at variant locus",
    },
    "was_mixed": {"Description": "Variant type was mixed"},
    "has_star": {
        "Description": "Variant locus coincides with a spanning deletion (represented by a star) observed elsewhere in the callset"
    },
    "pab_max": {
        "Number": "A",
        "Description": "Maximum p-value over callset for binomial test of observed allele balance for a heterozygous genotype, given expectation of AB=0.5",
    },
}

# FORMAT dict used in VCF export for MatrixTable entries
FORMAT_DICT = {
    "GT": {"Description": "Genotype", "Number": "1", "Type": "String"},
    "AD": {
        "Description": "Allelic depths for the ref and alt alleles in the order listed",
        "Number": "R",
        "Type": "Integer",
    },
    "DP": {
        "Description": "Approximate read depth (reads with MQ=255 or with bad mates are filtered)",
        "Number": "1",
        "Type": "Integer",
    },
    "GQ": {"Description": "Genotype Quality", "Number": "1", "Type": "Integer"},
    "MIN_DP": {
        "Description": "Minimum DP observed within the GVCF block",
        "Number": "1",
        "Type": "Integer",
    },
    "PGT": {
        "Description": "Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another",
        "Number": "1",
        "Type": "String",
    },
    "PID": {
        "Description": "Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group",
        "Number": "1",
        "Type": "String",
    },
    "PL": {
        "Description": "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification",
        "Number": "G",
        "Type": "Integer",
    },
    "SB": {
        "Description": "Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.",
        "Number": "4",
        "Type": "Integer",
    },
}
