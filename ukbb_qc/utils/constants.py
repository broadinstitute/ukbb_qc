INTERVAL_QC_PARAMETERS = (85, 20)
"""
Parameters used for interval QC (sample percentage cutoff, coverage cutoff).
"""

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
        "Description": "Allele-specific Z-score from Wilcoxon rank sum test of alternate vs. reference base qualities"
    },
    "ClippingRankSum": {
        "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference number of hard clipped bases"
    },
    "VarDP": {
        "Description": "Depth over variant genotypes (does not include depth of reference samples)"
    },
    "AS_VarDP": {
        "Description": "Allele-specific depth over variant genotypes (does not include depth of reference samples)"
    },
    "AS_VQSLOD": {
        "Description": "Allele-specific log-odds ratio of being a true variant versus being a false positive under the trained VQSR Gaussian mixture model"
    },
    "AS_VQSR_culprit": {
        "Description": "Allele-specific worst-performing annotation in the VQSR Gaussian mixture model"
    },
    "lcr": {"Description": "Variant falls within a low complexity region"},
    "fail_interval_qc": {
        "Description": f"Variant falls within a region where less than {INTERVAL_QC_PARAMETERS[0]}% of samples had a mean coverage of {INTERVAL_QC_PARAMETERS[1]}X"
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
    "rf_label": {"Description": "Random forest training label"},
    "rf_train": {"Description": "Variant was used in training random forest model"},
    "rf_tp_probability": {
        "Description": "Probability of a called variant being a true variant"
    },
    "transmitted_singleton": {
        "Description": "Variant was a callset-wide doubleton that was transmitted within a family from a parent to a child (i.e., a singleton amongst unrelated samples in cohort)"
    },
    "original_alleles": {"Description": "Alleles before splitting multiallelics"},
    "variant_type": {
        "Description": "Variant type (snv, indel, multi-snv, multi-indel, or mixed)"
    },
    "allele_type": {
        "Number": "A",
        "Description": "Allele type (snv, insertion, deletion, or mixed)",
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
"""
Dictionary used during VCF export to export variant annotations
"""

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
    "GQ": {
        "Description": "Phred-scaled confidence that the genotype assignment is correct. Derived from the PLs of the two most likely genotypes",
        "Number": "1",
        "Type": "Integer",
    },
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
"""
Dictionary used during VCF export to export MatrixTable entries
"""
