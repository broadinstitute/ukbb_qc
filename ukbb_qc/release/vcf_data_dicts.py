VCF_ANNOTATIONS = {
    "FORMAT_ANNOTATIONS": {
        "Description": "VCF FORMAT field annotations.",
        "sub_annotations": {
            "GT": {"Description": "Genotype."},
            "GQ": {
                "Description": "Phred-scaled confidence that the genotype assignment is correct. Value is the difference between the second lowest PL and the lowest PL (always normalized to 0)."
            },
            "DP": {
                "Description": "Approximate read depth (reads with MQ=255 or with bad mates are filtered)."
            },
            "AD": {
                "Description": "Allelic depths for the ref and alt alleles in the order listed."
            },
            "MIN_DP": {"Description": "Minimum DP observed within the GVCF block."},
            "PGT": {
                "Description": "Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another."
            },
            "PID": {
                "Description": "Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group."
            },
            "PL": {
                "Description": "Normalized, phred-scaled likelihoods for genotypes as defined in the VCF specification."
            },
            "SB": {
                "Description": "Per-sample component statistics which comprise the Fisher's exact test to detect strand bias. Values are: depth of reference allele on forward strand, depth of reference allele on reverse strand, depth of alternate allele on forward strand, depth of alternate allele on reverse strand."
            },
        },
    },
    "FILTER_ANNOTATIONS": {
        "Description": "VCF FILTER field annotations.",
        "sub_annotations": {
            "AC0": {
                "Description": "Allele count is zero after filtering out low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)"
            },
            "InbreedingCoeff": {"Description": "InbreedingCoeff < -0.3."},
            "MonoAllelic": {
                "Description": "Samples are all homozygous reference or all homozygous alternate for the variant."
            },
            "RF": {
                "Description": "Failed random forest filtering thresholds of 0.061026613170892155 for SNPs and 0.06405981151588196 for indels (probabilities of being a true positive variant."
            },
            "PASS": {"Description": "Passed all variant filters."},
        },
    },
    "INFO_ANNOTATIONS": {
        "Description": "VCF INFO field annotations.",
        "sub_annotations": {
            "AC": {"Description": "Alternate allele count."},
            "AN": {"Description": "Total number of alleles."},
            "AF": {"Description": "Alternate allele frequency."},
            "nhomalt": {"Description": "Count of homozygous individuals.",},
            "rf_tp_probability": {
                "Description": "Probability of a called variant being a true variant as determined by random forest model."
            },
            "FS": {
                "Description": "Phred-scaled p-value of Fisher's exact test for strand bias."
            },
            "InbreedingCoeff": {
                "Description": "Inbreeding coefficient, the excess heterozygosity at a variant site, computed as 1 - (the number of heterozygous genotypes)/(the number of heterozygous genotypes expected under Hardy-Weinberg equilibrium).",
            },
            "MQ": {
                "Description": "Root mean square of the mapping quality of reads across all samples."
            },
            "MQRankSum": {
                "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities."
            },
            "QD": {
                "Description": "Variant call confidence normalized by depth of sample reads supporting a variant."
            },
            "ReadPosRankSum": {
                "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference read position bias."
            },
            "SOR": {
                "Description": "Strand bias estimated by the symmetric odds ratio test."
            },
            "POSITIVE_TRAIN_SITE": {
                "Description": "Variant was used to build the positive training set of high-quality variants for VQSR."
            },
            "NEGATIVE_TRAIN_SITE": {
                "Description": "Variant was used to build the negative training set of low-quality variants for VQSR."
            },
            "BaseQRankSum": {
                "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference base qualities."
            },
            "VarDP": {
                "Description": "Depth over variant genotypes (does not include depth of reference samples)."
            },
            "VQSLOD": {
                "Description": "Log-odds ratio of being a true variant versus being a false positive under the trained VQSR Gaussian mixture model."
            },
            "culprit": {
                "Description": "Worst-performing annotation in the VQSR Gaussian mixture model."
            },
            "lcr": {"Description": "Variant falls within a low complexity region."},
            "nonpar": {
                "Description": "Variant (on sex chromosome) falls outside a pseudoautosomal region."
            },
            "segdup": {
                "Description": "Variant falls within a segmental duplication region."
            },
            "rf_positive_label": {
                "Description": "Variant was labelled as a positive example for training of random forest model."
            },
            "rf_negative_label": {
                "Description": "Variant was labelled as a negative example for training of random forest model."
            },
            "rf_label": {"Description": "Random forest training label."},
            "rf_train": {
                "Description": "Variant was used in training random forest model."
            },
            "transmitted_singleton": {
                "Description": "Variant was a callset-wide doubleton that was transmitted within a family from a parent to a child (i.e., a singleton amongst unrelated samples in cohort)."
            },
            "original_alleles": {
                "Description": "Alleles before splitting multiallelics."
            },
            "variant_type": {
                "Description": "Variant type (snv, indel, multi-snv, multi-indel, or mixed)."
            },
            "allele_type": {
                "Description": "Allele type (snv, insertion, deletion, or mixed)."
            },
            "n_alt_alleles": {
                "Description": "Total number of alternate alleles observed at variant locus.",
            },
            "was_mixed": {"Description": "Variant type was mixed."},
            "has_star": {
                "Description": "Variant locus coincides with a spanning deletion (represented by a star) observed elsewhere in the callset."
            },
            "AS_pab_max": {
                "Description": "Maximum p-value over callset for binomial test of observed allele balance for a heterozygous genotype, given expectation of 0.5",
            },
            "outside_capture_region": {
                "Description": "Variant falls outside exome capture regions."
            },
            "sibling_singleton": {
                "Description": "Variant was a callset-wide doubleton that was present only within a sibling pair."
            },
            "AS_FS": {
                "Description": "Allele-specific phred-scaled p-value of Fisher's exact test for strand bias.",
            },
            "AS_MQ": {
                "Description": "Allele-specific root mean square of the mapping quality of reads across all samples.",
            },
            "AS_MQRankSum": {
                "Description": "Allele-specific z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities.",
            },
            "AS_QD": {
                "Description": "Allele-specific variant call confidence normalized by depth of sample reads supporting a variant.",
            },
            "AS_ReadPosRankSum": {
                "Description": "Allele-specific z-score from Wilcoxon rank sum test of alternate vs. reference read position bias.",
            },
            "AS_SOR": {
                "Description": "Allele-specific strand bias estimated by the symmetric odds ratio test.",
            },
            "AS_VarDP": {
                "Description": "Allele-specific depth over variant genotypes (does not include depth of reference samples).",
            },
            "AS_culprit": {
                "Description": "Allele-specific worst-performing annotation in the VQSR Gaussian mixture model.",
            },
            "AS_VQSLOD": {
                "Description": "Allele-specific log-odds ratio of being a true variant versus being a false positive under the trained VQSR Gaussian mixture model.",
            },
            "AC_raw": {
                "Description": "Alternate allele count for samples, before removing low-confidence genotypes.",
            },
            "AN_raw": {
                "Description": "Total number of alleles in samples, before removing low-confidence genotypes.",
            },
            "AF_raw": {
                "Description": "Alternate allele frequency in samples, before removing low-confidence genotypes.",
            },
            "nhomalt_raw": {
                "Description": "Count of homozygous individuals in samples, before removing low-confidence genotypes.",
            },
            "AC_XX": {"Description": "Alternate allele count for XX samples.",},
            "AN_XX": {"Description": "Total number of alleles in XX samples.",},
            "AF_XX": {"Description": "Alternate allele frequency in XX samples.",},
            "nhomalt_XX": {
                "Description": "Count of homozygous individuals in XX samples.",
            },
            "AC_XY": {"Description": "Alternate allele count for XY samples.",},
            "AN_XY": {"Description": "Total number of alleles in XY samples.",},
            "AF_XY": {"Description": "Alternate allele frequency in XY samples.",},
            "nhomalt_XY": {
                "Description": "Count of homozygous individuals in XY samples.",
            },
            "AC_EUR": {
                "Description": "Alternate allele count for samples of European ancestry.",
            },
            "AN_EUR": {
                "Description": "Total number of alleles in samples of European ancestry.",
            },
            "AF_EUR": {
                "Description": "Alternate allele frequency in samples of European ancestry.",
            },
            "nhomalt_EUR": {
                "Description": "Count of homozygous individuals in samples of European ancestry.",
            },
            "AC_CSA": {
                "Description": "Alternate allele count for samples of Central/South Asian ancestry.",
            },
            "AN_CSA": {
                "Description": "Total number of alleles in samples of Central/South Asian ancestry.",
            },
            "AF_CSA": {
                "Description": "Alternate allele frequency in samples of Central/South Asian ancestry.",
            },
            "nhomalt_CSA": {
                "Description": "Count of homozygous individuals in samples of Central/South Asian ancestry.",
            },
            "AC_AFR": {
                "Description": "Alternate allele count for samples of African ancestry.",
            },
            "AN_AFR": {
                "Description": "Total number of alleles in samples of African ancestry.",
            },
            "AF_AFR": {
                "Description": "Alternate allele frequency in samples of African ancestry.",
            },
            "nhomalt_AFR": {
                "Description": "Count of homozygous individuals in samples of African ancestry.",
            },
            "AC_EAS": {
                "Description": "Alternate allele count for samples of East Asian ancestry.",
            },
            "AN_EAS": {
                "Description": "Total number of alleles in samples of East Asian ancestry.",
            },
            "AF_EAS": {
                "Description": "Alternate allele frequency in samples of East Asian ancestry.",
            },
            "nhomalt_EAS": {
                "Description": "Count of homozygous individuals in samples of East Asian ancestry.",
            },
            "AC_MID": {
                "Description": "Alternate allele count for samples of Middle Eastern ancestry.",
            },
            "AN_MID": {
                "Description": "Total number of alleles in samples of Middle Eastern ancestry.",
            },
            "AF_MID": {
                "Description": "Alternate allele frequency in samples of Middle Eastern ancestry.",
            },
            "nhomalt_MID": {
                "Description": "Count of homozygous individuals in samples of Middle Eastern ancestry.",
            },
            "AC_AMR": {
                "Description": "Alternate allele count for samples of Admixed American ancestry.",
            },
            "AN_AMR": {
                "Description": "Total number of alleles in samples of Admixed American ancestry.",
            },
            "AF_AMR": {
                "Description": "Alternate allele frequency in samples of Admixed American ancestry.",
            },
            "nhomalt_AMR": {
                "Description": "Count of homozygous individuals in samples of Admixed American ancestry.",
            },
            "AC_EUR_XX": {
                "Description": "Alternate allele count for XX samples of European ancestry.",
            },
            "AN_EUR_XX": {
                "Description": "Total number of alleles in XX samples of European ancestry.",
            },
            "AF_EUR_XX": {
                "Description": "Alternate allele frequency in XX samples of European ancestry.",
            },
            "nhomalt_EUR_XX": {
                "Description": "Count of homozygous individuals in XX samples of European ancestry.",
            },
            "AC_EUR_XY": {
                "Description": "Alternate allele count for XY samples of European ancestry.",
            },
            "AN_EUR_XY": {
                "Description": "Total number of alleles in XY samples of European ancestry.",
            },
            "AF_EUR_XY": {
                "Description": "Alternate allele frequency in XY samples of European ancestry.",
            },
            "nhomalt_EUR_XY": {
                "Description": "Count of homozygous individuals in XY samples of European ancestry.",
            },
            "AC_CSA_XX": {
                "Description": "Alternate allele count for XX samples of Central/South Asian ancestry.",
            },
            "AN_CSA_XX": {
                "Description": "Total number of alleles in XX samples of Central/South Asian ancestry.",
            },
            "AF_CSA_XX": {
                "Description": "Alternate allele frequency in XX samples of Central/South Asian ancestry.",
            },
            "nhomalt_CSA_XX": {
                "Description": "Count of homozygous individuals in XX samples of Central/South Asian ancestry.",
            },
            "AC_CSA_XY": {
                "Description": "Alternate allele count for XY samples of Central/South Asian ancestry.",
            },
            "AN_CSA_XY": {
                "Description": "Total number of alleles in XY samples of Central/South Asian ancestry.",
            },
            "AF_CSA_XY": {
                "Description": "Alternate allele frequency in XY samples of Central/South Asian ancestry.",
            },
            "nhomalt_CSA_XY": {
                "Description": "Count of homozygous individuals in XY samples of Central/South Asian ancestry.",
            },
            "AC_AFR_XX": {
                "Description": "Alternate allele count for XX samples of African ancestry.",
            },
            "AN_AFR_XX": {
                "Description": "Total number of alleles in XX samples of African ancestry.",
            },
            "AF_AFR_XX": {
                "Description": "Alternate allele frequency in XX samples of African ancestry.",
            },
            "nhomalt_AFR_XX": {
                "Description": "Count of homozygous individuals in XX samples of African ancestry.",
            },
            "AC_AFR_XY": {
                "Description": "Alternate allele count for XY samples of African ancestry.",
            },
            "AN_AFR_XY": {
                "Description": "Total number of alleles in XY samples of African ancestry.",
            },
            "AF_AFR_XY": {
                "Description": "Alternate allele frequency in XY samples of African ancestry.",
            },
            "nhomalt_AFR_XY": {
                "Description": "Count of homozygous individuals in XY samples of African ancestry.",
            },
            "AC_EAS_XX": {
                "Description": "Alternate allele count for XX samples of East Asian ancestry.",
            },
            "AN_EAS_XX": {
                "Description": "Total number of alleles in XX samples of East Asian ancestry.",
            },
            "AF_EAS_XX": {
                "Description": "Alternate allele frequency in XX samples of East Asian ancestry.",
            },
            "nhomalt_EAS_XX": {
                "Description": "Count of homozygous individuals in XX samples of East Asian ancestry.",
            },
            "AC_EAS_XY": {
                "Description": "Alternate allele count for XY samples of East Asian ancestry.",
            },
            "AN_EAS_XY": {
                "Description": "Total number of alleles in XY samples of East Asian ancestry.",
            },
            "AF_EAS_XY": {
                "Description": "Alternate allele frequency in XY samples of East Asian ancestry.",
            },
            "nhomalt_EAS_XY": {
                "Description": "Count of homozygous individuals in XY samples of East Asian ancestry.",
            },
            "AC_MID_XX": {
                "Description": "Alternate allele count for XX samples of Middle Eastern ancestry.",
            },
            "AN_MID_XX": {
                "Description": "Total number of alleles in XX samples of Middle Eastern ancestry.",
            },
            "AF_MID_XX": {
                "Description": "Alternate allele frequency in XX samples of Middle Eastern ancestry.",
            },
            "nhomalt_MID_XX": {
                "Description": "Count of homozygous individuals in XX samples of Middle Eastern ancestry.",
            },
            "AC_MID_XY": {
                "Description": "Alternate allele count for XY samples of Middle Eastern ancestry.",
            },
            "AN_MID_XY": {
                "Description": "Total number of alleles in XY samples of Middle Eastern ancestry.",
            },
            "AF_MID_XY": {
                "Description": "Alternate allele frequency in XY samples of Middle Eastern ancestry.",
            },
            "nhomalt_MID_XY": {
                "Description": "Count of homozygous individuals in XY samples of Middle Eastern ancestry.",
            },
            "AC_AMR_XX": {
                "Description": "Alternate allele count for XX samples of Admixed American ancestry.",
            },
            "AN_AMR_XX": {
                "Description": "Total number of alleles in XX samples of Admixed American ancestry.",
            },
            "AF_AMR_XX": {
                "Description": "Alternate allele frequency in XX samples of Admixed American ancestry.",
            },
            "nhomalt_AMR_XX": {
                "Description": "Count of homozygous individuals in XX samples of Admixed American ancestry.",
            },
            "AC_AMR_XY": {
                "Description": "Alternate allele count for XY samples of Admixed American ancestry.",
            },
            "AN_AMR_XY": {
                "Description": "Total number of alleles in XY samples of Admixed American ancestry.",
            },
            "AF_AMR_XY": {
                "Description": "Alternate allele frequency in XY samples of Admixed American ancestry.",
            },
            "nhomalt_AMR_XY": {
                "Description": "Count of homozygous individuals in XY samples of Admixed American ancestry.",
            },
            "faf95": {
                "Description": "Filtering allele frequency (using Poisson 95% CI).",
            },
            "faf99": {
                "Description": "Filtering allele frequency (using Poisson 99% CI).",
            },
            "faf95_raw": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for samples, before removing low-confidence genotypes.",
            },
            "faf99_raw": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for samples, before removing low-confidence genotypes.",
            },
            "faf95_XX": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XX samples.",
            },
            "faf99_XX": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XX samples.",
            },
            "faf95_XY": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XY samples.",
            },
            "faf99_XY": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XY samples.",
            },
            "faf95_EUR": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for samples of European ancestry.",
            },
            "faf99_EUR": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for samples of European ancestry.",
            },
            "faf95_CSA": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for samples of Central/South Asian ancestry.",
            },
            "faf99_CSA": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for samples of Central/South Asian ancestry.",
            },
            "faf95_AFR": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for samples of African ancestry.",
            },
            "faf99_AFR": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for samples of African ancestry.",
            },
            "faf95_EAS": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for samples of East Asian ancestry.",
            },
            "faf99_EAS": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for samples of East Asian ancestry.",
            },
            "faf95_MID": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for samples of Middle Eastern ancestry.",
            },
            "faf99_MID": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for samples of Middle Eastern ancestry.",
            },
            "faf95_AMR": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for samples of Admixed American ancestry.",
            },
            "faf99_AMR": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for samples of Admixed American ancestry.",
            },
            "faf95_EUR_XX": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XX samples of European ancestry.",
            },
            "faf99_EUR_XX": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XX samples of European ancestry.",
            },
            "faf95_EUR_XY": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XY samples of European ancestry.",
            },
            "faf99_EUR_XY": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XY samples of European ancestry.",
            },
            "faf95_CSA_XX": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XX samples of Central/South Asian ancestry.",
            },
            "faf99_CSA_XX": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XX samples of Central/South Asian ancestry.",
            },
            "faf95_CSA_XY": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XY samples of Central/South Asian ancestry.",
            },
            "faf99_CSA_XY": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XY samples of Central/South Asian ancestry.",
            },
            "faf95_AFR_XX": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XX samples of African ancestry.",
            },
            "faf99_AFR_XX": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XX samples of African ancestry.",
            },
            "faf95_AFR_XY": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XY samples of African ancestry.",
            },
            "faf99_AFR_XY": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XY samples of African ancestry.",
            },
            "faf95_EAS_XX": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XX samples of East Asian ancestry.",
            },
            "faf99_EAS_XX": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XX samples of East Asian ancestry.",
            },
            "faf95_EAS_XY": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XY samples of East Asian ancestry.",
            },
            "faf99_EAS_XY": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XY samples of East Asian ancestry.",
            },
            "faf95_MID_XX": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XX samples of Middle Eastern ancestry.",
            },
            "faf99_MID_XX": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XX samples of Middle Eastern ancestry.",
            },
            "faf95_MID_XY": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XY samples of Middle Eastern ancestry.",
            },
            "faf99_MID_XY": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XY samples of Middle Eastern ancestry.",
            },
            "faf95_AMR_XX": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XX samples of Admixed American ancestry.",
            },
            "faf99_AMR_XX": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XX samples of Admixed American ancestry.",
            },
            "faf95_AMR_XY": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XY samples of Admixed American ancestry.",
            },
            "faf99_AMR_XY": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XY samples of Admixed American ancestry.",
            },
            "age_hist_het_bin_freq": {
                "Description": "Histogram of ages of heterozygous individuals; bin edges are: 30.0|35.0|40.0|45.0|50.0|55.0|60.0|65.0|70.0|75.0|80.0; total number of individuals of any genotype bin: 0|0|7|51796|66083|76337|90827|121484|93585|2424|0|0.",
            },
            "age_hist_het_n_smaller": {
                "Description": "Count of age values falling below lowest histogram bin edge for heterozygous individuals.",
            },
            "age_hist_het_n_larger": {
                "Description": "Count of age values falling above highest histogram bin edge for heterozygous individuals.",
            },
            "age_hist_hom_bin_freq": {
                "Description": "Histogram of ages of homozygous alternate individuals; bin edges are: 30.0|35.0|40.0|45.0|50.0|55.0|60.0|65.0|70.0|75.0|80.0; total number of individuals of any genotype bin: 0|0|7|51796|66083|76337|90827|121484|93585|2424|0|0.",
            },
            "age_hist_hom_n_smaller": {
                "Description": "Count of age values falling below lowest histogram bin edge for homozygous alternate individuals.",
            },
            "age_hist_hom_n_larger": {
                "Description": "Count of age values falling above highest histogram bin edge for homozygous alternate individuals.",
            },
            "popmax": {"Description": "Population with maximum AF."},
            "AC_popmax": {
                "Description": "Allele count in the population with the maximum AF.",
            },
            "AN_popmax": {
                "Description": "Total number of alleles in the population with the maximum AF.",
            },
            "AF_popmax": {
                "Description": "Maximum allele frequency across populations.",
            },
            "nhomalt_popmax": {
                "Description": "Count of homozygous individuals in the population with the maximum allele frequency.",
            },
            "gq_hist_alt_bin_freq": {
                "Description": "Histogram for GQ in heterozygous individuals calculated on high quality genotypes; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.",
            },
            "gq_hist_alt_n_smaller": {
                "Description": "Count of GQ values falling below lowest histogram bin edge GQ in heterozygous individuals calculated on high quality genotypes.",
            },
            "gq_hist_alt_n_larger": {
                "Description": "Count of GQ values falling above highest histogram bin edge GQ in heterozygous individuals calculated on high quality genotypes.",
            },
            "gq_hist_all_bin_freq": {
                "Description": "Histogram for GQ calculated on high quality genotypes; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.",
            },
            "gq_hist_all_n_smaller": {
                "Description": "Count of GQ values falling below lowest histogram bin edge GQ calculated on high quality genotypes.",
            },
            "gq_hist_all_n_larger": {
                "Description": "Count of GQ values falling above highest histogram bin edge GQ calculated on high quality genotypes.",
            },
            "dp_hist_alt_bin_freq": {
                "Description": "Histogram for DP in heterozygous individuals calculated on high quality genotypes; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.",
            },
            "dp_hist_alt_n_smaller": {
                "Description": "Count of DP values falling below lowest histogram bin edge DP in heterozygous individuals calculated on high quality genotypes.",
            },
            "dp_hist_alt_n_larger": {
                "Description": "Count of DP values falling above highest histogram bin edge DP in heterozygous individuals calculated on high quality genotypes.",
            },
            "dp_hist_all_bin_freq": {
                "Description": "Histogram for DP calculated on high quality genotypes; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.",
            },
            "dp_hist_all_n_smaller": {
                "Description": "Count of DP values falling below lowest histogram bin edge DP calculated on high quality genotypes.",
            },
            "dp_hist_all_n_larger": {
                "Description": "Count of DP values falling above highest histogram bin edge DP calculated on high quality genotypes.",
            },
            "ab_hist_alt_bin_freq": {
                "Description": "Histogram for AB in heterozygous individuals calculated on high quality genotypes; bin edges are: 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00.",
            },
            "ab_hist_alt_n_smaller": {
                "Description": "Count of AB values falling below lowest histogram bin edge AB in heterozygous individuals calculated on high quality genotypes.",
            },
            "ab_hist_alt_n_larger": {
                "Description": "Count of AB values falling above highest histogram bin edge AB in heterozygous individuals calculated on high quality genotypes.",
            },
            "gq_hist_alt_raw_bin_freq": {
                "Description": "Histogram for GQ in heterozygous individuals; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.",
            },
            "gq_hist_alt_raw_n_smaller": {
                "Description": "Count of GQ values falling below lowest histogram bin edge GQ in heterozygous individuals.",
            },
            "gq_hist_alt_raw_n_larger": {
                "Description": "Count of GQ values falling above highest histogram bin edge GQ in heterozygous individuals.",
            },
            "gq_hist_all_raw_bin_freq": {
                "Description": "Histogram for GQ; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.",
            },
            "gq_hist_all_raw_n_smaller": {
                "Description": "Count of GQ values falling below lowest histogram bin edge GQ.",
            },
            "gq_hist_all_raw_n_larger": {
                "Description": "Count of GQ values falling above highest histogram bin edge GQ.",
            },
            "dp_hist_alt_raw_bin_freq": {
                "Description": "Histogram for DP in heterozygous individuals; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.",
            },
            "dp_hist_alt_raw_n_smaller": {
                "Description": "Count of DP values falling below lowest histogram bin edge DP in heterozygous individuals.",
            },
            "dp_hist_alt_raw_n_larger": {
                "Description": "Count of DP values falling above highest histogram bin edge DP in heterozygous individuals.",
            },
            "dp_hist_all_raw_bin_freq": {
                "Description": "Histogram for DP; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.",
            },
            "dp_hist_all_raw_n_smaller": {
                "Description": "Count of DP values falling below lowest histogram bin edge DP.",
            },
            "dp_hist_all_raw_n_larger": {
                "Description": "Count of DP values falling above highest histogram bin edge DP.",
            },
            "ab_hist_alt_raw_bin_freq": {
                "Description": "Histogram for AB in heterozygous individuals; bin edges are: 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00.",
            },
            "ab_hist_alt_raw_n_smaller": {
                "Description": "Count of AB values falling below lowest histogram bin edge AB in heterozygous individuals.",
            },
            "ab_hist_alt_raw_n_larger": {
                "Description": "Count of AB values falling above highest histogram bin edge AB in heterozygous individuals.",
            },
            "fail_interval_qc": {
                "Description": "Variant falls within a region where less than 85.0% of samples had a mean coverage of 20X on autosomes and 10X on sex chromosomes."
            },
            "vep": {
                "Description": "Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info."
            },
        },
    },
}
