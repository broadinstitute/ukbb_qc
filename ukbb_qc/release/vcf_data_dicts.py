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
                "Description": "Failed random forest filtering thresholds of 0.4 for SNPs and 0.2 for indels (probabilities of being a true positive variant."
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
            "rf_tp_probability": {
                "Description": "Probability of a called variant being a true variant as determined by random forest model."
            },
            "allele_type": {
                "Description": "Allele type (snv, insertion, deletion, or mixed).",
            },
            "has_star": {
                "Description": "Variant locus coincides with a spanning deletion (represented by a star) observed elsewhere in the callset",
            },
            "n_alt_alleles": {
                "Description": "Total number of alternate alleles observed at variant locus.",
            },
            "original_alleles": {
                "Description": "Alleles before splitting multiallelics."
            },
            "variant_type": {
                "Description": "Variant type (snv, indel, multi-snv, multi-indel, or mixed).",
            },
            "was_mixed": {"Description": "Variant type was mixed."},
            "lcr": {"Description": "Variant falls within a low complexity region."},
            "nonpar": {
                "Description": "Variant (on sex chromosome) falls outside a pseudoautosomal region."
            },
            "segdup": {
                "Description": "Variant falls within a segmental duplication region."
            },
            "fail_interval_qc": {
                "Description": "Variant falls within a region where less than 85.0% of samples had a mean coverage of 20X on autosomes and 10X on sex chromosomes."
            },
            "in_capture_region": {
                "Description": "Variant falls within an exome capture region."
            },
            "FS": {
                "Description": "Phred-scaled p-value of Fisher's exact test for strand bias."
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
            "VarDP": {
                "Description": "Depth over variant genotypes (does not include depth of reference samples)."
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
            "AS_VQSLOD": {
                "Description": "Allele-specific log-odds ratio of being a true variant versus being a false positive under the trained VQSR Gaussian mixture model."
            },
            "AS_culprit": {
                "Description": "Allele-specific worst-performing annotation in the VQSR Gaussian mixture model."
            },
            "NEGATIVE_TRAIN_SITE": {
                "Description": "Variant was used to build the negative training set of low-quality variants for VQSR."
            },
            "POSITIVE_TRAIN_SITE": {
                "Description": "Variant was used to build the positive training set of high-quality variants for VQSR."
            },
            "AS_FS": {
                "Description": "Allele-specific phred-scaled p-value of Fisher's exact test for strand bias."
            },
            "AS_MQ": {
                "Description": "Allele-specific root mean square of the mapping quality of reads across all samples."
            },
            "AS_MQRankSum": {
                "Description": "Allele-specific z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities."
            },
            "AS_pab_max": {
                "Description": (
                    "Maximum p-value over callset for binomial test of observed allele balance for a heterozygous "
                    "genotype, given expectation of 0.5."
                )
            },
            "AS_QD": {
                "Description": "Allele-specific variant call confidence normalized by depth of sample reads supporting a variant."
            },
            "AS_ReadPosRankSum": {
                "Description": "Allele-specific z-score from Wilcoxon rank sum test of alternate vs. reference read position bias."
            },
            "AS_SOR": {
                "Description": "Allele-specific strand bias estimated by the symmetric odds ratio test."
            },
            "AS_VarDP": {
                "Description": "Allele-specific depth over variant genotypes (does not include depth of reference samples)."
            },
            "InbreedingCoeff": {
                "Description": (
                    "Inbreeding coefficient, the excess heterozygosity at a variant site, computed as 1 - (the number "
                    "of heterozygous genotypes)/(the number of heterozygous genotypes expected under Hardy-Weinberg "
                    "equilibrium)."
                )
            },
            "sibling_singleton": {
                "Description": "Variant was a callset-wide doubleton that was present only within a sibling pair."
            },
            "gq_hist_alt_bin_freq": {
                "Description": (
                    "Histogram for GQ in heterozygous individuals calculated on high quality genotypes; "
                    "bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                )
            },
            "gq_hist_alt_raw_bin_freq": {
                "Description": (
                    "Histogram for GQ in heterozygous individuals; bin edges are: "
                    "0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                )
            },
            "gq_hist_all_bin_freq": {
                "Description": (
                    "Histogram for GQ calculated on high quality genotypes; "
                    "bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                )
            },
            "gq_hist_all_raw_bin_freq": {
                "Description": (
                    "Histogram for GQ; bin edges are: "
                    "0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                )
            },
            "gq_hist_all_raw_bin_freq": {
                "Description": (
                    "Histogram for GQ; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                )
            },
            "dp_hist_alt_bin_freq": {
                "Description": (
                    "Histogram for DP in heterozygous individuals calculated on high quality genotypes; "
                    "bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                )
            },
            "dp_hist_alt_n_larger": {
                "Description": (
                    "Count of DP values falling above highest histogram bin edge DP in heterozygous individuals "
                    "calculated on high quality genotypes."
                )
            },
            "dp_hist_alt_raw_bin_freq": {
                "Description": (
                    "Histogram for DP in heterozygous individuals; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                )
            },
            "dp_hist_alt_raw_n_larger": {
                "Description": (
                    "Count of DP values falling above highest histogram bin edge DP in heterozygous individuals."
                )
            },
            "dp_hist_all_bin_freq": {
                "Description": (
                    "Histogram for DP calculated on high quality genotypes; bin edges are: "
                    "0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                )
            },
            "dp_hist_all_n_larger": {
                "Description": (
                    "Count of DP values falling above highest histogram bin edge DP calculated on high quality genotypes."
                )
            },
            "dp_hist_all_raw_bin_freq": {
                "Description": (
                    "Histogram for DP; bin edges are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                )
            },
            "dp_hist_all_raw_n_larger": {
                "Description": (
                    "Count of DP values falling above highest histogram bin edge DP."
                )
            },
            "ab_hist_alt_bin_freq": {
                "Description": (
                    "Histogram for AB in heterozygous individuals calculated on high quality genotypes; bin edges are: "
                    "0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00."
                )
            },
            "ab_hist_alt_raw_bin_freq": {
                "Description": (
                    "Histogram for AB in heterozygous individuals; bin edges are: "
                    "0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00."
                )
            },
            "nhomalt": {"Description": "Count of homozygous individuals."},
            "AC_raw": {
                "Description": "Alternate allele count for samples before removing low-confidence genotypes."
            },
            "AN": {
                "Description": "Total number of alleles in samples before removing low-confidence genotypes."
            },
            "AF": {
                "Description": "Alternate allele frequency in samples before removing low-confidence genotypes."
            },
            "nhomalt_raw": {
                "Description": "Count of homozygous individuals in samples before removing low-confidence genotypes."
            },
            "AC_afr": {
                "Description": "Alternate allele count for samples of African/African-American ancestry."
            },
            "AN_afr": {
                "Description": "Total number of alleles in samples of African/African-American ancestry."
            },
            "AF_afr": {
                "Description": "Alternate allele frequency in samples of African/African-American ancestry."
            },
            "nhomalt_afr": {
                "Description": "Count of homozygous individuals in samples of African/African-American ancestry."
            },
            "AC_amr": {
                "Description": "Alternate allele count for samples of Latino ancestry."
            },
            "AN_amr": {
                "Description": "Total number of alleles in samples of Latino ancestry."
            },
            "AF_amr": {
                "Description": "Alternate allele frequency in samples of Latino ancestry."
            },
            "nhomalt_amr": {
                "Description": "Count of homozygous individuals in samples of Latino ancestry."
            },
            "AC_asj": {
                "Description": "Alternate allele count for samples of Ashkenazi Jewish ancestry."
            },
            "AN_asj": {
                "Description": "Total number of alleles in samples of Ashkenazi Jewish ancestry."
            },
            "AF_asj": {
                "Description": "Alternate allele frequency in samples of Ashkenazi Jewish ancestry."
            },
            "nhomalt_asj": {
                "Description": "Count of homozygous individuals in samples of Ashkenazi Jewish ancestry."
            },
            "AC_eas": {
                "Description": "Alternate allele count for samples of East Asian ancestry."
            },
            "AN_eas": {
                "Description": "Total number of alleles in samples of East Asian ancestry."
            },
            "AF_eas": {
                "Description": "Alternate allele frequency in samples of East Asian ancestry."
            },
            "nhomalt_eas": {
                "Description": "Count of homozygous individuals in samples of East Asian ancestry."
            },
            "AC_fin": {
                "Description": "Alternate allele count for samples of Finnish ancestry."
            },
            "AN_fin": {
                "Description": "Total number of alleles in samples of Finnish ancestry."
            },
            "AF_fin": {
                "Description": "Alternate allele frequency in samples of Finnish ancestry."
            },
            "nhomalt_fin": {
                "Description": "Count of homozygous individuals in samples of Finnish ancestry."
            },
            "AC_nfe": {
                "Description": "Alternate allele count for samples of Non-Finnish European ancestry."
            },
            "AN_nfe": {
                "Description": "Total number of alleles in samples of Non-Finnish European ancestry."
            },
            "AF_nfe": {
                "Description": "Alternate allele frequency in samples of Non-Finnish European ancestry."
            },
            "nhomalt_nfe": {
                "Description": "Count of homozygous individuals in samples of Non-Finnish European ancestry."
            },
            "AC_oth": {
                "Description": "Alternate allele count for samples of Other ancestry."
            },
            "AN_oth": {
                "Description": "Total number of alleles in samples of Other ancestry."
            },
            "AF_oth": {
                "Description": "Alternate allele frequency in samples of Other ancestry."
            },
            "nhomalt_oth": {
                "Description": "Count of homozygous individuals in samples of Other ancestry."
            },
            "AC_sas": {
                "Description": "Alternate allele count for samples of South Asian ancestry."
            },
            "AN_sas": {
                "Description": "Total number of alleles in samples of South Asian ancestry."
            },
            "AF_sas": {
                "Description": "Alternate allele frequency in samples of South Asian ancestry."
            },
            "nhomalt_sas": {
                "Description": "Count of homozygous individuals in samples of South Asian ancestry."
            },
            "AC_XX": {"Description": "Alternate allele count for XX samples."},
            "AN_XX": {"Description": "Total number of alleles in XX samples."},
            "AF_XX": {"Description": "Alternate allele frequency in XX samples."},
            "nhomalt_XX": {
                "Description": "Count of homozygous individuals in XX samples."
            },
            "AC_XY": {"Description": "Alternate allele count for XY samples."},
            "AN_XY": {"Description": "Total number of alleles in XY samples."},
            "AF_XY": {"Description": "Alternate allele frequency in XY samples."},
            "nhomalt_XY": {
                "Description": "Count of homozygous individuals in XY samples."
            },
            "AC_afr_XX": {
                "Description": "Alternate allele count for XX samples of African/African-American ancestry."
            },
            "AN_afr_XX": {
                "Description": "Total number of alleles in XX samples of African/African-American ancestry."
            },
            "AF_afr_XX": {
                "Description": "Alternate allele frequency in XX samples of African/African-American ancestry."
            },
            "nhomalt_afr_XX": {
                "Description": "Count of homozygous individuals in XX samples of African/African-American ancestry."
            },
            "AC_afr_XY": {
                "Description": "Alternate allele count for XY samples of African/African-American ancestry."
            },
            "AN_afr_XY": {
                "Description": "Total number of alleles in XY samples of African/African-American ancestry."
            },
            "AF_afr_XY": {
                "Description": "Alternate allele frequency in XY samples of African/African-American ancestry."
            },
            "nhomalt_afr_XY": {
                "Description": "Count of homozygous individuals in XY samples of African/African-American ancestry."
            },
            "AC_amr_XX": {
                "Description": "Alternate allele count for XX samples of Latino ancestry."
            },
            "AN_amr_XX": {
                "Description": "Total number of alleles in XX samples of Latino ancestry."
            },
            "AF_amr_XX": {
                "Description": "Alternate allele frequency in XX samples of Latino ancestry."
            },
            "nhomalt_amr_XX": {
                "Description": "Count of homozygous individuals in XX samples of Latino ancestry."
            },
            "AC_amr_XY": {
                "Description": "Alternate allele count for XY samples of Latino ancestry."
            },
            "AN_amr_XY": {
                "Description": "Total number of alleles in XY samples of Latino ancestry."
            },
            "AF_amr_XY": {
                "Description": "Alternate allele frequency in XY samples of Latino ancestry."
            },
            "nhomalt_amr_XY": {
                "Description": "Count of homozygous individuals in XY samples of Latino ancestry."
            },
            "AC_asj_XX": {
                "Description": "Alternate allele count for XX samples of Ashkenazi Jewish ancestry."
            },
            "AN_asj_XX": {
                "Description": "Total number of alleles in XX samples of Ashkenazi Jewish ancestry."
            },
            "AF_asj_XX": {
                "Description": "Alternate allele frequency in XX samples of Ashkenazi Jewish ancestry."
            },
            "nhomalt_asj_XX": {
                "Description": "Count of homozygous individuals in XX samples of Ashkenazi Jewish ancestry."
            },
            "AC_asj_XY": {
                "Description": "Alternate allele count for XY samples of Ashkenazi Jewish ancestry."
            },
            "AN_asj_XY": {
                "Description": "Total number of alleles in XY samples of Ashkenazi Jewish ancestry."
            },
            "AF_asj_XY": {
                "Description": "Alternate allele frequency in XY samples of Ashkenazi Jewish ancestry."
            },
            "nhomalt_asj_XY": {
                "Description": "Count of homozygous individuals in XY samples of Ashkenazi Jewish ancestry."
            },
            "AC_eas_XX": {
                "Description": "Alternate allele count for XX samples of East Asian ancestry."
            },
            "AN_eas_XX": {
                "Description": "Total number of alleles in XX samples of East Asian ancestry."
            },
            "AF_eas_XX": {
                "Description": "Alternate allele frequency in XX samples of East Asian ancestry."
            },
            "nhomalt_eas_XX": {
                "Description": "Count of homozygous individuals in XX samples of East Asian ancestry."
            },
            "AC_eas_XY": {
                "Description": "Alternate allele count for XY samples of East Asian ancestry."
            },
            "AN_eas_XY": {
                "Description": "Total number of alleles in XY samples of East Asian ancestry."
            },
            "AF_eas_XY": {
                "Description": "Alternate allele frequency in XY samples of East Asian ancestry."
            },
            "nhomalt_eas_XY": {
                "Description": "Count of homozygous individuals in XY samples of East Asian ancestry."
            },
            "AC_fin_XX": {
                "Description": "Alternate allele count for XX samples of Finnish ancestry."
            },
            "AN_fin_XX": {
                "Description": "Total number of alleles in XX samples of Finnish ancestry."
            },
            "AF_fin_XX": {
                "Description": "Alternate allele frequency in XX amples of Finnish ancestry."
            },
            "nhomalt_fin_XX": {
                "Description": "Count of homozygous individuals in XX samples of Finnish ancestry."
            },
            "AC_fin_XY": {
                "Description": "Alternate allele count for XY samples of Finnish ancestry."
            },
            "AN_fin_XY": {
                "Description": "Total number of alleles in XY samples of Finnish ancestry."
            },
            "AF_fin_XY": {
                "Description": "Alternate allele frequency in XY samples of Finnish ancestry."
            },
            "nhomalt_fin_XY": {
                "Description": "Count of homozygous individuals in XY samples of Finnish ancestry."
            },
            "AC_nfe_XX": {
                "Description": "Alternate allele count for XX samples of Non-Finnish European ancestry."
            },
            "AN_nfe_XX": {
                "Description": "Total number of alleles in XX samples of Non-Finnish European ancestry."
            },
            "AF_nfe_XX": {
                "Description": "Alternate allele frequency in XX samples of Non-Finnish European ancestry."
            },
            "nhomalt_nfe_XX": {
                "Description": "Count of homozygous individuals in XX samples of Non-Finnish European ancestry."
            },
            "AC_nfe_XY": {
                "Description": "Alternate allele count for XY samples of Non-Finnish European ancestry."
            },
            "AN_nfe_XY": {
                "Description": "Total number of alleles in XY samples of Non-Finnish European ancestry."
            },
            "AF_nfe_XY": {
                "Description": "Alternate allele frequency in XY samples of Non-Finnish European ancestry."
            },
            "nhomalt_nfe_XY": {
                "Description": "Count of homozygous individuals in XY samples of Non-Finnish European ancestry."
            },
            "AC_oth_XX": {
                "Description": "Alternate allele count for XX samples of Other ancestry."
            },
            "AN_oth_XX": {
                "Description": "Total number of alleles in XX samples of Other ancestry."
            },
            "AF_oth_XX": {
                "Description": "Alternate allele frequency in XX samples of Other ancestry."
            },
            "nhomalt_oth_XX": {
                "Description": "Count of homozygous individuals in XX samples of Other ancestry."
            },
            "AC_oth_XY": {
                "Description": "Alternate allele count for XY samples of Other ancestry."
            },
            "AN_oth_XY": {
                "Description": "Total number of alleles in XY samples of Other ancestry."
            },
            "AF_oth_XY": {
                "Description": "Alternate allele frequency in XY samples of Other ancestry."
            },
            "nhomalt_oth_XY": {
                "Description": "Count of homozygous individuals in XY samples of Other ancestry."
            },
            "AC_sas_XX": {
                "Description": "Alternate allele count for XX samples of South Asian ancestry."
            },
            "AN_sas_XX": {
                "Description": "Total number of alleles in XX samples of South Asian ancestry."
            },
            "AF_sas_XX": {
                "Description": "Alternate allele frequency in XX samples of South Asian ancestry."
            },
            "nhomalt_sas_XX": {
                "Description": "Count of homozygous individuals in XX samples of South Asian ancestry."
            },
            "AC_sas_XY": {
                "Description": "Alternate allele count for XY samples of South Asian ancestry."
            },
            "AN_sas_XY": {
                "Description": "Total number of alleles in XY samples of South Asian ancestry."
            },
            "AF_sas_XY": {
                "Description": "Alternate allele frequency in XY samples of South Asian ancestry."
            },
            "nhomalt_sas_XY": {
                "Description": "Count of homozygous individuals in XY samples of South Asian ancestry."
            },
            "faf95": {
                "Description": "Filtering allele frequency (using Poisson 95% CI)."
            },
            "faf99": {
                "Description": "Filtering allele frequency (using Poisson 99% CI)."
            },
            "faf95_afr": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for samples of African/African-American ancestry."
            },
            "faf99_afr": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for samples of African/African-American ancestry."
            },
            "faf95_amr": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for samples of Latino ancestry."
            },
            "faf99_amr": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for samples of Latino ancestry."
            },
            "faf95_eas": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for samples of East Asian ancestry."
            },
            "faf99_eas": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for samples of East Asian ancestry."
            },
            "faf95_nfe": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for samples of Non-Finnish European ancestry."
            },
            "faf99_nfe": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for samples of Non-Finnish European ancestry."
            },
            "faf95_sas": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for samples of South Asian ancestry."
            },
            "faf99_sas": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for samples of South Asian ancestry."
            },
            "faf95_XX": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XX samples."
            },
            "faf99_XX": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XX samples."
            },
            "faf95_XY": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XY samples."
            },
            "faf99_XY": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XY samples."
            },
            "faf95_afr_XX": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XX samples of African/African-American ancestry."
            },
            "faf99_afr_XX": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XX samples of African/African-American ancestry."
            },
            "faf95_afr_XY": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XY samples of African/African-American ancestry."
            },
            "faf99_afr_XY": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XY samples of African/African-American ancestry."
            },
            "faf95_amr_XX": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XX samples of Latino ancestry."
            },
            "faf99_amr_XX": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XX samples of Latino ancestry."
            },
            "faf95_amr_XY": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XY samples of Latino ancestry."
            },
            "faf99_amr_XY": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XY samples of Latino ancestry."
            },
            "faf95_eas_XX": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XX samples of East Asian ancestry."
            },
            "faf99_eas_XX": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XX samples of East Asian ancestry."
            },
            "faf95_eas_XY": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XY samples of East Asian ancestry."
            },
            "faf99_eas_XY": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XY samples of East Asian ancestry."
            },
            "faf95_nfe_XX": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XX samples of Non-Finnish European ancestry."
            },
            "faf99_nfe_XX": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XX samples of Non-Finnish European ancestry."
            },
            "faf95_nfe_XY": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XY samples of Non-Finnish European ancestry."
            },
            "faf99_nfe_XY": {
                "Description": "Filtering allele frequency (using Poisson 99% CI) for XY samples of Non-Finnish European ancestry."
            },
            "faf95_sas_XX": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XX samples of South Asian ancestry."
            },
            "faf99_sas_XX": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XX samples of South Asian ancestry."
            },
            "faf95_sas_XY": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XY samples of South Asian ancestry."
            },
            "faf99_sas_XY": {
                "Description": "Filtering allele frequency (using Poisson 95% CI) for XY samples of South Asian ancestry."
            },
            "popmax": {"Description": "Population with maximum AF."},
            "AC_popmax": {
                "Description": "Allele count in the population with maximum allele frequency."
            },
            "AN_popmax": {
                "Description": "Total number of alleles in the population with maximum allele frequency."
            },
            "AF_popmax": {
                "Description": "Maximum allele frequency across populations."
            },
            "nhomalt_popmax": {
                "Description": "Count of homozygous individuals in the population with the maximum allele frequency."
            },
            "age_hist_het_bin_freq": {
                "Description": (
                    "Histogram of ages of heterozygous individuals; bin edges are: 30.0|35.0|40.0|45.0|50.0|55.0|60.0|65.0|70.0|75.0|80.0; total number of individuals of any genotype bin: 0|0|7|51796|66083|76337|90827|121484|93585|2424|0|0."
                )
            },
            "age_hist_het_n_smaller": {
                "Description": (
                    "Count of age values falling below lowest histogram bin edge for heterozygous individuals."
                )
            },
            "age_hist_het_n_larger": {
                "Description": (
                    "Count of age values falling above lowest histogram bin edge for heterozygous individuals."
                )
            },
            "age_hist_hom_bin_freq": {
                "Description": (
                    "Histogram of ages of homozygous alternate individuals; bin edges are: 30.0|35.0|40.0|45.0|50.0|55.0|60.0|65.0|70.0|75.0|80.0; total number of individuals of any genotype bin: 0|0|7|51796|66083|76337|90827|121484|93585|2424|0|0."
                )
            },
            "age_hist_hom_n_smaller": {
                "Description": (
                    "Count of age values falling below lowest histogram bin edge for homozygous alternate individuals."
                )
            },
            "age_hist_hom_n_larger": {
                "Description": (
                    "Count of age values falling above lowest histogram bin edge for homozygous alternate individuals."
                )
            },
            "vep": {
                "Description": (
                    "Consequence annotations from Ensembl VEP. More details about VEP output is described here: "
                    "https://uswest.ensembl.org/info/docs/tools/vep/vep_formats.html#output. VEP was run using the LOFTEE "
                    "plugin and information about the additional LOFTEE annotations can be found here: "
                    "https://github.com/konradjk/loftee. Format:"
                    "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info"
                )
            },
            "vep_version": {"Description": "VEP version."},
            "dbsnp_version": {"Description": "dbSNP version."},
        },
    },
}
