GLOBAL_SAMPLE_ANNOTATIONS = {
    "array_concordance_sites_cutoffs": {
        "Description": (
            "Contains information about the call rate and allele frequency cutoffs used to determine which sites "
            "to include when checking array and exome concordance."
        )
    },
    "sex_imputation_ploidy_cutoffs": {
        "Description": (
            "Contains sex chromosome ploidy cutoffs used when determining sex chromosome karyotypes for the gnomAD sex "
            "imputation. Format: (upper cutoff for single X, (lower cutoff for double X, upper cutoff for double X), "
            "lower cutoff for triple X) and (lower cutoff for single Y, upper cutoff for single Y), lower cutoff for "
            "double Y)."
        )
    },
    "platform_inference_hdbscan_parameters": {
        "Description": (
            "Contains HDBSCAN parameters (minimum cluster size and minimum number of samples) "
            "used in inferring platform labels."
        )
    },
    "population_inference_pca_metrics": {
        "Description": (
            "Contains the number of principal components (PCs) used when running PC-project, the minimum cutoff "
            "probability of belonging to a given population for assignment, and number of PCs used when performing "
            "the from-scratch ancestry principal component analysis (PCA)."
        )
    },
    "relatedness_inference_cutoffs": {
        "Description": (
            "Contains the minor allele frequency cutoff for sites included in running ancestry imputation, "
            "the minimum kinship threshold for emitting a pair of samples in PC relate, the IBD0 cutoff used "
            "to distinguish parent-child relationships from full sibling relationships, the minimum kinship "
            "threshold for determining second-degree relatedness, and the kinship thresholds for determining "
            "whether a pair had first-degree relatedness."
        )
    },
    "outlier_detection_metrics": {
        "Description": (
            "Contains the medians, the median absolute deviations (MADs), the upper MAD threshold, and the "
            "lower MAD threshold for each sample QC metric used in outlier filtering. The metrics are: number "
            "of SNPs (n_SNP), transition to transversion ratio (r_ti_tv), insertion to deletion ratio "
            "(r_insertion_deletion), number of insertions (n_insertion), number of deletions (n_deletion), "
            "and ratio of heterozygous to homozygous variants (r_het_hom_var)."
        )
    },
}

SAMPLE_ANNOTATIONS = {
    "s": {"Description": "Sample ID."},
    "ukbb_meta": {
        "Description": "Struct containing metadata information provided by the UK Biobank (UKBB).",
        "sub_annotations": {
            "batch": {"Description": "Sample's data tranche (150K, 100K, 200K, 300K)."},
            "batch_num": {
                "Description": (
                    "Numeric version of batch. 0: 150K, 1: 100K, 2: 200K, 3: 300K."
                )
            },
            "ukbb_app_26041_id": {"Description": "UKBB sample ID."},
            "withdrawn_consent": {
                "Description": ("Whether sample has withdrawn consent.")
            },
            "age": {"Description": ("Sample's age at recruitment.")},
            "array_sex": {
                "Description": ("Sample sex as inferred from genotype data.")
            },
            "self_reported_ancestry": {
                "Description": "Sample's self reported ancestry."
            },
        },
    },
    "array_concordance": {
        "Description": "Struct containing sample concordance between array and exome data.",
        "sub_annotations": {
            "concordance": {
                "Description": (
                    "List of lists containing concordance as calculated using `hl.concordance` for: "
                    "No Data (missing variant), No Call (missing genotype call), Homozygous Reference "
                    "Heterozygous, Homozygous Variant"
                )
            },
            "n_discordant": {
                "Description": (
                    "Number of genotypes which were called (homozygous reference, heterozygous, or homozygous "
                    "variant) in both datasets, but where the call did not match between the two. Provided by "
                    "`hl.concordance`."
                )
            },
            "num_gt_con_non_ref": {
                "Description": "Number of non-reference genotypes that were concordant between the exome and array datasets."
            },
            "num_gt_non_ref": {
                "Description": (
                    "Total number of non-ref genotypes between the array and exome data."
                )
            },
            "prop_gt_con_non_ref": {
                "Description": (
                    "Total proportion of concordant non-reference genotypes. "
                )
            },
        },
    },
    "sex_imputation": {
        "Description": "Struct containing sex imputation information.",
        "sub_annotations": {
            "is_female": {
                "Description": "True if the imputed sex (from `hl.impute_sex`) is female, false if male, missing if undetermined."
            },
            "f_stat": {
                "Description": "Inbreeding coefficient (excess heterozygosity) on chromosome X."
            },
            "n_called": {"Description": "Number of variants with a genotype call."},
            "expected_homs": {"Description": "Expected number of homozygotes."},
            "observed_homs": {"Description": "Observed number of homozygotes."},
            "chr20_mean_dp": {
                "Description": "Sample's mean depth across chromosome 20."
            },
            "chrX_mean_dp": {"Description": "Sample's mean depth across chromosome X."},
            "chrY_mean_dp": {"Description": "Sample's mean depth across chromosome Y."},
            "chrX_ploidy": {
                "Description": "Sample's chromosome X ploidy (chrX_mean_dp normalized using chr20_mean_dp)."
            },
            "chrY_ploidy": {
                "Description": "Sample's chromosome Y ploidy (chrY_mean_dp normalized using chr20_mean_dp)."
            },
            "X_karyotype": {"Description": "Sample's chromosome X karyotype."},
            "Y_karyotype": {"Description": "Sample's chromosome Y karyotype."},
            "sex_karyotype": {
                "Description": "Sample's sex karyotype (combined X and Y karyotype)."
            },
        },
    },
    "sample_qc": {
        "Description": "Struct containing sample QC metrics calculated using hl.sample_qc().",
        "sub_annotations": {
            "call_rate": {"Description": "Fraction of calls not missing or filtered."},
            "dp_mean": {"Description": "Sample mean depth."},
            "dp_stdev": {"Description": "Standard deviation for sample depth."},
        },
    },
    "platform_inference": {
        "Description": "Struct containing platform inference information.",
        "sub_annotations": {
            "callrate_pcs": {"Description": "Sample's scores for each call rate PC."},
            "qc_platform": {"Description": "Sample's inferred platform."},
        },
    },
    "hybrid_pop_data": {
        "Description": "Struct containing information from hybrid ancestry assignment.",
        "sub_annotations": {
            "scores": {
                "Description": "Sample scores for each PC used in from-scratch ancestry PCA."
            },
            "cluster": {
                "Description": "Sample's cluster assignment from the from-scratch ancestry PCA."
            },
            "pop": {"Description": "Sample's inferred hybrid ancestry assignment."},
        },
    },
    "sample_filters": {
        "Description": "Struct containing Booleans for any sample flags or filters.",
        "sub_annotations": {
            "low_callrate": {
                "Description": "Whether the sample was removed from QC for having a call rate less than 0.99."
            },
            "ambiguous_sex": {
                "Description": "Whether the sample was removed from QC for having ambiguous sex."
            },
            "sex_aneuploidy": {
                "Description": "Whether the sample was removed from QC for having a sex aneuploidy."
            },
            "related": {
                "Description": (
                    "Whether the sample was removed from QC when filtering related samples. True if sample had a second degree "
                    "relationship with a second sample."
                )
            },
            "duplicate": {
                "Description": (
                    "Whether the sample was removed from QC when filtering related samples. True if sample had a duplicate or "
                    "monozygotic twin relationship with a second sample."
                )
            },
            "parent_child": {
                "Description": (
                    "Whether the sample was removed from QC when filtering related samples. True if sample had a parent-child "
                    "relationship with a second sample."
                )
            },
            "sibling": {
                "Description": (
                    "Whether the sample was removed from QC when filtering related samples. True if sample had a full sibling "
                    "relationship with a second sample."
                )
            },
            "control": {
                "Description": "True if sample is a control sample (NA12878, syndip)."
            },
            "fail_n_snp": {
                "Description": "Whether the sample was removed from QC for having a n_snp value larger upper MAD or smaller than the lower MAD threshold."
            },
            "fail_r_ti_tv": {
                "Description": "Whether the sample was removed from QC for having a r_ti_tv value larger upper MAD or smaller than the lower MAD threshold."
            },
            "fail_r_insertion_deletion": {
                "Description": "Whether the sample was removed from QC for having a r_insertion_deletion value larger upper MAD or smaller than the lower MAD threshold."
            },
            "fail_n_insertion": {
                "Description": "Whether the sample was removed from QC for having a n_insertion value larger upper MAD or smaller than the lower MAD threshold."
            },
            "fail_n_deletion": {
                "Description": "Whether the sample was removed from QC for having a n_deletion value larger upper MAD or smaller than the lower MAD threshold."
            },
            "fail_r_het_hom_var": {
                "Description": "Whether the sample was removed from QC for having a r_het_hom_var value larger upper MAD or smaller than the lower MAD threshold."
            },
            "high_quality": {
                "Description": "Whether the sample is high quality (True if sample was not hard filtered or removed during outlier detection)."
            },
            "release": {
                "Description": "Whether the sample should be included in the final release (true if the sample has a defined `batch` and is high quality)."
            },
            "qc_metrics_filters": {
                "Description": "Set summarizing any outlier detection flags for this sample."
            },
        },
    },
    "relatedness_inference": {
        "Description": "Struct containing relatedness inference information.",
        "sub_annotations": {
            "relationships": {
                "Description": (
                    "Set containing all relationship strings (e.g., parent-child) for this sample."
                )
            },
        },
    },
    "gnomad_pc_project_pop_data": {
        "Description": (
            "Struct containing ancestry information derived from projecting the UKBB data onto the gnomAD population PCs and using a random forest classifier trained on gnomAD ancestry labels to assign ancestry."
        ),
        "sub_annotations": {
            "scores": {"Description": "Sample's scores for each gnomAD population PC."},
            "pop": {
                "Description": "Sample's inferred gnomAD PC project population label."
            },
        },
    },
}


TEMP = {
    "gnomad_age_distribution": {
        "Description": "gnomAD callset-wide age histogram calculated on release samples.",
        "sub_globals": {
            "bin_edges": {"Description": "Bin edges for the age histogram."},
            "bin_freq": {
                "Description": "Bin frequencies for the age histogram. This is the number of records found in each bin."
            },
            "n_smaller": {
                "Description": "Count of age values falling below lowest histogram bin edge."
            },
            "n_larger": {
                "Description": "Count of age values falling above highest histogram bin edge."
            },
        },
    },
}
GLOBAL_VARIANT_ANNOTATIONS = {
    "hgdp_tgp_freq_meta": {
        "Description": (
            "HGDP and 1KG frequency metadata. An ordered list containing the frequency aggregation group for each "
            "element of the hgdp_tgp_freq array row annotation."
        )
    },
    "gnomad_freq_meta": {
        "Description": (
            "gnomAD frequency metadata. An ordered list containing the frequency aggregation group for each element of "
            "the gnomad_freq array row annotation."
        )
    },
    "hgdp_tgp_freq_index_dict": {
        "Description": (
            "Dictionary keyed by specified label grouping combinations (group: adj/raw, pop: HGDP or 1KG subpopulation, "
            "sex: sex karyotype), with values describing the corresponding index of each grouping entry in the "
            "HGDP + 1KG frequency array annotation."
        )
    },
    "gnomad_freq_index_dict": {
        "Description": (
            "Dictionary keyed by specified label grouping combinations (group: adj/raw, pop: gnomAD inferred global "
            "population sex: sex karyotype), with values describing the corresponding index of each grouping entry in "
            "the gnomAD frequency array annotation."
        )
    },
    "gnomad_faf_meta": {
        "Description": (
            "gnomAD filtering allele frequency metadata. An ordered list containing the frequency aggregation group "
            "for each element of the gnomad_faf array row annotation."
        )
    },
    "gnomad_faf_index_dict": {
        "Description": (
            "Dictionary keyed by specified label grouping combinations (group: adj/raw, pop: gnomAD inferred global "
            "population sex: sex karyotype), with values describing the corresponding index of each grouping entry in "
            "the filtering allele frequency (using Poisson 99% CI) annotation."
        )
    },
    "variant_filtering_model": {
        "Description": {"The variant filtering model used and its specific cutoffs."},
        "sub_globals": {
            "model_name": {
                "Description": (
                    "Variant filtering model name used in the 'filters' row annotation to indicate the variant was "
                    "filtered by the model during variant QC."
                )
            },
            "score_name": {"Description": "Name of score used for variant filtering."},
            "snv_cutoff": {
                "Description": "SNV filtering cutoff information.",
                "sub_globals": {
                    "bin": {"Description": "Filtering percentile cutoff for SNVs."},
                    "min_score": {
                        "Description": "Minimum score at SNV filtering percentile cutoff."
                    },
                },
            },
            "indel_cutoff": {
                "Description": "Information about cutoff used for indel filtering.",
                "sub_globals": {
                    "bin": {"Description": "Filtering percentile cutoff for indels."},
                    "min_score": {
                        "Description": "Minimum score at indel filtering percentile cutoff."
                    },
                },
            },
            "snv_training_variables": {
                "Description": "Variant annotations used as features in SNV filtering model."
            },
            "indel_training_variables": {
                "Description": "Variant annotations used as features in indel filtering model."
            },
        },
    },
    "variant_inbreeding_coeff_cutoff": {
        "Description": "Hard-filter cutoff for InbreedingCoeff on variants."
    },
    "vep_version": {"Description": "VEP version."},
    "vep_csq_header": {"Description": "VEP header for VCF export."},
    "dbsnp_version": {"Description": "dbSNP version."},
}
GLOBAL_ANNOTATIONS = {**GLOBAL_SAMPLE_ANNOTATIONS, **GLOBAL_VARIANT_ANNOTATIONS}


VARIANT_ANNOTATIONS = {
    "locus": {
        "Description": "Variant locus. Contains contig and position information.",
    },
    "alleles": {"Description": "Variant alleles.",},
    "rsid": {"Description": "dbSNP reference SNP identification (rsID) numbers.",},
    "a_index": {
        "Description": (
            "The original index of this alternate allele in the multiallelic representation (1 is the first alternate "
            "allele or the only alternate allele in a biallelic variant)."
        )
    },
    "was_split": {
        "Description": "True if this variant was originally multiallelic, otherwise False."
    },
    "hgdp_tgp_freq": {
        "Description": (
            "Allele frequency information (AC, AN, AF, homozygote count) in HGDP + 1KG samples that pass the "
            "high_quality sample annotation and are inferred as unrelated (False in relatedness_inference.related "
            "annotation)."
        ),
        "sub_annotations": {
            "AC": {
                "Description": "Alternate allele count  in HGDP + 1KG samples that pass the high_quality sample annotation."
            },
            "AF": {
                "Description": "Alternate allele frequency  in HGDP + 1KG samples that pass the high_quality sample annotation."
            },
            "AN": {
                "Description": "Total number of alleles in HGDP + 1KG samples that pass the high_quality sample annotation."
            },
            "homozygote_count": {
                "Description": "Count of homozygous individuals in HGDP + 1KG samples that pass the high_quality sample annotation."
            },
        },
    },
    "gnomad_freq": {
        "Description": "Allele frequency information (AC, AN, AF, homozygote count) in gnomAD release.",
        "sub_annotations": {
            "AC": {"Description": "Alternate allele count in gnomAD release."},
            "AF": {"Description": "Alternate allele frequency in gnomAD release."},
            "AN": {"Description": "Total number of alleles in gnomAD release."},
            "homozygote_count": {
                "Description": "Count of homozygous individuals in gnomAD release."
            },
        },
    },
    "gnomad_popmax": {
        "Description": "Allele frequency information (AC, AN, AF, homozygote count) for the population with maximum AF in gnomAD.",
        "sub_annotations": {
            "AC": {
                "Description": "Allele count in the population with the maximum AF in gnomAD."
            },
            "AF": {
                "Description": "Maximum allele frequency across populations in gnomAD."
            },
            "AN": {
                "Description": "Total number of alleles in the population with the maximum AF in gnomAD."
            },
            "homozygote_count": {
                "Description": "Count of homozygous individuals in the population with the maximum allele frequency in gnomAD."
            },
            "pop": {"Description": "Population with maximum AF in gnomAD."},
            "faf95": {
                "Description": (
                    "Filtering allele frequency (using Poisson 95% CI) for the population with the maximum allele "
                    "frequency in gnomAD."
                )
            },
        },
    },
    "gnomad_faf": {
        "Description": "Filtering allele frequency in gnomAD release.",
        "sub_annotations": {
            "faf95": {
                "Description": "Filtering allele frequency in gnomAD release (using Poisson 95% CI)."
            },
            "faf99": {
                "Description": "Filtering allele frequency in gnomAD release (using Poisson 99% CI)."
            },
        },
    },
    "gnomad_qual_hists": {
        "Description": "gnomAD genotype quality metric histograms for high quality genotypes.",
        "sub_annotations": {
            "gq_hist_all": {
                "Description": "Histogram for GQ calculated on high quality genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the GQ histogram calculated on high quality genotypes are: "
                            "0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100"
                        ),
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the GQ histogram calculated on high quality genotypes. The number of "
                            "records found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": (
                            "Count of GQ values falling below lowest histogram bin edge, for GQ calculated on high "
                            "quality genotypes"
                        ),
                    },
                    "n_larger": {
                        "Description": (
                            "Count of GQ values falling above highest histogram bin edge, for GQ calculated on high "
                            "quality genotypes"
                        ),
                    },
                },
            },
            "dp_hist_all": {
                "Description": "Histogram for DP calculated on high quality genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the DP histogram calculated on high quality genotypes are: "
                            "0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                        )
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the DP histogram calculated on high quality genotypes. The number of "
                            "records found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": (
                            "Count of DP values falling below lowest histogram bin edge, for DP calculated on high "
                            "quality genotypes."
                        )
                    },
                    "n_larger": {
                        "Description": (
                            "Count of DP values falling above highest histogram bin edge, for DP calculated on high "
                            "quality genotypes."
                        )
                    },
                },
            },
            "gq_hist_alt": {
                "Description": "Histogram for GQ in heterozygous individuals calculated on high quality genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the histogram of GQ in heterozygous individuals calculated on high quality "
                            "genotypes are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                        ),
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the histogram of GQ in heterozygous individuals calculated on high "
                            "quality genotypes. The number of records found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": (
                            "Count of GQ values falling below lowest histogram bin edge, for GQ in heterozygous "
                            "individuals calculated on high quality genotypes."
                        )
                    },
                    "n_larger": {
                        "Description": (
                            "Count of GQ values falling above highest histogram bin edge, for GQ in heterozygous "
                            "individuals calculated on high quality genotypes."
                        )
                    },
                },
            },
            "dp_hist_alt": {
                "Description": "Histogram for DP in heterozygous individuals calculated on high quality genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the histogram of DP in heterozygous individuals calculated on high quality "
                            "genotypes are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                        )
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the histogram of DP in heterozygous individuals calculated on high "
                            "quality genotypes. The number of records found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": (
                            "Count of DP values falling below lowest histogram bin edge, for DP in heterozygous "
                            "individuals calculated on high quality genotypes."
                        )
                    },
                    "n_larger": {
                        "Description": (
                            "Count of DP values falling above highest histogram bin edge, for DP in heterozygous "
                            "individuals calculated on high quality genotypes."
                        )
                    },
                },
            },
            "ab_hist_alt": {
                "Description": "Histogram for AB in heterozygous individuals calculated on high quality genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the histogram of AB in heterozygous individuals calculated on high quality "
                            "genotypes are: 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00."
                        )
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the histogram of AB in heterozygous individuals calculated on high "
                            "quality genotypes. The number of records found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": (
                            "Count of AB values falling below lowest histogram bin edge, for AB in heterozygous "
                            "individuals calculated on high quality genotypes."
                        )
                    },
                    "n_larger": {
                        "Description": (
                            "Count of AB values falling above highest histogram bin edge, for AB in heterozygous "
                            "individuals calculated on high quality genotypes."
                        )
                    },
                },
            },
        },
    },
    "gnomad_raw_qual_hists": {
        "Description": "gnomAD genotype quality metric histograms.",
        "sub_annotations": {
            "gq_hist_all": {
                "Description": "Histogram for GQ calculated on all genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the GQ histogram calculated on all genotypes are: "
                            "0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                        )
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the GQ histogram calculated on all genotypes. The number of records "
                            "found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": (
                            "Count of GQ values falling below lowest histogram bin edge, for GQ calculated on all "
                            "genotypes."
                        )
                    },
                    "n_larger": {
                        "Description": (
                            "Count of GQ values falling above highest histogram bin edge, for GQ calculated on all "
                            "genotypes."
                        )
                    },
                },
            },
            "dp_hist_all": {
                "Description": "Histogram for DP calculated on all genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the DP histogram calculated on all genotypes are: "
                            "0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100"
                        )
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the DP histogram calculated on all genotypes. The number of records "
                            "found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": (
                            "Count of DP values falling below lowest histogram bin edge, for DP calculated on all "
                            "genotypes."
                        )
                    },
                    "n_larger": {
                        "Description": (
                            "Count of DP values falling above highest histogram bin edge, for DP calculated on all "
                            "genotypes."
                        )
                    },
                },
            },
            "gq_hist_alt": {
                "Description": "Histogram for GQ in heterozygous individuals calculated on all genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the histogram of GQ in heterozygous individuals calculated on all genotypes "
                            "are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                        )
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the histogram of GQ in heterozygous individuals calculated on all "
                            "genotypes. The number of records found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": (
                            "Count of GQ values falling below lowest histogram bin edge, for GQ in heterozygous "
                            "individuals calculated on all genotypes."
                        )
                    },
                    "n_larger": {
                        "Description": (
                            "Count of GQ values falling above highest histogram bin edge, for GQ in heterozygous "
                            "individuals calculated on all genotypes."
                        )
                    },
                },
            },
            "dp_hist_alt": {
                "Description": "Histogram for DP in heterozygous individuals calculated on all genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the histogram of DP in heterozygous individuals calculated on all genotypes "
                            "are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                        )
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the histogram of DP in heterozygous individuals calculated on all "
                            "genotypes. The number of records found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": "Count of DP values falling below lowest histogram bin edge, for DP in heterozygous individuals calculated on all genotypes."
                    },
                    "n_larger": {
                        "Description": "Count of DP values falling above highest histogram bin edge, for DP in heterozygous individuals calculated on all genotypes."
                    },
                },
            },
            "ab_hist_alt": {
                "Description": "Histogram for AB in heterozygous individuals calculated on all genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the histogram of AB in heterozygous individuals calculated on all genotypes "
                            "are: 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00."
                        )
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the histogram of AB in heterozygous individuals calculated on all "
                            "genotypes. The number of records found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": (
                            "Count of AB values falling below lowest histogram bin edge, for AB in heterozygous "
                            "individuals calculated on all genotypes."
                        )
                    },
                    "n_larger": {
                        "Description": (
                            "Count of AB values falling above highest histogram bin edge, for AB in heterozygous "
                            "individuals calculated on all genotypes."
                        )
                    },
                },
            },
        },
    },
    "gnomad_age_hist_het": {
        "Description": "Histogram for age in all heterozygous gnomAD release samples calculated on high quality genotypes.",
        "sub_annotations": {
            "bin_edges": {"Description": "Bin edges for the age histogram."},
            "bin_freq": {
                "Description": "Bin frequencies for the age histogram. This is the number of records found in each bin."
            },
            "n_smaller": {
                "Description": "Count of age values falling below lowest histogram bin edge."
            },
            "n_larger": {
                "Description": "Count of age values falling above highest histogram bin edge."
            },
        },
    },
    "gnomad_age_hist_hom": {
        "Description": "Histogram for age in all homozygous gnomAD release samples calculated on high quality genotypes.",
        "sub_annotations": {
            "bin_edges": {"Description": "Bin edges for the age histogram."},
            "bin_freq": {
                "Description": "Bin frequencies for the age histogram. This is the number of records found in each bin."
            },
            "n_smaller": {
                "Description": "Count of age values falling below lowest histogram bin edge."
            },
            "n_larger": {
                "Description": "Count of age values falling above highest histogram bin edge."
            },
        },
    },
    "filters": {
        "Description": (
            "Variant filters; AC0: Allele count is zero after filtering out low-confidence genotypes (GQ < 20; DP < 10;"
            " and AB < 0.2 for het calls), AS_VQSR: Failed VQSR filtering thresholds of -2.7739 for SNPs and -1.0606 "
            "for indels, InbreedingCoeff: GATK InbreedingCoeff < -0.3, PASS: Passed all variant filters."
        )
    },
    "info": {
        "Description": "Struct containing typical GATK allele-specific (AS) info fields and additional variant QC fields.",
        "sub_annotations": {
            "QUALapprox": {
                "Description": "Sum of PL[0] values; used to approximate the QUAL score."
            },
            "SB": {
                "Description": (
                    "Per-sample component statistics which comprise the Fisher's exact test to detect strand bias. "
                    "Values are: depth of reference allele on forward strand, depth of reference allele on reverse "
                    "strand, depth of alternate allele on forward strand, depth of alternate allele on reverse strand."
                )
            },
            "MQ": {
                "Description": "Root mean square of the mapping quality of reads across all samples."
            },
            "MQRankSum": {
                "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities."
            },
            "VarDP": {
                "Description": "Depth over variant genotypes (does not include depth of reference samples)."
            },
            "AS_ReadPosRankSum": {
                "Description": "Allele-specific z-score from Wilcoxon rank sum test of alternate vs. reference read position bias."
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
            "AS_MQ": {
                "Description": "Allele-specific root mean square of the mapping quality of reads across all samples."
            },
            "QD": {
                "Description": "Variant call confidence normalized by depth of sample reads supporting a variant."
            },
            "AS_MQRankSum": {
                "Description": "Allele-specific z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities."
            },
            "FS": {
                "Description": "Phred-scaled p-value of Fisher's exact test for strand bias."
            },
            "AS_FS": {
                "Description": "Allele-specific phred-scaled p-value of Fisher's exact test for strand bias."
            },
            "ReadPosRankSum": {
                "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference read position bias."
            },
            "AS_QUALapprox": {
                "Description": "Allele-specific sum of PL[0] values; used to approximate the QUAL score."
            },
            "AS_SB_TABLE": {
                "Description": "Allele-specific forward/reverse read counts for strand bias tests."
            },
            "AS_VarDP": {
                "Description": "Allele-specific depth over variant genotypes (does not include depth of reference samples)."
            },
            "AS_SOR": {
                "Description": "Allele-specific strand bias estimated by the symmetric odds ratio test."
            },
            "SOR": {
                "Description": "Strand bias estimated by the symmetric odds ratio test."
            },
            "transmitted_singleton": {
                "Description": (
                    "Variant was a callset-wide doubleton that was transmitted within a family from a parent to a "
                    "child (i.e., a singleton amongst unrelated samples in cohort)."
                )
            },
            "omni": {
                "Description": "Variant is present on the Omni 2.5 genotyping array and found in 1000 Genomes data."
            },
            "mills": {"Description": "Indel is present in the Mills and Devine data."},
            "monoallelic": {
                "Description": "All samples are all homozygous alternate for the variant."
            },
            "InbreedingCoeff": {
                "Description": (
                    "Inbreeding coefficient, the excess heterozygosity at a variant site, computed as 1 - (the number "
                    "of heterozygous genotypes)/(the number of heterozygous genotypes expected under Hardy-Weinberg "
                    "equilibrium)."
                )
            },
        },
    },
    "vep": {
        "Description": (
            "Consequence annotations from Ensembl VEP. More details about VEP output is described here: "
            "https://uswest.ensembl.org/info/docs/tools/vep/vep_formats.html#output. VEP was run using the LOFTEE "
            "plugin and information about the additional LOFTEE annotations can be found here: "
            "https://github.com/konradjk/loftee."
        )
    },
    "vqsr": {
        "Description": "VQSR related variant annotations.",
        "sub_annotations": {
            "AS_VQSLOD": {
                "Description": (
                    "Allele-specific log-odds ratio of being a true variant versus being a false positive under the "
                    "trained VQSR Gaussian mixture model."
                )
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
        },
    },
    "region_flag": {
        "Description": "Struct containing flags for problematic regions.",
        "sub_annotations": {
            "lcr": {"Description": "Variant falls within a low complexity region."},
            "segdup": {
                "Description": "Variant falls within a segmental duplication region."
            },
        },
    },
    "allele_info": {
        "Description": "Allele information.",
        "sub_annotations": {
            "variant_type": {
                "Description": "Variant type (snv, indel, multi-snv, multi-indel, or mixed).",
            },
            "allele_type": {
                "Description": "Allele type (snv, insertion, deletion, or mixed).",
            },
            "n_alt_alleles": {
                "Description": "Total number of alternate alleles observed at variant locus.",
            },
        },
    },
    "was_mixed": {"Description": "Variant type was mixed."},
    "cadd": {
        "sub_annotations": {
            "raw_score": {
                "Description": (
                    "Raw CADD scores are interpretable as the extent to which the annotation profile for a given "
                    "variant suggests that the variant is likely to be 'observed' (negative values) vs 'simulated' "
                    "(positive values); higher values indicate that a variant is more likely to be simulated (or 'not "
                    "observed') and therefore more likely to have deleterious effects. More information can be found "
                    "on the CADD website: https://cadd.gs.washington.edu/info."
                )
            },
            "phred": {
                "Description": (
                    "CADD Phred-like scores ('scaled C-scores') ranging from 1 to 99, based on the rank of each "
                    "variant relative to all possible 8.6 billion substitutions in the human reference genome. Larger "
                    "values are more deleterious. More information can be found on the CADD website: "
                    "https://cadd.gs.washington.edu/info."
                )
            },
            "has_duplicate": {
                "Description": (
                    "A True/False flag that indicates whether the variant has more than one CADD score associated with "
                    "it. For a small set of variants, the in silico predictors calculated multiple scores per variant "
                    "based on additional information. For example, if a variant is found in multiple transcripts or if "
                    "it has multiple trinucleotide contexts, an in silico predictor may report scores for multiple "
                    "scenarios. The highest score was taken for each variant in the cases where the in silico predictor "
                    "calculates multiple scores, and we flag variants with multiple scores."
                )
            },
        },
    },
    "revel": {
        "Description": (
            "dbNSFP's Revel score, ranging from 0 to 1. Variants with higher scores are predicted to be more likely to "
            "be deleterious."
        ),
        "sub_annotations": {
            "revel_score": {"Description": "Revel’s numerical score from 0 to 1."},
            "has_duplicate": {
                "Description": (
                    "A True/False flag that indicates whether the variant has more than one revel_score associated "
                    "with it. For a small set of variants, the in silico predictors calculated multiple scores per "
                    "variant based on additional information. For example, if a variant is found in multiple "
                    "transcripts or if it has multiple trinucleotide contexts, an in silico predictor may report "
                    "scores for multiple scenarios. The highest score was taken for each variant in the cases where "
                    "the in silico predictor calculates multiple scores, and we flag variants with multiple scores."
                )
            },
        },
    },
    "splice_ai": {
        "sub_annotations": {
            "splice_ai": {
                "Description": "The maximum delta score, interpreted as the probability of the variant being splice-altering."
            },
            "splice_consequence": {
                "Description": "The consequence term associated with the max delta score in 'splice_ai’."
            },
            "has_duplicate": {
                "Description": (
                    "A True/False flag that indicates whether the variant has more than one splice_ai score associated "
                    "with it. For a small set of variants, the in silico predictors calculated multiple scores per "
                    "variant based on additional information. For example, if a variant is found in multiple "
                    "transcripts or if it has multiple trinucleotide contexts, an in silico predictor may report "
                    "scores for multiple scenarios. The highest score was taken for each variant in the cases where "
                    "the in silico predictor calculates multiple scores, and we flag variants with multiple scores."
                )
            },
        },
    },
    "primate_ai": {
        "sub_annotations": {
            "primate_ai_score": {
                "Description": "PrimateAI's deleteriousness score from 0 (less deleterious) to 1 (more deleterious)."
            },
            "has_duplicate": {
                "Description": (
                    "A True/False flag that indicates whether the variant has more than one primate_ai_score associated "
                    "with it. For a small set of variants, the in silico predictors calculated multiple scores per "
                    "variant based on additional information. For example, if a variant is found in multiple "
                    "transcripts or if it has multiple trinucleotide contexts, an in silico predictor may report "
                    "scores for multiple scenarios. The highest score was taken for each variant in the cases where "
                    "the in silico predictor calculates multiple scores, and we flag variants with multiple scores."
                )
            },
        },
    },
    "AS_lowqual": {
        "Description": (
            "Whether the variant falls below a low quality threshold and was excluded from the gnomAD dataset. We "
            "recommend filtering all such variants. This is similar to the GATK LowQual filter, but is allele-specific. "
            "GATK computes this annotation at the site level, which uses the least stringent prior for mixed sites."
        )
    },
    "telomere_or_centromere": {
        "Description": (
            "Whether the variant falls within a telomere or centromere region. These variants were excluded from the "
            "gnomAD dataset. We recommend filtering all such variants."
        )
    },
}
