GLOBAL_SAMPLE_ANNOTATIONS = {
    "array_concordance_sites_cutoffs": {
        "Description": (
            "Contains information about the call rate and allele frequency cutoffs used to determine which sites "
            "to include when checking array and exome concordance."
        )
    },
    "sex_imputation_ploidy_cutoffs": {
        "Description": (
            "Contains sex chromosome ploidy cutoffs used when determining sex chromosome karyotypes for sex "
            "imputation. Format: (upper cutoff for single X, (lower cutoff for double X, upper cutoff for double X), "
            "lower cutoff for triple X) and (lower cutoff for single Y, upper cutoff for single Y), lower cutoff for "
            "double Y)."
        )
    },
    "platform_inference_hdbscan_parameters": {
        "Description": (
            "Contains HDBSCAN parameters (minimum cluster size and minimum number of samples) "
            "used to infer platform labels."
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
    "s": {"Description": "UKBB application 26041 sample ID."},
    "ukbb_meta": {
        "Description": "Struct containing metadata information provided by the UK Biobank (UKBB).",
        "sub_annotations": {
            "batch": {"Description": "Sample's data tranche (150K, 100K, 200K, 300K)."},
            "batch_num": {
                "Description": (
                    "Numeric version of batch. 0: 150K, 1: 100K, 2: 200K, 3: 300K."
                )
            },
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
                "Description": "Whether the sample should be included in the final release (True if the sample has a defined `batch` and is high quality)."
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
