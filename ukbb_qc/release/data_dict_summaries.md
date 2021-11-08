# Sample and variant Quality Control (QC) summary
## Sample QC
### Concordance between arrays and exomes
We examined the concordance of the exome data with the UK Biobank array data lifted over to GRCh38. The minimum proportion concordance of non-reference sample genotypes was <b>0.96</b>.
### Sex inference
We used the [`impute_sex`](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.impute_sex) method in Hail to infer sex using <b>common (allele frequency > 0.1%), non-pseudoautosomal, bi-allelic single nucleotide variants (SNVs)</b> on chromosome X. We then aggregated non-pseudoautosomal chromosome X and chromosome Y coverage for all samples and normalized these values by coverage across chromosome 20. We then checked the distribution of chromosome X and Y ploidies for XX and XY karyotypes to determine each karyotype's ploidy cutoffs. The adjusted ploidy cutoffs helped us add additional granularity to the inferred sex, differentiating between X0, XX, XXX, XY, XXY, XYY, and XXYY karyotypes.
### Interval QC
Due to the Hail package's [sparse](https://hail.is/docs/0.2/experimental/vcf_combiner.html#working-with-sparse-matrix-tables) calling method, this dataset included every observed variant from each sample [gVCF](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format). We decided to apply interval QC to determine which intervals had adequate coverage and were therefore most appropriate for sample QC. Our analysis showed that padding 100 base pairs into the introns and using Broad-specific calling intervals reduced sample call rate significantly. We decided to use 50 base pair padded UK Biobank intervals for sample QC analysis and selected intervals where 85% of all samples had a mean coverage of at least 20X.

### Hard filters
The hard filters applied to the raw data and used in downstream QC were: sex imputation filters, low call rate, and low coverage:
- Sex imputation filters (filtered all karyotypes except XX and XY):
- Low call rate (cutoff <b>0.99</b>)
- Low mean coverage (cutoff <b>20X</b>)
We excluded these hard-filtered samples from subsequent sample QC analyses (relatedness inference, ancestry imputation, and outlier detection) so that these low quality samples would not influence our downstream results.

### Relatedness inference
We used the [`pc_relate`](https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.pc_relate) method in Hail to infer relatedness on variants that are <b>autosomal, bi-allelic single nucleotide variants (SNVs), common (allele frequency > 0.1%), high call rate (> 99%), and LD-pruned with a cutoff of r2 = 0.1</b>. We then used the [`maximal_independent_set`](https://hail.is/docs/0.2/methods/misc.html#hail.methods.maximal_independent_set) method in Hail to get the largest set of samples with no pair of samples related at the 2nd degree or closer, prioritizing samples with greater mean depth. The related samples that are not in the maximal independent set are flagged as related.

### Ancestry inference
We used a hybrid method to infer ancestry. We projected the UK Biobank data onto the [gnomAD population principal components](https://gnomad.broadinstitute.org/downloads#v2-ancestry-classification) (PCs) and then used the gnomAD random forest classifier to assign ancestry to the UK Biobank samples. We observed that many samples being labeled as African using this method were being flagged as outliers by our population stratified outlier detection method. This seemed to be due to the fact that there is a cluster of samples being labeled as African that look highly admixed. To account for this, we ran a principal component analysis (PCA) on the UK Biobank samples and applied a clustering method ([HDBSCAN](https://hdbscan.readthedocs.io/en/latest/index.html)). We found that this clustering method split the African labeled samples into unique clusters and reduced the number of samples flagged as outliers while recapturing most of the same global population clusters observed in gnomAD. As a result, we chose to assign ancestry using a hybrid of the projection onto gnomAD PCs and the UK Biobank specific PCA clustering: for any sample that was assigned to a cluster using the UK Biobank-specific PCA, the sample was given that cluster as their ancestry assignment in order to preserve the sub-structure observed using the clustering approach. Any sample that was not assigned to a cluster was given the label from the initial (gnomAD) PCA projection and random forest classification.

This hybrid ancestry label was used in sample QC to detect outlier samples. However, we switched to using each sample's ancestry label from the Pan-UK Biobank project (https://pan.ukbb.broadinstitute.org/) to calculate genotype callstats for maximum compatibility with downstream analyses that used the Pan-UK Biobank labels.

### Outlier detection
We flagged any sample falling outside <b>4</b> median absolute deviations (MAD) from the median of any of the following metrics (stratified by hybrid ancestry label and tranche as a proxy for platform):
- Number of deletions
- Number of insertions
- Number of SNVs
- Insertion : deletion ratio
- Transition : transversion (TiTv) ratio
- Heterozygous : homozygous ratio

Note that a sample's tranche refers to the data freeze when a sample was released: 150K (the initial 50K exome public release), 100K, 200K, 300K, or 455K.
We ran outlier detection on all samples that passed hard filters and included all related samples.

### Sample QC filter summary
| Filter type | Filtered | Not Filtered | Percent Filtered |
| ----------- | -------- | ------------ | ---------------- |
| Hard filters | 683 | 453988 | 0.15 |
| Relatedness filters | 29899 | 424772 | 6.58 |
| Outlier filters | 3061 | 450927 | 0.67 |
| All filters | 33459 | 421212 | 7.36 | 

The majority of the samples removed from QC were samples inferred to have a second-degree or greater relationship with other samples in the dataset. This finding is not unexpected, as [previous studies](https://doi.org/10.1038/s41586-018-0579-z) have shown that the UK Biobank data contains a large number (approximately 30%) of related samples (Bycroft et al. 2018). 

## Variant Quality Control
For variant QC, we used a combination of a random forest (RF) classifier and hard filters.

### Random forest (RF) training sets
The following training sets were used as true positives for training the random forest model:
- Omni: SNVs present on the Omni 2.5 genotyping array and found in 1000 Genomes data
- Mills Indels present in the Mills and Devine data
- Transmitted singletons: Variants found in two and only two individuals, which were a parent-offspring pair
- Sibling singletons: Variants found in two and only two individuals, which were a sibling pair
- Common (AF > 0.1%) and concordant (> 90% non-reference concordance) array variants
Variants failing traditional GATK hard filters: <b>QD < 2 | FS > 60 | MQ < 30</b> were used as false positives for the random forest training. 

To balance the number of variants used as true positive and false positives in the training, we randomly downsample the false positives training set to match the number of training examples in the true positive set.

### Random forest training features
We used the following allele and site annotations as features in the random forest model:
- Allele type: SNV, indel
- Number of alleles: Total number of alleles present at the site
- Variant type: SNV, indel, multi-SNV, multi-indel, mixed
- Mixed site: Whether more than one allele type is present at the site
- Spanning deletion: Whether there is a spanning deletion (STAR_AC > 0) that overlaps the site
- Quality by depth: Sum of the non-reference genotype likelihoods divided by the sum of the depth in all carriers of that allele
- Read position RankSum: Rank Sum Test for relative positioning of reference versus alternate alleles within reads
- Mapping quality RankSum: Rank Sum Test for mapping qualities of reference versus alternate reads
- Strand bias odds ratio: Symmetric Odds Ratio test of 2x2 contingency table to detect strand bias
- Max probability of allele balance: Highest p-value for sampling the observed allele balance under a binomial model with p=0.5 (maximum across heterozygotes)

### Random forest filtering threshold
In order to set a threshold for the variants that PASS the RF filter, we looked at a few metrics:
- Number of de novo mutations in the 672 trios
- Sensitivity to ClinVar variants
- Precision recall curves in two truth samples present in our data: NA12878 and a pseudo-diploid sample (syndip; https://www.nature.com/articles/s41592-018-0054-7?WT.feed_name=subjects_standards;%20https://github.com/lh3/CHM-eval). Note: syndip was sequenced at Broad, not with the UK Biobank cohort.
- Ti/Tv ratio
- Proportion singletons
- Proportion bi-allelic variants

The thresholds were chosen to maximize sensitivity and specificity based on the metrics listed above. Final thresholds were RF true positive bin of <b>90</b> in well covered intervals for single nucleotide variants and <b>75</b> in well covered intervals for indels.

### Variant hard filters
Two hard filters were also used to exclude variants:
- Excess heterozygotes defined by an inbreeding coefficient <b>-0.3</b>
- Variants where no sample had a high quality genotype <b>(depth >= 10, genotype quality >= 20 and minor allele balance > 0.2 for heterozygous genotypes)</b>