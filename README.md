# Homalt branch
This branch was created to release a patch for the 455K UKBB release. The fix we implemented to correct a technical artifact creating the depletion of homozygous alternate genotypes created situations where samples that were individually heterozygous non-reference had both a heterozygous and homozygous alternate call at the same locus.

For example: 
In the raw MT, chr22:19233125 has the alleles G (reference) and A, C, GTGTGGT, T (alternate). UKB_4046155_339844576's original genotype call at this locus is A/T.
However, the G>A variant call in this sample matches the criteria for the homalt hotfix (AF > 1% and AB > 0.9). Thus, in the release, UKB_4046155_339844576 has a heterozygous call for the G>T variant and a homozygous alternate call for the G>A variant.

See this slack thread for more details: https://atgu.slack.com/archives/CNNTF8Z46/p1620401614102100.
