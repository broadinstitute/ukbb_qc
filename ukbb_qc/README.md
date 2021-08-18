# 300K Return of results branch
This branch was created to re-export the 300K UKBB VCFs in order to return data to the UKBB.

The VCF re-export was necessary to fix a bug with the 300K frequencies. In the original 300K frequency code, we accidentally appended missing "cohort" frequency structs (frequencies calculated on all NFE samples, including relateds) to our "release" frequency struct (frequencies calculated on high quality, unrelated samples only), which resulted in empty/missing frequency structs for impacted loci.

Other fixes in this branch include:
- Code to remove samples that withdrew consent after we completed 300K QC (in `prepare_vcf_release.py`).
- Code to remove duplicate sample IDs from the VCF (in `prepare_vcf_release.py`).
- Code to remove frequencies stratified by subpopulation (to remove hybrid population labels from VCF export; in `prepare_vcf_release.py` and `utils.py`).
- Code to flag segdup regions (in `prepare_vcf_release.py`).
- Code to update rsids to allow for multiple rsids per locus (in `prepare_vcf_release.py`).
- Script to reformat VCF header to allow UKBB sample ID rebadging (`ukbb_header_reformat.sh`).
- Script to get first and last positions in VCF shards (`get_shard_positions.py`).
