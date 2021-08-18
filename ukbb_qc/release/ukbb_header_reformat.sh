#!/bin/bash

set -e -o pipefail

if [[ "$1" == "" ]]
then
  echo "Usage: $0 file.vcf.gz"
  echo ""
  echo "Repackages the VCF so that the samples in the header are contained in a separate, 0-compressed bgzip block." echo "Requires bcftools, bgzip, tabix."
  echo "Outputs file.repackaged.vcf.gz+tbi (to disk), and header coordinates (lo,hi) to stdout."
  exit 1
fi

vcfgz_path="$1"
vcfgz_prefix="${vcfgz_path%.vcf.gz}"
outfile="${vcfgz_prefix.repackaged.vcf.gz}"

# Get header up until the sample names, bgzip regularly
( bcftools view --no-version -h "$vcfgz_path" | head -n -1; bcftools view --no-version -h "$vcfgz_path" | tail -n 1 | awk '{for (i=1; i<=9; i++) {printf("%s\t", $i) }}' ) | bgzip | head -c -28 > "$outfile"

# Mark the starting byte range 
lo=$(stat -c '%s' "$outfile")
let lo++

# Reprint sample ids eith 7 digits and bgzip with 0-compression
bcftools view --no-version -h "$vcfgz_path" | tail -n 1 | awk '{for (i=10; i<=NF; i++) {printf "%07d%c", $i, (i==NF?"\n":"\t")}}' | bgzip -l 0 | head -c -28 >> "$outfile"

# Mark the ending byte range
hi=$(stat -c '%s' "$outfile")

# Bgzip the rest of the file
bcftools view -H "$vcfgz_path" | bgzip -@ `nproc` >> "$outfile"

# Tabix the output file 
tabix -p vcf "$outfile"

# Output the byte range of the sample names: 
echo $lo $hi

