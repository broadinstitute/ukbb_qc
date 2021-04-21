#!/bin/bash
usage() {
cat << EOF
    This script extracts variants from gnomAD v2.1 (exomes) and v3.1 (genomes) that have fewer than the max number of samples available for readviz.

    This script outputs TSVs which will be imported into hail Tables in a downstream step.
    
    Script assumes that all sqlite databases are on the same machine that script is being run (i.e., local computer or cloud VM).

    Inputs:
        -n    Maximum number of samples displayed for readviz. Default is 3.
        -d    Data type (should be one of 'exomes' or 'genomes').
        -p    Absolute path to directory that contains sqlite databases.
        -t    Absolute path to directory to store output (temporary) TSV files.
EOF
}

# Check num of command line args
if [[ $# -lt 8 ]]; then
    usage
    exit 0
fi

# Parse command line args
# Set default for max n samples
max_n_samples=3
while getopts "n:d:p:t:h" opt; do
    case $opt in
        n)
            max_n_samples=$OPTARG
        ;;
        d)
            data_type=$OPTARG
        ;;
        p)
            db_path=$OPTARG
        ;;
        t)
            tsv_path=$OPTARG
        ;;
        h)
            usage
            exit 0
        ;;
        \?)
            usage
            exit 0
        ;;
    esac
done

# Prepare sqlite commands
# Cycle through contigs
for i in {{1..22},X,Y};
do
  if [ $data_type == "exomes" ]; then
    # Fields in exomes databases are:
    # 'id', 'chrom', 'pos', 'ref', 'alt', 'het_or_hom_or_hemi',
    # 'variant_id', 'n_expected_samples', 'n_available_samples'
    # gnomAD v2.1 has 1000 sqlite databases per chromosome
    files=(
        ${tsv_path}/combined_chr${i}_*.db
    )
    for j in "${files[@]}"
    do
      part=$(basename $j | awk -F '[_.]' '{print $3}')
      output_tsv=${tsv_path}/chr${i}_${part}.tsv.gz
      set -x
      (
        echo -e "chrom\tpos\tref\talt\thet_or_hom_or_hemi\tn_available_samples" | gzip -c - > $output_tsv
        sqlite3 -separator $'\t' $j "SELECT chrom, pos, ref, alt, het_or_hom_or_hemi, n_available_samples FROM t WHERE n_available_samples < ${max_n_samples}" | gzip -c -  >> $output_tsv
      ) &> ${tsv_path}/exomes.log &
      set +x
    done
  else
    # Base for gnomAD genomes databases is 'all_variants_s42811_gs50_gn857'
    # Fields in genomes databases are:
    # 'id', 'chrom', 'pos', 'ref', 'alt', 'zygosity',
    # 'qual', 'combined_bamout_id', 'read_group_id'
    files=(
        ${tsv_path}/all_variants_s42811_gs50_gn857.chr${i}_*.db
    )
    for j in "${files[@]}"; do
      output_tsv=${tsv_path}/chr${i}.tsv.gz
      set -x
      (
        echo -e "chrom\tpos\tref\talt\tzygosity\tn_available_samples" | gzip -c - > $output_tsv
        sqlite3 -separator $'\t' $j "SELECT * FROM (SELECT chrom, pos, ref, alt, zygosity, count(*) AS c FROM variants GROUP BY chrom, pos, ref, alt, zygosity)  WHERE c < ${max_n_samples}" | gzip -c -  >> $output_tsv
      ) &> ${tsv_path}/genomes.log &
    done
  fi
done