import hail as hl

hl.init(log='/combiner.log', default_reference='GRCh38')

try:
    # Input file
    # 2 columns: sample and gvcf path
    path_to_input_list = 'gs://broad-ukbb/broad.freeze_7/temp/duplicate_sample_map_no_UU.tsv' 

    inputs = []
    with hl.hadoop_open(path_to_input_list, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            inputs.append(line[1])

    # path to output MT
    output_file = 'gs://broad-ukbb/broad.freeze_7/temp/most_recent_dup.mt'

    # temp bucket
    temp_bucket = 'gs://broad-ukbb/broad.freeze_7/temp/combiner_temp/'
    hl.experimental.run_combiner(
        inputs, 
        out_file=output_file, 
        tmp_path=temp_bucket, 
        key_by_locus_and_alleles=True,
        reference_genome='GRCh38',
        use_exome_default_intervals=True,
    )

finally:
    hl.copy_log(temp_bucket)
