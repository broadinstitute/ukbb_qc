import hail as hl

from ukbb_qc.resources.basics import dup_gvcf_path, dup_mt_path


hl.init(log="/combiner.log", default_reference="GRCh38")

"""
We noticed that the freeze 7/450k callset had duplicate sample IDs when starting the data loading steps.

The duplicate samples had two different gVCF files, and the pharma partners requested that we use the most recent gVCF version for each duplicate.

We used the hail gVCF combiner to create a MatrixTable containing only the most recent version of each duplicate sample to resolve
the duplicates in the 450k MT.
"""
try:
    freeze = 7

    # Input file
    # 2 columns: sample and gvcf path
    path_to_input_list = dup_gvcf_path(freeze)

    inputs = []
    with hl.hadoop_open(path_to_input_list, "r") as f:
        for line in f:
            # line has two columns:
            # line[0] is the sample ID
            # line[1] is the gVCF path
            # the combiner only requires the gVCF path as input
            line = line.strip().split("\t")
            inputs.append(line[1])

    # temp bucket for combiner and log
    temp_bucket = f"gs://broad-ukbb/broad.freeze_{freeze}/temp/combiner_temp/"

    hl.experimental.run_combiner(
        inputs,
        out_file=dup_mt_path(freeze),
        tmp_path=f"gs://broad-ukbb/broad.freeze_{freeze}/temp/combiner_temp/",
        key_by_locus_and_alleles=True,
        reference_genome="GRCh38",
        use_exome_default_intervals=True,
    )

finally:
    hl.copy_log(temp_bucket)
