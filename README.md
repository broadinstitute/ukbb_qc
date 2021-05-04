
# UKBB QC: Densified sparse
This branch was created in November 2019 for the 200K tranche (data freeze 5). The 200K tranche was the first data tranche to be joint called using the hail [vcf combiner](https://hail.is/docs/0.2/experimental/vcf_combiner.html#vcf-combiner) and was therefore the first sparse UKBB dataset.

Due to time pressure from the external partners, we decided to [densify](https://hail.is/docs/0.2/experimental/vcf_combiner.html#hail.experimental.densify) the full 200K dataset and run the dense dataset through our production code rather than re-engineer the pipeline to be sparse-aware. This branch contains the code used to create the 200K data release.

## Data Loading

### Densify
The first step in processing the 200K was densifying the sparse MatrixTable (MT).

Cluster configuration:
```
hailctl dataproc --beta start --autoscaling-policy=autoscale_densify kc1 --master-machine-type n1-highmem-8 --worker-machine-type n1-highmem-8 --init gs://gnomad-public/tools/inits/master-init.sh --max-idle=30m --project maclab-ukbb --properties=spark:spark.executor-memory=25g,spark:spark.speculation=true,spark:spark.speculation.quantile=0.9,spark:spark.speculation.multiplier=3 --master-boot-disk-size 400 --preemptible-worker-boot-disk-size 400 --worker-boot-disk-size 400
```

Command:
```
hailctl dataproc submit kc1 densify_sparse_mt_1.py -i gs://broad-pharma5-ukbb-outputs/hail_dataproc_20191108115937 -f 5 -n 25000 --overwrite
```
Job ID: `9a79da459e8d4516b9060aba3aea38c7`
Job page: https://console.cloud.google.com/dataproc/jobs/9a79da459e8d4516b9060aba3aea38c7?region=us-central1&project=maclab-ukbb
Job completed in 22 hours 16 minutes.

### Sites VCF
The second step of data loading was to generate a site-only VCF to pass to DSP (input necessary for VQSR).

Cluster:
```
hailctl dataproc --beta start --autoscaling-policy=autoscale_workers kc1 --master-machine-type n1-highmem-8 --worker-machine-type n1-highmem-8 --init gs://gnomad-public/tools/inits/master-init.sh --max-idle=60m --worker-boot-disk-size=100 --project maclab-ukbb --properties=spark:spark.executor-memory=25g,spark:spark.speculation=true,spark:spark.speculation.quantile=0.9,spark:spark.speculation.multiplier=3 --master-boot-disk-size 400
```

Command (script has been renamed to `create_sites_vcf_from_sparse.py`): 
```
hailctl dataproc submit kc1 create_sites_vcf.py --input gs://broad-pharma5-ukbb-outputs/hail_dataproc_20191108115937 -f 5
```
Job ID: `48c5bd7cc0904f0cbecb331c49244d25`
Job page: https://console.cloud.google.com/dataproc/jobs/48c5bd7cc0904f0cbecb331c49244d25?region=us-central1&project=maclab-ukbb
Job completed in 1 hour 40 minutes.


## Sample QC
### Impute sex
Cluster:
```
hailctl dataproc start kc2 --master-machine-type n1-highmem-8 --worker-machine-type n1-highmem-8 --init gs://broad-ukbb/broad.freeze_5/temp/temp_init.sh --max-age=480m --worker-boot-disk-size=100 --project maclab-ukbb --properties=spark:spark.executor-memory=25g,spark:spark.speculation=true,spark:spark.speculation.quantile=0.9,spark:spark.speculation.multiplier=3 --packages=holoviews,statsmodels,matplotlib --num-workers=60
```

Command:
```
hailctl dataproc submit kc2 generate_hardcalls.py --impute_sex
```

Job ID: `f1e36aecb06c437dab15522e44fb70fd`
Job page: https://console.cloud.google.com/dataproc/jobs/f1e36aecb06c437dab15522e44fb70fd?region=us-central1&project=maclab-ukbb
Job completed in 45 minutes 43 seconds.

### Generate hardcalls
Cluster:
```
hailctl dataproc start kc1 --master-machine-type n1-highmem-8 --worker-machine-type n1-highmem-8 --init gs://broad-ukbb/broad.freeze_5/temp/temp_init.sh --max-idle=30m --worker-boot-disk-size=100 --project maclab-ukbb --properties=spark:spark.executor-memory=25g,spark:spark.speculation=true,spark:spark.speculation.quantile=0.9,spark:spark.speculation.multiplier=3 --packages=holoviews,statsmodels,matplotlib,seaborn --num-workers=60
hailctl dataproc modify kc1 --num-workers=100
```

Command: 
```
hailctl dataproc submit kc1 generate_hardcalls.py --split_hardcalls --overwrite
```

Job ID: `4bd0dcd7644f4083bef3ef6ef7da28f2`
Job page: https://console.cloud.google.com/dataproc/jobs/4bd0dcd7644f4083bef3ef6ef7da28f2?region=us-central1&project=maclab-ukbb
Job completed in 7 hours 12 minutes.

### Hard filters 
Cluster:
```
hailctl dataproc start kc --master-machine-type n1-highmem-8 --worker-machine-type n1-highmem-8 --init gs://broad-ukbb/broad.freeze_5/temp/temp_init.sh --max-idle=30m --worker-boot-disk-size=100 --project maclab-ukbb --properties=spark:spark.executor-memory=25g,spark:spark.speculation=true,spark:spark.speculation.quantile=0.9,spark:spark.speculation.multiplier=3 --packages=holoviews,statsmodels,matplotlib,seaborn --num-workers=100
```

Command:
```
hailctl dataproc submit kc apply_hard_filters.py --overwrite
```
Job ID: `8a8a59d429114749842671efb8e602c8`
Job page: https://console.cloud.google.com/dataproc/jobs/8a8a59d429114749842671efb8e602c8?region=us-central1&project=maclab-ukbb
Job completed in 36 seconds.

### Platform PCA
Same cluster as hard filters above.

#### Run platform PCA
Command:
```
hailctl dataproc submit kc platform_pca.py --compute_callrate_mt --overwrite
```

Job ID: `1390ea55e00f491dbe7d306335a63888`
Job page: https://console.cloud.google.com/dataproc/jobs/1390ea55e00f491dbe7d306335a63888?region=us-central1&project=maclab-ukbb&ref=https:%2F%2Fconsole.cloud.google.com%2Fdataproc%2Fjobs
Job completed in 8 minutes 39 seconds.

#### Assign platforms
Same cluster as hard filters above.
Command:
```
hailctl dataproc submit kc platform_pca.py --assign_platforms --overwrite
```

Job ID: ``

### Relatedness inference
Cluster:
```
hailctl dataproc start kc1 --master-machine-type n1-highmem-8 --worker-machine-type n1-highmem-8 --init gs://broad-ukbb/broad.freeze_5/temp/temp_init.sh --max-idle=30m --worker-boot-disk-size=100 --project maclab-ukbb --properties=spark:spark.executor-memory=25g,spark:spark.speculation=true,spark:spark.speculation.quantile=0.9,spark:spark.speculation.multiplier=3 --packages=holoviews,statsmodels,matplotlib,seaborn --num-workers=60 --num-worker-local-ssds 1 --properties=spark:spark.task.maxFailures=20,spark:spark.driver.memory=40g,hdfs:dfs.replication=1,spark:spark.executor.extraJavaOptions=-Xss4M,spark:spark.driver.extraJavaOptions=-Xss4M,spark:spark.kryoserializer.buffer.max=1g,spark:spark.driver.maxResultSize=0,spark:spark.executor.memory=15g
```

Command:
```
hailctl dataproc submit kc1 relatedness_inference.py --overwrite
```

Job ID: `4363b65b76cb4ae1bf269bd4695159e7`
Crashed after 7 minutes 56 seconds: ran PCA for PC relate but crashed on typo.

New cluster:
```
hailctl dataproc start kc1 --master-machine-type n1-highmem-8 --worker-machine-type n1-highmem-8 --init gs://broad-ukbb/broad.freeze_5/temp/temp_init.sh --max-idle=30m --worker-boot-disk-size=100 --project maclab-ukbb --properties=spark:spark.executor-memory=25g,spark:spark.speculation=true,spark:spark.speculation.quantile=0.9,spark:spark.speculation.multiplier=3 --packages=holoviews,statsmodels,matplotlib,seaborn --num-workers=60
```

Command: 
```
hailctl dataproc submit kc1 relatedness_inference.py --skip_pc_relate --skip_filter_dups --overwrite
```
Job ID:  `481dea32f76b4198a48b69a2a72d77d4`
Job page: https://console.cloud.google.com/dataproc/jobs/481dea32f76b4198a48b69a2a72d77d4?region=us-central1&project=maclab-ukbb
Job completed in 1 minute 43 seconds.

### Population inference
Same cluster as relatedness above.

Commands:
```
hailctl dataproc submit kc1 population_pca.py --run_pca --overwrite
hailctl dataproc submit kc1 population_pca.py --assign_clusters --overwrite
hailctl dataproc submit kc1 population_pca.py --assign_clusters_joint_scratch_array_pcs --overwrite
hailctl dataproc submit kc1 population_pca.py --run_pc_project --overwrite
hailctl dataproc submit kc1 population_pca.py --run_rf --overwrite
hailctl dataproc submit kc1 population_pca.py --assign_hybrid_ancestry --overwrite
```

Job IDs:
[60a461db28d44f9db5cedefab3d569aa](https://console.cloud.google.com/dataproc/jobs/60a461db28d44f9db5cedefab3d569aa?region=us-central1&project=maclab-ukbb): Completed in 12 minutes 10 seconds
[3dfe43f66a14400490e11ccf13b8ddc0](https://console.cloud.google.com/dataproc/jobs/3dfe43f66a14400490e11ccf13b8ddc0?region=us-central1&project=maclab-ukbb): Completed in 2 minutes 30 seconds
[93645dd56ecc4bfd8937dd16b493e90f](https://console.cloud.google.com/dataproc/jobs/93645dd56ecc4bfd8937dd16b493e90f?region=us-central1&project=maclab-ukbb): Completed in 5 minutes 46 seconds
[2abcc51a4fc246dcb5c690c3a3e887fe](https://console.cloud.google.com/dataproc/jobs/2abcc51a4fc246dcb5c690c3a3e887fe?region=us-central1&project=maclab-ukbb): Completed in 1 hour 6 minutes
[e09177c09ba14aeb81f6c37b8a1322f6](https://console.cloud.google.com/dataproc/jobs/e09177c09ba14aeb81f6c37b8a1322f6?region=us-central1&project=maclab-ukbb): Completed in 2 minutes 59 seconds
[434cc2eb5e3647d28f21d7274f5ee018](https://console.cloud.google.com/dataproc/jobs/434cc2eb5e3647d28f21d7274f5ee018?region=us-central1&project=maclab-ukbb): Completed in 45 seconds

### Outlier detection
Same cluster as relatedness above.

Commands:
```
hailctl dataproc submit kc1 outlier_filter.py --run_mini_qc --overwrite
hailctl dataproc submit kc1 outlier_filter.py --skip_platform_filter --overwrite
hailctl dataproc submit kc1 outlier_filter.py --overwrite
```

Job IDs:
 02b26328594742be976cd45d27f3cb6f
62e33a0cc2a64d6db9eb4978d408ae7e
a098d8b3cd074f28bf539cdf8f74c795

### Meta HT