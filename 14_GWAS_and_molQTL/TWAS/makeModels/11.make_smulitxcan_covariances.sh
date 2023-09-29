#!/bin/bash
module load python/3.7.6-gcc-4.8.5-anaconda
software_path=$1
meta_covariance=${software_path}/src/meta_covariance_for_models_current.py
type=$2
parquet_genotype_folder=$3
model_path=$4
result_path=$5

##
python3 ${meta_covariance} \
-parquet_genotype_folder ${parquet_genotype_folder} \
-parquet_genotype_pattern "PGRP_no_MAF_monoallelic_variants.chr(.*).variants.parquet" \
-model_db_folder ${model_path}/dbs \
-model_db_file_pattern "PigGTEx_V1_(.*)_${type}_ElasticNet_models_filtered_signif.db" \
-output ${result_path}/${type}.smultixcan_covariances.txt.gz \
-parsimony 9