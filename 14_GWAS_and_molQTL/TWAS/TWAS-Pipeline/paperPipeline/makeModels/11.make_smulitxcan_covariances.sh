#!/bin/bash
module load python/3.7.6-gcc-4.8.5-anaconda

meta_covariance="~/TWAS/GWAS_imputation/summary-gwas-imputation-master/src/meta_covariance_for_models_current.py"
type=splicing
##
python3 ${meta_covariance} \
-parquet_genotype_folder ~/TWAS/GWAS_imputation/genotype \
-parquet_genotype_pattern "PGRP_no_MAF_monoallelic_variants.chr(.*).variants.parquet" \
-model_db_folder ~/TWAS/GTEx_models/type/${type}/dbs \
-model_db_file_pattern "PigGTEx_V1_(.*)_${type}_ElasticNet_models_filtered_signif.db" \
-output ~/TWAS/GWAS_imputation/smultixcan_covariances/${type}.smultixcan_covariances.txt.gz \
-parsimony 9