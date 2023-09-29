# PigGTEx: Pipeline for TWAS

Pipeline for TWAS using MetaXcan.

#### **Descriptions:**

1. Prepare models of TWAS
2. Run S-Predixcan and S-MultiXcan

#### Software: ==Please confirm the versions and download links.==
- `MetaXcan` : https://github.com/xqwen/torus
- `summary-gwas-imputation` :  https://github.com/hakyimlab/summary-gwas-imputation

#### Attached scripts and files:

- 1.makeAnnotation.R
- 2.tissues_job.sh
- 3.tissues_chr_jobs.sh
- 4.tissues_make_models.sh
- 5.parallel_tiss_chrom_training.R
- 6.parallel_nested_cv_elnet.R
- 7.combined.sh
- 8.make_db.R
- 9.make_allchr_covariances_gz.sh
- 10.meta_covariance_for_models_current.py
- 11.make_smulitxcan_covariances.sh
- TWAS_SPrediXcan_pipeline_V1.sh
- TWAS_SMulTiXcan_pipeline_V1.sh



## Manuals and Codes

### ###============================Code of model file================================###

#### 1. Extract the corresponding information from the data on Ensembl and prepare the gene annotation file

```R
#==================================================#
# >> Input
# work_dir  # Original comment file directory
# out_dir  # Directory of newly generated comment files
# << Output
# *_annot.txt  # New annotation file in software format
#==================================================#
R CMD BATCH "--args ${work_dir} ${out_dir}" 1.makeAnnotation.R
```

#### 2. Prepare corresponding data

```shell
#==================================================#
# >> Input
# work_dir  # Work Path
# type  # Type of QTLs
# geno_dir # Genotype file paths by tissue
# tmm_dir # Transcriptome data file path for each tissue
# << Output
# ${tissue}.snp_annot.chr${chr}.txt  # SNP annotation files for each chromosome in each tissue
# ${tissue}.${type}.transformed_expression.txt # Various types of transcriptome files from various tissues
# ${tissue}.genotype.chr${chr}.txt # Genotype files for each chromosome in each tissue
#==================================================#
shell 2.tissues_job.sh ${work_dir} ${type} ${geno_dir} ${tmm_dir}
####  When this script is run, it will automatically run 3.tissues_chr_jobs.sh
```

#### 3. Building models

##### 3.1 Parallel submission of tasks to the server

```shell
#==================================================#
# >> Input
# work_dir  # Work Path
# type  # Type of QTLs
# covar_dir # Paths to covariate files for each type of tissue
# software_dir # https://github.com/hakyimlab/PredictDBPipeline
# << Output
# The set of files required to build the model are detailed in the corresponding GitHub URLs
#==================================================#
shell 4.tissues_make_models.sh ${work_dir} ${type} ${covar_dir}
```

```R
#==================================================#
# The author's source code is modified so that it can run in parallel
#==================================================#
5.parallel_tiss_chrom_training.R # This script will automatically run in 4.tissues_make_models.sh
6.parallel_nested_cv_elnet.R #Modified the author's source code, and added the transfer parameters of 'tissue' and 'type' (to determine the location of the result file)
#For more information, see https://github.com/hakyimlab/PredictDBPipeline
```

##### 3.2 Merge the result files of 3.1

```shell
#==================================================#
# >> Input
# work_dir  # The output location of the result file from 3.1
# type  # Type of QTLs
# tissue # Name of individual tissue
# << Output
# Merge all chromosome files
#==================================================#
shell 7.combined.sh ${work_dir} ${type} ${tissue}
```

##### 3.3 Build the model DB file

```R
#==================================================#
# >> Input
# type  # Type of QTLs
# tissue # Name of individual tissue
# << Output
# PigGTEx_V1_${tissue}_${type}_ElasticNet_models 
# PigGTEx_V1_${tissue}_${type}_ElasticNet_models_filtered_signif # Various types of model DB files for each tissue
#==================================================#
R CMD BATCH "--args ${type} ${tissue}" 8.make_db.R
```

#### 4. Combine covariate files for all chromosomes of the model

```shell
#==================================================#
# >> Input
# work_dir  # The output location of the result file from 3.1
# type  # Type of QTLs
# tissueList # A list with all tissue names
# << Output
# ${type}.${tissue}_Model_training_covariances.txt.gz # Covariate files for each type of model for each tissue after consolidation
#==================================================#
shell 9.make_allchr_covariances_gz.sh ${work_dir} ${type} ${tissueList}
```

#### 5. Construction of multi-tissue covariate files

```shell
#==================================================#
# >> Input
# software_dir  # https://github.com/hakyimlab/summary-gwas-imputation
# type  # Type of QTLs
# parquet_genotype_folder # Path to parquet_genotype
# model_dir # Path to the model file
# result_dir # Result Path
# << Output
# ${type}.smultixcan_covariances.txt.gz # Multi-tissue covariate files
#==================================================#
shell 11.make_smulitxcan_covariances.sh ${software_dir} ${type} ${parquet_genotype_folder} ${model_dir} ${result_dir}
### 10.meta_covariance_for_models_current.py #Modified the author's source code to fit the existing data
```

### ###================================Run TWAS===================================###

#### 1. Run S-PrediXcan

```shell
#==================================================#
# work_dir  # Work Path
# tissue # Name of individual tissue
# trait # Name of individual trait
# type # Type of QTLs
# summary_dir # GWAS results file
# tmp_dir # Temporary file path
# snp_dir # Genotype file path
# covar_dir # Model covariate files
# model_dir # Model file path
# result_dir # Result files path
# software_dir # https://github.com/hakyimlab/MetaXcan
# << Output
# ${type}.${trait}.${tissue}.csv # The result file of S-PrediXcan
#==================================================#
shell TWAS_SPrediXcan_pipeline_V1.sh ${work_dir} ${tissue} ${trait} ${type} ${summary_dir} ${tmp_dir} ${snp_dir} ${covar_dir} ${model_dir} ${result_dir} ${software_dir}
```

#### 2. Run S-MulTiXcan

```shell
#==================================================#
# work_dir  # Work Path
# trait_list # A list of trait names
# type # Type of QTLs
# covar_dir # Multi-tissue covariate files
# SPrediXcan_dir # Path to SPrediXcan's results file
# model_dir # Model file path
# summary_dir # GWAS results file
# result_dir # Result files path
# software_dir # https://github.com/hakyimlab/MetaXcan
# << Output
# ${trait}/${type}.${trait}.SMultixcan.txt# The result file of S-Multixcan
#==================================================#
shell TWAS_SMulTiXcan_pipeline_V1.sh ${work_dir} ${trait_list} ${type} ${covar_dir} ${SPrediXcan_dir} ${model_dir} ${summary_dir} ${result_dir} ${software_dir}
```

