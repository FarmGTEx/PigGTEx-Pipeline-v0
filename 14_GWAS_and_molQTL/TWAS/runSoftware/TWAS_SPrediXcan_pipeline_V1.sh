#!/bin/sh
userpath=$1
cd ${userpath}
###1. Receive external parameters
tissue=$2
trait=$3
type=$4


###2. Set the corresponding path
summary_path=$5
tmp_path=$6
snp_path=$7
covariance_path=$8
model_path=$9
result_path=$10
SPrediXcan=$11


###3. Extract the coincident part from the summary file
## The sample data is a compressed file harmonized by software(https://github.com/hakyimlab/summary-gwas-imputation),
## and please adjust the code according to the actual data.
# The extracted column name: variant_id(snpID) effect_allele(alt) non_effect_allele(ref) pvalue effect_size(beta)
zcat ${summary_path}/${trait}.txt.gz | awk '{print $1,$5,$6,$8,$10}' > ${tmp_path}/${tissue}.${trait}.new.summary.tmp 

# Save original file column name
sed -n '1p' ${tmp_path}/${tissue}.${trait}.new.summary.tmp > ${tmp_path}/${tissue}.${trait}.colnames.tmp

## In order to ensure the subsequent calculation, the snpID of the summary file 
#  needs to be consistent with the snpID of the corresponding tissue site information file(${snp_path}),
## and please replace it according to chr and pos.

# Extract the coincident part of two files
awk -F "\t" 'FNR==NR {a[$1]=$1; next}($1 in a){print $0}' ${snp_path}/${tissue}.snp_chr_pos.txt ${tmp_path}/${tissue}.${trait}.new.summary.tmp > ${tmp_path}/${tissue}.${trait}.overlap.tmp

# Splice column names back into the coincident files
cat ${tmp_path}/${tissue}.${trait}.colnames.tmp ${tmp_path}/${tissue}.${trait}.overlap.tmp > ${tmp_path}/${tissue}.${trait}.overlap.txt

# Delete temporary files
rm ${tmp_path}/${tissue}.${trait}.*.tmp 


###4. Run S-PrediXcan software
# Call Python program
module load python/3.7.6-gcc-4.8.5-anaconda

## Please adjust the input according to the actual column
## name of the corresponding file(${tissue}.${trait}.overlap.txt).
# Run S-PrediXcan
python3 ${SPrediXcan}/SPrediXcan.py \
--model_db_path ${model_path}/PigGTEx_V1_${tissue}_${type}_ElasticNet_models.db \
--covariance ${covariance_path}/${type}.${tissue}_Model_training_covariances.txt \
--gwas_file ${tmp_path}/${tissue}.${trait}.overlap.txt \
--snp_column variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--beta_column effect_size \
--pvalue_column pvalue \
--keep_non_rsid \
--overwrite \
--output_file ${result_path}/${type}.${trait}.${tissue}.csv

# Delete temporary files
rm ${tmp_path}/${tissue}.${trait}.overlap.txt
