#!/bin/bash
userpath=$1
cd ${userpath}
###1. Call Python program
module load python/3.7.6-gcc-4.8.5-anaconda

###2. Set the corresponding path
Traits=$2 # Extract trait_name
type=$3
snp_covariance=$4
SPrediXcan_result=$5
model_path=$6
summary_path=$7
ouput=$8
METAXCAN=$9
###3. Parallel running program
fifoFile="test_fifo"			
rm -f ${fifoFile}				
mkfifo ${fifoFile}						
#				
exec 9<> ${fifoFile}			
rm -f ${fifoFile}						
#		
for ((i=0;i<18;i++))					
do										
    echo "" >&9								
done										
echo "wait all task finish,then exit!!!"

for trait in ${Traits[*]}
do
read -u9        
{                
	python $METAXCAN/SMulTiXcan.py \
	--models_folder ${model_path} \
	--models_name_filter "PigGTEx_V1_(.*)_${type}_ElasticNet_models_filtered_signif.db" \
	--models_name_pattern "PigGTEx_V1_(.*)_${type}_ElasticNet_models_filtered_signif.db" \
	--snp_covariance ${snp_covariance}/${type}.smultixcan_covariances.txt.gz \
	--metaxcan_folder ${SPrediXcan_result}/${trait}/ \
	--metaxcan_filter "${type}.${trait}.(.*).csv" \
	--metaxcan_file_name_parse_pattern "${type}.${trait}.(.*).csv" \
	--gwas_file ${summary_path}/${trait}.txt.gz \
	--snp_column panel_variant_id \
	--effect_allele_column effect_allele \
	--non_effect_allele_column non_effect_allele \
	--beta_column effect_size \
	--pvalue_column pvalue \
	--keep_non_rsid \
	--model_db_snp_key varID \
	--cutoff_condition_number 30 \
	--verbosity 7 \
	--throw \
	--ouput $ouput/${trait}/${type}.${trait}.SMultixcan.txt
echo "" >&9    
} & 
done

wait            
exec 9>&-		 
echo			 
echo "success"   