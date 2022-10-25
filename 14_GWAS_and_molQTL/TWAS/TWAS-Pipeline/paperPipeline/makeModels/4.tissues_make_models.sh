#!/bin/sh
###Call the corresponding R version
module load R/4.1.0-gcc-4.8.5
###Set workspace
work_path='~/TWAS/GTEx_models/software/code'
cd ${work_path}
###Extraction string
tissue_list=(`ls ~/TWAS/GTEx_models/type/eQTL/tissue_cover`)
Tissues=($(echo ${tissue_list[*]} | sed 's/.covariates.txt//g' | sort -u))
type='splicing'
###Cycle submit task
#The sample size of muscle is too large. Run it at another node
type=eQTL
for tissue in ${Tissues[*]}
do
	if [[ ${tissue} = "Muscle" ]] 
	then
		for chrom in {1..18}
		do
			yhbatch -N 1 -p MEM_128 -J eqtl${tissue} R CMD BATCH "--args ${chrom} ${tissue} ${type}" ./parallel_tiss_chrom_training.R
			sleep 10
		done
	else
		for chrom in {1..18}
		do
			yhbatch -N 1 -p bigdata -J eqtl${tissue} R CMD BATCH "--args ${chrom} ${tissue} ${type}" ./parallel_tiss_chrom_training.R
			sleep 10
		done
	fi
done
