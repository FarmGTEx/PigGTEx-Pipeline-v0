#!/bin/sh
###Set workspace
type=exon
work_path=~/TWAS/GTEx_models/type/${type}
cd ${work_path}
###Extraction organization name
tissue_list=(`ls ~/tissue_covar`)
Tissues=($(echo ${tissue_list[*]} | sed 's/.covariates.txt//g' | sort -u))

###Cycle submit task
for tissue in ${Tissues[*]}
do
	for ttt in {1..18}
	do
		awk 'NR>1{print $0}' ./covariances/${tissue}/${type}.${tissue}_Model_training_chr${ttt}_covariances.txt >> ./covariances/${tissue}/chr.tmp #Remove column names and merge files
	done
	sed -i '1i\GENE RSID1 RSID2 VALUE' ./covariances/${tissue}/chr.tmp #Add column name
	gzip -c ./covariances/${tissue}/chr.tmp > ./covariances/${tissue}/${type}.${tissue}_Model_training_covariances.txt.gz
	rm ./covariances/${tissue}/chr.tmp #Delete temporary files
done
