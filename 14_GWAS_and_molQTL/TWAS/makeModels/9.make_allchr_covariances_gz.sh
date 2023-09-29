#!/bin/sh
###Set workspace
work_path=$1
type=$2
Tissues=$3
cd ${work_path}

###Cycle submit task
for tissue in ${Tissues[*]}
do
	for ttt in {1..18}
	do
		awk 'NR>1{print $0}' ./${type}/covariances/${tissue}/${type}.${tissue}_Model_training_chr${ttt}_covariances.txt >> ./${type}/covariances/${tissue}/chr.tmp #Remove column names and merge files
	done
	sed -i '1i\GENE RSID1 RSID2 VALUE' ./${type}/covariances/${tissue}/chr.tmp #Add column name
	gzip -c ./${type}/covariances/${tissue}/chr.tmp > ./${type}/covariances/${tissue}/${type}.${tissue}_Model_training_covariances.txt.gz
	rm ./${type}/covariances/${tissue}/chr.tmp #Delete temporary files
done
