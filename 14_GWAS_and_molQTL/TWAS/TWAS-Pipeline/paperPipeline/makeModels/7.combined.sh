#!/bin/bash
tissue=$1
type=$2

cd ~/TWAS/GTEx_models/type/${type}/summary/${tissue}
for tt in {1..18}
do
	cd ${tt}
	cat ${type}.${tissue}_Model_training_chr${tt}.*.model_summaries.txt > ../${type}.${tissue}_Model_training_chr${tt}_model_summaries.txt
	sed -i '1i\gene_id\tgene_name\tgene_type\talpha\tn_snps_in_window\tn_snps_in_model\tlambda_min_mse\ttest_R2_avg\ttest_R2_sd\tcv_R2_avg\tcv_R2_sd\tin_sample_R2\tnested_cv_fisher_pval\trho_avg\trho_se\trho_zscore\trho_avg_squared\tzscore_pval\tcv_rho_avg\tcv_rho_se\tcv_rho_avg_squared\tcv_zscore_est\tcv_zscore_pval\tcv_pval_est' ../${type}.${tissue}_Model_training_chr${tt}_model_summaries.txt
	rsync -av ${type}.${tissue}_Model_training_chr${tt}_summary.txt ../
	cd ..
done

cd ~/TWAS/GTEx_models/type/${type}/covariances/${tissue}
for ttt in {1..18}
do
	cd ${ttt}
	cat ${type}.${tissue}_Model_training_chr${ttt}.*.covariances.txt > ../${type}.${tissue}_Model_training_chr${ttt}_covariances.txt
	sed -i '1i\GENE RSID1 RSID2 VALUE' ../${type}.${tissue}_Model_training_chr${ttt}_covariances.txt
	cd ..
done

cd ~/TWAS/GTEx_models/type/${type}/weights/${tissue}
for tttt in {1..18}
do
	cd ${tttt}
	cat ${type}.${tissue}_Model_training_chr${tttt}.*.weights.txt > ../${type}.${tissue}_Model_training_chr${tttt}_weights.txt
	sed -i '1i\gene_id\trsid\tvarID\tref\talt\tbeta' ../${type}.${tissue}_Model_training_chr${tttt}_weights.txt
	cd ..
done
