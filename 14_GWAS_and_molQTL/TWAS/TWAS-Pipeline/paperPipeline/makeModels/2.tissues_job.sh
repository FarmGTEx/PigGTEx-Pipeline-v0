#!/bin/sh
###Set workspace
type=eQTL
work_path=~/TWAS/GTEx_models/type/${type} ## Other QTL were also treated
cd ${work_path}
###Extraction organization name
tissue_list=(`ls ../tissue_expression/raw_expression/*.gz`)
tissue_list=(${tissue_list[*]#*raw_expression/})
tissue_list=(${tissue_list[*]%.pcg_*})
Tissues=($( echo ${tissue_list[*]} | sed 's/ /\n/g' | sort -u ))
###Cycle submit task
for tissue in ${Tissues[*]} #Ensure that the number of server task nodes is sufficient
do
	##Extract the corresponding information columns from the VCF file to generate a whole chromosome SNP annotation file
	mkdir -p ../tissue_snp_annotation/${tissue} 
	bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%ID\n' ../raw_tissue_gene_data/raw_QC_vcf_data/${tissue}.filtered_maf0.05_mac6.geno.vcf.gz > ../tissue_snp_annotation/${tissue}/${tissue}.snp_annotation.txt
	##Prepare corresponding gene expression files
	#Extract the corresponding files and delete unnecessary information
	gunzip -c ../tissue_expression/raw_expression/${tissue}.pcg_expr_tmm_inv.bed.gz | awk '{for(i=4;i<NF;i++)printf("%s ",$i);print $NF}' > ../tissue_expression/${tissue}.expr.tmp1.txt
	#File transposition
	awk '{i=1;while(i <= NF){col[i]=col[i] $i "\t";i=i+1}} END {i=1;while(i<=NF){print col[i];i=i+1}}' ../tissue_expression/${tissue}.expr.tmp1.txt | sed 's/[ \t]*$//g' > ../tissue_expression/${tissue}.expr.tmp2.txt
	#Delete column name：gene_id
	sed 's/gene_id//' ../tissue_expression/${tissue}.expr.tmp2.txt > ../tissue_expression/${tissue}.${type}.transformed_expression.txt
	rm ../tissue_expression/${tissue}.expr.tmp*
	###Processing files requiring chromosome sorting
	for chr in {1..18}
	do
		yhbatch -N 1 -p bigdata -J ${chr}${tissue} ./tissues_chr_jobs.sh ${chr} ${tissue} #Transfer the tissue name and chromosome parameters into the script
	done
done
