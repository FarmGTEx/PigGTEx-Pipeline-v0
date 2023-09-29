#!/bin/sh
###Set workspace
work_path=$1## Other QTL were also treated
type=$2
geno_path=$3
tmm_path=$4

cd ${work_path}
mkdir ./type/${type}/tissue_expression
mkdir ./type/${type}/tissue_snp_annotation
###Extraction organization name
tissue_list=(`ls ${tmm_path}/*.gz`)
tissue_list=(${tissue_list[*]#*raw_expression/})
tissue_list=(${tissue_list[*]%.pcg_*})
Tissues=($( echo ${tissue_list[*]} | sed 's/ /\n/g' | sort -u ))
###Cycle submit task
for tissue in ${Tissues[*]} #Ensure that the number of server task nodes is sufficient
do
	##Extract the corresponding information columns from the VCF file to generate a whole chromosome SNP annotation file
	mkdir -p ./type/${type}/tissue_snp_annotation/${tissue} 
	bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%ID\n' ${geno_path}/${tissue}.filtered_maf0.05_mac6.geno.vcf.gz > ./type/${type}/tissue_snp_annotation/${tissue}/${tissue}.snp_annotation.txt
	##Prepare corresponding gene expression files
	#Extract the corresponding files and delete unnecessary information
	gunzip -c ${tmm_path}/${tissue}.pcg_expr_tmm_inv.bed.gz | awk '{for(i=4;i<NF;i++)printf("%s ",$i);print $NF}' > ./type/${type}/tissue_expression/${tissue}.expr.tmp1.txt
	#File transposition
	awk '{i=1;while(i <= NF){col[i]=col[i] $i "\t";i=i+1}} END {i=1;while(i<=NF){print col[i];i=i+1}}' ./type/${type}/tissue_expression/${tissue}.expr.tmp1.txt | sed 's/[ \t]*$//g' > ./type/${type}/tissue_expression/${tissue}.expr.tmp2.txt
	#Delete column name：gene_id
	sed 's/gene_id//' ./type/${type}/tissue_expression/${tissue}.expr.tmp2.txt > ./type/${type}/tissue_expression/${tissue}.${type}.transformed_expression.txt
	rm ./type/${type}/tissue_expression/${tissue}.expr.tmp*
	###Processing files requiring chromosome sorting
	for chr in {1..18}
	do
		yhbatch -N 1 -p bigdata -J ${chr}${tissue} 2.tissues_chr_jobs.sh ${work_path} ${chr} ${tissue} ${geno_path} ${type}
	done
done
