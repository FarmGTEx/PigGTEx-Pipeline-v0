#!/bin/sh
###Set workspace
work_path=$1
chr=$2
tissue=$3
geno_path=$4
type=$5
cd ${work_path}
###The corresponding SNP annotation files were processed by chromosome division
#mkdir -p ./type/${type}/tissue_snp_annotation/${tissue}
##Extract the information column of the corresponding chromosome and set the interval to tab
awk '{if($1=='${chr}') print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ./type/${type}/tissue_snp_annotation/${tissue}/${tissue}.snp_annotation.txt > ./type/${type}/tissue_snp_annotation/${tissue}/${tissue}.snp_annot.chr${chr}.txt
##Add the corresponding column name
sed -i '1i\chromosome\tpos\tvarID\tref_vcf\talt_vcf\trsid' ./type/${type}/tissue_snp_annotation/${tissue}/${tissue}.snp_annot.chr${chr}.txt

###Processing genotype files by chromosome
mkdir -p ${geno_path}/tissue_chr_data/${tissue}
mkdir -p ./tissue_gene_matrix/${tissue}
##The genotype 0 1 2 data of the corresponding chromosomes were extracted and encoded with ref
plink --vcf ${geno_path}/${tissue}.filtered_maf0.05_mac6.geno.vcf.gz --const-fid --keep-allele-order --chr ${chr} --recode A --out ${geno_path}/tissue_chr_data/${tissue}/${tissue}.ref.chr${chr}
##Delete unwanted columns
cat ${geno_path}/tissue_chr_data/${tissue}/${tissue}.ref.chr${chr}.raw | cut -d" " -f2,7- > ./tissue_gene_matrix/${tissue}/${tissue}.genotype_chr${chr}.txt
##Delete unnecessary endnotes (Plink's raw file will generate an endnote of _ [A-Z])
sed -i '1s/_[A-Z] / /g' ./tissue_gene_matrix/${tissue}/${tissue}.genotype_chr${chr}.txt
sed -i '1s/_[A-Z]$//' ./tissue_gene_matrix/${tissue}/${tissue}.genotype_chr${chr}.txt
##Replace IID with ID
sed -i '1s/IID/Id/' ./tissue_gene_matrix/${tissue}/${tissue}.genotype_chr${chr}.txt
##Transpose file
awk '{i=1;while(i <= NF){col[i]=col[i] $i "\t";i=i+1}} END {i=1;while(i<=NF){print col[i];i=i+1}}' ./tissue_gene_matrix/${tissue}/${tissue}.genotype_chr${chr}.txt | sed 's/[ \t]*$//g' > ./tissue_gene_matrix/${tissue}/${tissue}.genotype.chr${chr}.txt
rm ./tissue_gene_matrix/${tissue}/${tissue}.genotype_chr${chr}.txt #Delete temporary files
