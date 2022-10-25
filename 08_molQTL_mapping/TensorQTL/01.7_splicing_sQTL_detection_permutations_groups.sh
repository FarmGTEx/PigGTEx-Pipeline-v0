#!/bin/bash
# conda activate deepmd
cd /BIGDATA2/scau_hzhang_1/USER/tengjy/Pig_GTEx_splicing/sqtl_mapping
#watch -n 2 nvidia-smi
# module load icu4c/60.1-gcc-4.8.5
mkdir output_splicing

tis_main_final=("Adipose" "Large_intestine" "Small_intestine" "Blood" "Fetal_thymus" "Lymph_node" "Milk" "Spleen" "Cartilage" "Synovial_membrane" "Brain" "Embryo" "Heart" "Kidney" "Liver" "Lung" "Testis" "Muscle" "Oocyte" "Artery" "Ovary" "Uterus" "Placenta")

tis_sub_final=("Colon" "Duodenum" "Ileum" "Jejunum" "Macrophage" "Frontal_cortex" "Hypothalamus" "Pituitary" "Blastocyst" "Blastomere" "Morula")

tis_names=(${tis_main_final[*]} ${tis_sub_final[*]})

# for tis_i in {0..33}
for tis_i in $1
do
{
    tissue=${tis_names[tis_i]}
    mkdir "./output_splicing/$tissue"
    prefix="./output_splicing/$tissue/$tissue"

    # Input files 
    # Phenotypes
    expression_bed="../expression/splicing_bed/${tissue}.splicing_qqnorm.bed.gz"
    
    # if [ -f "./output_splicing/${tissue}/${tissue}.cis_qtl_pairs.1.txt.gz" ]; then
        # echo "Stop."
        # exit 0
    # fi

    # Covariates can be provided as a tab-delimited text file (covariates x samples) or dataframe (samples x covariates), with row and column headers.
    covariates_file="../pca_peer/covar_sqtl/${tissue}.covariates.txt"
    
    if [ ! -f ${covariates_file} ]; then
        echo "Can not find covariates_file."
        exit 1
    fi

    # Genotypes must be in PLINK format. Minor allele frequencies >=0.05; with the minor allele observed in at least 6 samples.
    plink_prefix_path="/BIGDATA2/scau_hzhang_1/USER/tengjy/Pig_GTEx_eQTL/eqtl_mapping/geno/${tissue}.filtered_maf0.05_mac6"

    # phenotype groups
    groups=../expression/splicing_groups/${tissue}.phenotype_groups.tsv

    # # Running tensorQTL >>> ${prefix} specifies the output file name.
    # # cis-QTL mapping: permutations  [CPU]
    # python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
        # --covariates ${covariates_file} \
        # --phenotype_groups ${groups} \
        # --mode cis \
        # --seed 9823
    
    # # eGene
    # Rscript eGene_detection.R ${prefix}.cis_qtl.txt.gz ${prefix}.cis_qtl_fdr0.05.txt.gz 0.05
    
    # independent eQTLs
    python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
        --covariates ${covariates_file} \
        --phenotype_groups ${groups} \
        --cis_output ${prefix}.cis_qtl_fdr0.05.txt.gz \
        --mode cis_independent \
        --seed 9823
    
    echo "Done."
}
done
wait

# for tis_i in {0..33}
# do
# {
    # while (( 1 ))
    # do
      # nl_yhq=$(yhq | wc -l)
      # if (( $nl_yhq > 145 ))
      # then
        # sleep 5
      # else
        # break
      # fi
    # done
    # yhbatch -p bigdata -J tjy${tis_names[tis_i]} pcg_eQTL_detection_permutations.sh $tis_i
# }
# done
# wait

