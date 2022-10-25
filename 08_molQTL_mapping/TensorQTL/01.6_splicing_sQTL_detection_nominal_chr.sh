#!/bin/bash
# conda activate deepmd
cd /BIGDATA2/scau_hzhang_1/USER/tengjy/Pig_GTEx_splicing/sqtl_mapping
#watch -n 2 nvidia-smi
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

    for CHR in $2
    do
    {
        # Input files 
        # Phenotypes
        expression_bed="../expression/splicing_bed/${tissue}/${tissue}.Chr${CHR}.splicing_qqnorm.bed.gz"
        
        if [ -f "./output_splicing/${tissue}/${tissue}.cis_qtl_pairs.${CHR}.txt.gz" ]; then
            echo "Stop."
            exit 0
        fi

        # Covariates can be provided as a tab-delimited text file (covariates x samples) or dataframe (samples x covariates), with row and column headers.
        covariates_file="../pca_peer/covar_sqtl/${tissue}.covariates.txt"

        if [ ! -f ${covariates_file} ]; then
            echo "Can not find covariates_file."
            exit 1
        fi

        # Genotypes must be in PLINK format. Minor allele frequencies >=0.05; with the minor allele observed in at least 6 samples.
        plink_prefix_path="/BIGDATA2/scau_hzhang_1/USER/tengjy/Pig_GTEx_eQTL/eqtl_mapping/geno/${tissue}.filtered_maf0.05_mac6"
        
        # phenotype groups
        groups=../expression/splicing_groups/${tissue}/${tissue}.Chr${CHR}.phenotype_groups.tsv

        # Running tensorQTL >>> ${prefix} specifies the output file name.
        # cis-QTL mapping: cis_nominal  [CPU]
        python3 -m tensorqtl ${plink_prefix_path} ${expression_bed} ${prefix} \
            --covariates ${covariates_file} \
            --phenotype_groups ${groups} \
            --mode cis_nominal
        # Convert .parquet to .txt.gz
        input="./output_splicing/${tissue}/${tissue}.cis_qtl_pairs.${CHR}.parquet"
        output="./output_splicing/${tissue}/${tissue}.cis_qtl_pairs.${CHR}.txt.gz"
        python3 parquet2txt.py $input $output
        if [ -f ${output} ]; then
            rm ${input}
        fi
    }
    done
    wait
    echo "Done."
}
done
wait

