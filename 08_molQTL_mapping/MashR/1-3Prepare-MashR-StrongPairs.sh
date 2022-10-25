#!/bin/bash

dir_nominal=${1}
dir_output=${2}

# tis_main_final=("Adipose" "Large_intestine" "Small_intestine" "Blood" "Fetal_thymus" "Lymph_node" "Milk" "Spleen" "Cartilage" "Synovial_membrane" "Brain" "Embryo" "Heart" "Kidney" "Liver" "Lung" "Testis" "Muscle" "Oocyte" "Artery" "Ovary" "Uterus" "Placenta")
# tis_sub_final=("Colon" "Duodenum" "Ileum" "Jejunum" "Macrophage" "Frontal_cortex" "Hypothalamus" "Pituitary" "Blastocyst" "Blastomere" "Morula")
# tis_names=(${tis_main_final[*]} ${tis_sub_final[*]})

mkdir ${dir_output}

############################################################################################
#   1: combine_signif_pairs_tjy.py   2: extract_pairs_tjy.py   3: mashr_prepare_input.py   #
#   1: extract info of SNPs across all studies                                             #
#   2: prepare files using info from Step 1                                                #
#   3: prepare input file by combining files from Step 2                                   #
############################################################################################

### 3. prepare strong SNP-gene pairs for MashR
strong_pairs_files=(`ls ${dir_output}/*.extracted_pairs.txt.gz`)
rm -f ${dir_output}/strong_pairs_files.txt
for l in ${strong_pairs_files[*]}
do
{
    echo ${l} >> ${dir_output}/strong_pairs_files.txt
}
done
# MashR format file (z-score)
./mashr_prepare_input.py ${dir_output}/strong_pairs_files.txt strong_pairs -o ${dir_output} --only_zscore
zcat ${dir_output}/strong_pairs.MashR_input.txt.gz | sed -e 's/_zval//g' | gzip > ${dir_output}/strong_pairs.temp.txt.gz
rm -f ${dir_output}/strong_pairs.MashR_input.txt.gz
mv ${dir_output}/strong_pairs.temp.txt.gz ${dir_output}/strong_pairs.MashR_input.txt.gz



