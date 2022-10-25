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

### 1. output all top SNP-gene pairs information from permutation results to a .txt file
# (colnames: phenotype_id,variant_id,chr,pos)
# file list of permutation results
rm -f ${dir_output}/permutation_files.txt
# permutation results
perm_files=(`find ${dir_nominal} -name "*.cis_qtl.txt.gz"`)
for l in ${perm_files[*]}
do
{
    echo ${l} >> ${dir_output}/permutation_files.txt
}
done
./combine_signif_pairs_tjy.py ${dir_output}/permutation_files.txt strong_pairs -o ${dir_output}
#> output file: strong_pairs.combined_signifpairs.txt.gz

