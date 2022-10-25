#!/bin/bash

dir_nominal=${1}
dir_output=${2}
# tis_i=${3}
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

### 2. extract top pairs from nominal results for each tissue
nFile=`cat ${dir_output}/permutation_files.txt | wc -l`
# echo $nFile # 0..34
perm_files=(`find ${dir_nominal} -name "*.cis_qtl.txt.gz"`)
NAMEs=(${perm_files[@]/*\//})
NAMEs=(${NAMEs[@]/.cis_qtl.txt.gz/})
for tis_i in `seq 0 $[nFile-1]`
do
{
    name=${NAMEs[tis_i]}
    # nominal_files=(`ls ${dir_nominal}/${tissue}/${tissue}.cis_qtl_pairs.*.txt.gz`)
    nominal_files=(`find ${dir_nominal} -name "${name}.cis_qtl_pairs.*.txt.gz"`)
    rm -f ${dir_output}/${name}.nominal_files.txt
    for l in ${nominal_files[*]}
    do
    {
        echo ${l} >> ${dir_output}/${name}.nominal_files.txt
    }
    done
    # extract_pairs
    ./extract_pairs_tjy.py ${dir_output}/${name}.nominal_files.txt ${dir_output}/strong_pairs.combined_signifpairs.txt.gz ${name} -o ${dir_output}
    ### output file: *.extracted_pairs.txt.gz
    # rm -f ${dir_output}/${name}.nominal_files.txt
} &
done
wait

