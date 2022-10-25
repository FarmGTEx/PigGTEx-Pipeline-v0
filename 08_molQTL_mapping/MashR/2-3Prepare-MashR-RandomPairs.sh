#!/bin/bash

dir_output=${1}
CHR=${2}

mkdir ${dir_output}

############################################################################################
#   1: combine_signif_pairs_tjy.py   2: extract_pairs_tjy.py   3: mashr_prepare_input.py   #
#   1: extract info of SNPs across all studies                                             #
#   2: prepare files using info from Step 1                                                #
#   3: prepare input file by combining files from Step 2                                   #
############################################################################################

### 3. prepare random SNP-gene pairs for MashR
nominal_pairs_files=(`ls ${dir_output}/*_Chr${CHR}.nominal_pairs.extracted_pairs.txt.gz`)
rm ${dir_output}/Chr${CHR}.nominal_pairs_list_file.txt
for l in ${nominal_pairs_files[*]}
do
{
    echo ${l} >> ${dir_output}/Chr${CHR}.nominal_pairs_list_file.txt
}
done
./mashr_prepare_input.py ${dir_output}/Chr${CHR}.nominal_pairs_list_file.txt nominal_pairs_Chr${CHR} -o ${dir_output} --only_zscore --dropna
rm -f ${dir_output}/Chr${CHR}.nominal_pairs_list_file.txt


