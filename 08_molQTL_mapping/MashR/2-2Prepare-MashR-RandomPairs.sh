#!/bin/bash

dir_nominal=${1}
dir_output=${2}
tis_i=${3}
CHR=${4}

mkdir ${dir_output}

############################################################################################
#   1: combine_signif_pairs_tjy.py   2: extract_pairs_tjy.py   3: mashr_prepare_input.py   #
#   1: extract info of SNPs across all studies                                             #
#   2: prepare files using info from Step 1                                                #
#   3: prepare input file by combining files from Step 2                                   #
############################################################################################

### 2. extract all nominal pairs from nominal results for each tissue
# permutation results
perm_files=(`find ${dir_nominal} -name "*.cis_qtl.txt.gz"`)
NAMEs=(${perm_files[@]/*\//})
NAMEs=(${NAMEs[@]/.cis_qtl.txt.gz/})

name=${NAMEs[tis_i]}
nominal_files=(`find ${dir_nominal} -name "${name}.cis_qtl_pairs.${CHR}.txt.gz"`)
rm -f ${dir_output}/${name}_Chr${CHR}.nominal_files.txt
for l in ${nominal_files[*]}
do
{
    echo ${l} >> ${dir_output}/${name}_Chr${CHR}.nominal_files.txt
}
done
# extract_pairs
./extract_pairs_tjy.py ${dir_output}/${name}_Chr${CHR}.nominal_files.txt ${dir_output}/nominal_pairs.combined_signifpairs.txt.gz ${name}_Chr${CHR}.nominal_pairs --chrom ${CHR} -o ${dir_output}
# > output file: *_nominal_pairs.extracted_pairs.txt.gz
rm -f ${dir_output}/${name}_Chr${CHR}.nominal_files.txt

