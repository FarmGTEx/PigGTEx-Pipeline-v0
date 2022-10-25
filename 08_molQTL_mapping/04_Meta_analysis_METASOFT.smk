prefix = "name"  # prefix of output file

# 1. prepare SNP-gene pairs for Metasoft
rule pre_data:
    input:
        DIR_summary_Strong="~/breeds_eqtl_summary_Strong"
    output:
        dataDIR="./data/metasoft_data_{prefix}"
    shell:
        '''        
        ### 
        # Metasoft format file (beta and se)
        mkdir -p {output.dataDIR}
        
        # find Strong pairs file
        meta_files=(`ls {input.DIR_summary_Strong}/{prefix}.*.extracted_pairs.txt.gz`)
        rm -f {output.dataDIR}/meta_files.{prefix}.txt
        for l in ${{meta_files[*]}}
        do
            echo ${{l}} >> {output.dataDIR}/meta_files.{prefix}.txt
        done
        files_list={output.dataDIR}/meta_files.{prefix}.txt #nonimal file list the same as MashR
        
        ./metasoft_prepare_input_tjy.py ${{files_list}} {prefix} -o {output.dataDIR} 
        '''

# 2. Run METASOFT
rule run_metasoft:
    input:
        metasoft="~/software/Metasoft/Metasoft.jar",
        dataDIR="./data/metasoft_data_{prefix}"
    output:
        resultDIR="./result/output_meta_{prefix}"
    shell:
        '''
        mkdir -p {output.resultDIR}
        chunks_files=(`ls {input.dataDIR}/{prefix}*`)
        input=${{chunks_files[0]}}
        ./run_metasoft.py {input.metasoft} ${{input}} {prefix} -o {output.resultDIR}
        '''
