
dir_nominal = "../output_eqtl_breeds"
# 1. Prepare Strong pairs
rule pre_strong_pairs:
    output:
        dir_output_Strong="./breeds_eqtl_summary_Strong"
    shell:
        '''        
        # 1.1
        ./1-1Prepare-MashR-StrongPairs.sh {dir_nominal} {output.dir_output_Strong}
        # 1.2
        ./1-2Prepare-MashR-StrongPairs.sh {dir_nominal} {output.dir_output_Strong}
        # 1.3
        ./1-3Prepare-MashR-StrongPairs.sh {dir_nominal} {output.dir_output_Strong}
        '''

# 2. Prepare Random pairs
rule pre_random_pairs:
    output:
        dir_output_Random="./breeds_eqtl_summary_Random"
    params:
        subset_size=1000000
    shell:
        '''        
        # 2.1
        ./2-1Prepare-MashR-RandomPairs.sh {dir_nominal} {output.dir_output_Random} 
        # 2.2
        perm_files=(`find {dir_nominal} -name "*.cis_qtl.txt.gz"`)
        echo ${{#perm_files[@]}}   #0..34
        nFile=${{#perm_files[@]}}
        for i in `seq 0 $[nFile-1]`
        do
        {{
            for CHR in `seq 1 18`
            do
                ./2-2Prepare-MashR-RandomPairs.sh {dir_nominal} {output.dir_output_Random} ${{i}} ${{CHR}}
            done
        }}
        done
        wait
        
        # 2.3
        for CHR in `seq 1 18`
        do
            ./2-3Prepare-MashR-RandomPairs.sh {output.dir_output_Random} ${{CHR}}
        done
        wait
        
        # 2.4 
        ./MashR-random_subset.R {output.dir_output_Random} {params.subset_size}
        '''

# 3. Run MashR
rule run_mashr:
    input:
        prefix="name",
        dir_output_Strong="./breeds_eqtl_summary_Strong",
        dir_output_Random = "./breeds_eqtl_summary_Random"
    params:
        subset_size=1000000
    script:
        "./run_MashR.R {input.dir_output_Strong}/strong_pairs.MashR_input.txt.gz {input.dir_output_Random}/MashR.random_subset_{params.subset_size}.RDS 0 ./output_{input.prefix}_top_paris"
