tissue="Adipose"
# 1. Make GRM
rule make_grm:
    input:
        plink_prefix="../eqtl_mapping/geno/{tissue}.filtered_maf0.05_mac6"
    output:
        out_prefix="./GRM/{tissue}"
    shell:
        '''
        mkdir ./GRM
        gcta64 --bfile {input.plink_prefix} --make-grm --out ./GRM/{tissue}
        gcta64 --grm ./GRM/{tissue} --make-bK-sparse 0.05 --out ./GRM/{tissue} 
        '''

# 2. Run fastGWA
rule run_fastgwa:
    input:
        CHR=1,  # 1..18
        TIS={tissue}
    shell:
        '''
        julia cis_eqtl_mapping_fastGWA.jl {input.CHR} {input.TIS} 
        '''

# 3. Combine results
rule combine_results:
    input:
        CHR=1,# 1..18
        TIS={tissue}
    shell:
        '''
        julia combine_eqtls_fastGWA.jl {input.CHR} {input.TIS}
        '''
