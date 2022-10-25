# 1. run aFC
rule run_afc:
    input:
        tis="Adipose",  # Prefix of output file, like tissue name
        CHR=1  # 1..18
    output:
        outputdir="output_aFC"  # output dir
    params: boot=100
    shell:
        '''        
        geno=./genotypes/filtered/{input.tis}.filtered_maf0.05_mac6.geno.vcf.gz
        exp=./bed_TMM/{input.tis}.expr_tmm.bed.gz
        cov=./covar/{input.tis}.covariates.txt
        cis_qtl=./cis_qtl/{input.tis}.cis_qtl.txt.gz
        mkdir {output.outputdir}
        ./aFC.py --vcf ${{geno}} --pheno ${{exp}} --cov ${{cov}} --qtl ${{cis_qtl}} --chr {input.CHR} --log_xform 0 --output {output.outputdir}/{input.tis}.Chr{input.CHR}.log2aFC.txt --boot {params.boot}
        '''

# 2. Combine results from each chromosome
rule combine_results:
    input:
        tis="Adipose",
        dirRes="./output_log2aFC_TMM_top_cis_eqtl/"
    shell:
        '''
        awk 'FNR>1||NR==1' {input.dirRes}/{input.tis}.Chr*.log2aFC.txt | gzip > {input.dirRes}/{input.tis}.log2aFC.txt.gz
        rm {input.dirRes}/{input.tis}.Chr*.log2aFC.txt
        '''
