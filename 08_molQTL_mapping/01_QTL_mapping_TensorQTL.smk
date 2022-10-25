tissue="Adipose"
# 1. Data preparation for TensorQTL
rule pre_data:
    input:
        file_Counts="{tissue}_Counts.txt.gz", # Row is gene, column is sample; rowname is gene id, colname is sample id
        file_TPM="{tissue}_TPM.txt.gz", # Row is gene, column is sample; rowname is gene id, colname is sample id
        tss_annot_file="Sus_scrofa.Sscrofa11.1.100.tss.gz",
        vcf_file="{tissue}.geno.vcf.gz"  # from step 3.2 in first part.
    output:
        "{tissue}.covariates.txt"
    shell:
        '''        
        # Generate expression and covariates file
        Rscript prepare_data_for_tensorQTL.R {input.file_Counts} {input.file_TPM} {input.tss_annot_file} {input.vcf_file} {tissue} 
        '''

# 2. QTL mapping (eQTL, eeQTL, lncQTL, and enQTL) using tensorQTL
rule qtl_mapping:
    input:
        prefix="",
        expression_bed="./bed/{tissue}.expr_tmm_inv.bed.gz", # Phenotypes
        covariates_file="./covar/{tissue}.covariates.txt", # Covariates
        plink_prefix_path="{tissue}.geno" # Genotypes
    output:
        outputdir="output_eqtls"
    shell:
        '''
        mkdir {output.outputdir}
        mkdir "{output.outputdir}/{input.tissue}"
        prefix="{output.outputdir}/{input.tissue}/{input.tissue}"  # prefix specifies the output file name.
        
        # 2.2. permutation cis-eQTL mapping
        python3 -m tensorqtl {input.plink_prefix_path} {input.expression_bed} {input.prefix} \
            --covariates {input.covariates_file} \
            --mode cis \
            --seed 9823
            
        # Define eGene (FDR 0.05)
        Rscript eGene_detection.R {input.prefix}.cis_qtl.txt.gz {input.prefix}.cis_qtl_fdr0.05.txt.gz 0.05
    
        
        # 2.3. nominal cis-eQTL mapping
        python3 -m tensorqtl {input.plink_prefix_path} {input.expression_bed} {input.prefix} \
            --covariates {input.covariates_file} \
            --mode cis_nominal
        
        # 2.4. Convert .parquet to .txt.gz
        for CHR in `seq 1 18`
        do
        {{
            input="{input.prefix}.cis_qtl_pairs.${{CHR}}.parquet"
            output="{input.prefix}.cis_qtl_pairs.${{CHR}}.txt.gz"
            python3 parquet2txt.py $input $output
            if [ -f "${{output}}" ]; then
                rm ${{input}}
            fi
        }} &
        done
        
        # 2.5. Conditionally independent cis-eQTL mapping
        python3 -m tensorqtl {input.plink_prefix_path} {input.expression_bed} {input.prefix} \
            --covariates {input.covariates_file} \
            --cis_output {input.prefix}.cis_qtl.txt.gz \
            --mode cis_independent
        wait
        
        ### Done.
        '''

# sQTL...
