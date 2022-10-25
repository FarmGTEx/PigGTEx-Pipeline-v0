# 1. Data preparation for tensorQTL
rule pre_data:
    input:
        dir_deconvolution_results="./data_celltype/Multi-tissues_deconvolution",
        gsub_pattern="CIBERSORT(.)*", # to get tissue names
        dir_covariates="~/covariates_files/",
    output:
        dir_output="interaction_term_celltype"
    params: medianPro_cutoff=0.1
    shell:
        '''
        Rscript prepare_data_celltype.R {input.dir_deconvolution_results} {input.gsub_pattern} {input.dir_covariates} {params.medianPro_cutoff} {output.dir_output}
        '''

# 2. ieQTL mapping using tensorQTL
rule ieqtl_mapping:
    input:
        tissue="Blood",
        dir_exprssion="~/bed/",
        dir_covariates="~/covariates_files/",
        dir_interaction_term="~/interaction_term_celltype/",
        dir_geno="~/geno/"
    output:
        outputdir="output_celltype_maf0.1"
    params: maf_threshold_interaction=0.05
    shell:
        '''
        # Expresson levels (Inverse nominal transformed TMM)
        expression_bed="{input.dir_exprssion}/{input.tissue}.expr_tmm_inv.bed.gz"
        
        # Covariates can be provided as a tab-delimited text file (covariates x samples) or dataframe (samples x covariates), with row and column headers.
        covariates_file="{input.dir_covariates}/{input.tissue}.covariates.txt"
        
        # Inverse nominal transformed cell type enrichment/percentage
        interactions_file="{input.dir_interaction_term}/{input.tissue}.interactions_${pros[pro_i]}.txt"
        
        # Genotypes must be in PLINK format. Minor allele frequencies >=0.05; with the minor allele observed in at least 6 samples.
        plink_prefix_path="{input.dir_geno}/{input.tissue}.filtered_maf0.05_mac6"

        mkdir -p "./{output.outputdir}/{input.tissue}"
        
        files=(`ls ./interaction_term_celltype/{input.tissue}*.txt`)
        pros=(`echo ${{files[@]//*interactions_/}}`)
        pros=(`echo ${{pros[@]//.txt/}}`)
    
        for pro_i in `seq 0 $[${{#pros[@]}}-1]`
        do
        {{
            prefix="./{output.outputdir}/{input.tissue}/{input.tissue}.interactions_${{pros[pro_i]}}"
            python3 -m tensorqtl ${{plink_prefix_path}} ${{expression_bed}} ${{prefix}} \
                --covariates ${{covariates_file}} \
                --interaction ${{interactions_file}} \
                --maf_threshold_interaction {params.maf_threshold_interaction} \
                --best_only \
                --mode cis_nominal
        }}
        done
        wait
        '''

# 3. Define ieGene
rule define_iegene:
    input:
        wkdir="~/ieqtls/",
        ieqtl_dir="output_celltype_maf0.1"
    params:
        fdr_cutoff=0.05
    shell:
        '''
        Rscript Define_ieGene.R {input.wkdir} {input.ieqtl_dir} {params.fdr_cutoff}
        '''