src="~/software/QTLEnrich-master/src"
wkdir="QTLEnrich_analysis"
# 1. Combine all variant-gene pairs from eQTL results
rule combine_all_pairs:
    input:
        dir_eqtl_chrs="~/eqtl_mapping/output_pcg/"
    output:
        dir_all_pairs="./PreDATA/allpairs_pcg/"
    shell:
        '''
        cd {wkdir}
        mkdir -p {output.dir_all_pairs}
        tis_names=("Adipose" "Large_intestine" "Small_intestine" "Blood" "Fetal_thymus" "Lymph_node" "Milk" "Spleen" "Cartilage" "Synovial_membrane" "Brain" "Embryo" "Heart" "Kidney" "Liver" "Lung" "Testis" "Muscle" "Oocyte" "Artery" "Ovary" "Uterus" "Placenta" "Colon" "Duodenum" "Ileum" "Jejunum" "Macrophage" "Frontal_cortex" "Hypothalamus" "Pituitary" "Blastocyst" "Blastomere" "Morula")
        
        for tis_i in {{0..33}}
        do
        {{
            tissue=${{tis_names[tis_i]}}
            echo $tissue
            cp {input.dir_eqtl_chrs}/${{tissue}}/${{tissue}}.cis_qtl_pairs.1.txt.gz {output.dir_all_pairs}/${{tissue}}.eqtl_allpairs.txt.gz
            for CHR in {{2..18}}
            do
                zcat {input.dir_eqtl_chrs}/${{tissue}}/${{tissue}}.cis_qtl_pairs.${{CHR}}.txt.gz | sed '1d' | pigz >> {output.dir_all_pairs}/${{tissue}}.eqtl_allpairs.txt.gz
            done
            wait
        }} &
        done
        wait
        '''

# 2. Prepare input files for QTLEnrich
rule prepare_variant_list:
    input:
        dir_eqtl_chrs="~/eqtl_mapping/output_pcg/",
        dir_geno_vcf="~/eqtl_mapping/geno_vcf/"
    output:
        outdir="./PreDATA/covar/"
    shell:
        '''
        mkdir -p ./data/eQTL_top
        mkdir -p ./PreDATA/
        
        # 2.1. Copy top eQTL files to specific fold
        cp {input.dir_eqtl_chrs}/*/*.cis_qtl_fdr0.05.txt.gz ./data/eQTL_top
        
        # 2.2. Generate file with unique list of variant IDs of all variants tested in the QTL analysis:
        echo "variant_id" > ./PreDATA/allsnps.txt
        for CHR in {{1..18}}
        do
            bcftools query -f "%CHROM\_%POS\_%REF\_%ALT\n" ./{input.dir_geno_vcf}/SNP.chr${{CHR}}.imputed.vcf.gz >> ./PreDATA/all_variants.txt
        done
        wait
        
        # 2.3. Generate list of significant QTLs (e.g., eQTLs and sQTLs if both computed)
        python3 {src}/extract_unique_variants.py --directory ../../Pig_GTEx_pcg/eqtl_mapping/sig_cis_qtl_pairs_pcg --file_name .cis_qtl_pairs.significant.txt.gz --output_file ./PreDATA/eQTL_significant_variants.txt
        
        # 2.4. Generate list of null variants:
        python3 {src}/create_null_variants.py --all_variants_file ./PreDATA/all_variants.txt --significant_variants_file ./PreDATA/eQTL_significant_variants.txt --output_file ./PreDATA/eQTL_null_variants.txt
        
        # 2.5 Generate the null variants table: (only use the unique variant_id)
        python3 {src}/Generate_Null_Table.py --QTL_Directory ../../Pig_GTEx_pcg/eqtl_mapping/allpairs_pcg/ --File_Extension .eqtl_allpairs.txt.gz --Null_Variants_List ./PreDATA/eQTL_null_variants.txt --Output_File ./PreDATA/eQTL_null_variants_table.txt.gz
        '''

# 3. Prepare LD proxy variants
rule prepare_ld_proxy_variants:
    input:
        dir_geno_ref="~/Pig_GTEx/ref_WGS_panel_chr/"
    output:
        "./PreDATA/number_ld_proxy_per_variant.txt"
    shell:
        '''        
        # 3.1. Compute LD proxy variants using PLINLK 1.90
        mkdir -p ./PreDATA/PGRP_genotype
        for CHR in {{1..18}}
        do
        plink --bfile {input.dir_geno_ref}/PGRP.chr${{CHR}} --extract ./PreDATA/all_variants.txt --make-bed --out ./PreDATA/PGRP_genotype/PGRP_extracted_variants.chr${CHR} &
        done
        wait
        
        mkdir -p ./PreDATA/PGRP_LD
        for CHR in {{1..18}}
        do
        echo $CHR
        plink --bfile ./PreDATA/PGRP_genotype/PGRP_extracted_variants.chr${CHR} --r2 --ld-snp-list ./PreDATA/all_variants.txt --ld-window 99999 --ld-window-kb 1000 --ld-window-r2 0.5 --out ./PreDATA/PGRP_LD/PGRP_LD_0.5.chr${CHR}
        done
        wait
        
        # 3.2. Parse Plink LD proxy output file into a table that contains number of LD proxy variants per variant that is inputted into Generate_Confounders_Table.py:
        python3 {src}/create_ld_table.py --plink_output_directory ./PreDATA/PGRP_LD/ --output_file ./PreDATA/number_ld_proxy_per_variant.txt

        '''

# 4. Prepare gene annotation
rule prepare_gtf:
    input:
        gtf="~/Pig_GTEx/Sus_scrofa.Sscrofa11.1.100.gtf.gz"
    output:
        "./PreDATA/Sus_scrofa.Sscrofa11.1.100.gene.gtf"
    shell:
        '''
        # 4. Processing the GENCODE file
        cp {input.gtf} ./PreDATA/Sus_scrofa.Sscrofa11.1.100.gtf.gz
        gunzip ./PreDATA/Sus_scrofa.Sscrofa11.1.100.gtf.gz
        
        echo -e "chr\tfeature\tstart\tend\tstrand\tgene_id\tgene_type" > ./PreDATA/Sus_scrofa.Sscrofa11.1.100.gene.gtf
        awk -F"[\t\";]" '{{if($3=="gene")print $1"\t"$3"\t"$4"\t"$5"\t"$7"\t"$10"\t"$22}}' ./PreDATA/Sus_scrofa.Sscrofa11.1.100.gtf >> ./PreDATA/Sus_scrofa.Sscrofa11.1.100.gene.gtf
        '''

# 5. Generating the Confounders Table
rule generate_confounders_table:
    input:
        dir_all_pairs="./PreDATA/allpairs_pcg/"
    output:
        "./PreDATA/eQTL_output_confounder.txt.gz"
    shell:
        '''
        # 5. Generating the Confounders Table
        python3 {src}/Generate_Confounders_Table.py --QTL_Directory {input.dir_all_pairs} --File_Extension .eqtl_allpairs.txt.gz --Variants_List ./PreDATA/all_variants.txt --LD_Proxy ./PreDATA/number_ld_proxy_per_variant.txt --GENCODE_File ./PreDATA/Sus_scrofa.Sscrofa11.1.100.gene.gtf --Output_File ./PreDATA/eQTL_output_confounder.txt.gz
        '''

# 6. Processing the GWAS input file
rule process_gwas_input:
    shell:
        '''                
        # 6. Processing the GWAS input file
        mkdir PreDATA
        mkdir ./PreDATA/gwas_summary
        
        gwasfiles=(`ls ~/GWAS_Summary/*.txt.gz`)
        ntrait=${{#gwasfiles[@]}}
        echo ${{#gwasfiles[@]}} # 267
        for t_i in `seq 0 $[ntrait-1]`
        do
        {{
            gf=${{gwasfiles[t_i]}}
            trait=${{gf/*Overlap_GWAS_eQTL\//}}
            trait=${{trait/.txt.gz/}}
            echo ${{trait}}
            echo -e "variant\tchr\tpos\tnon_effect_allele\teffect_allele\tgwas_p_value" | gzip > ./PreDATA/gwas_summary/${{trait}}.txt.gz
            zcat ${{gf}} | awk -F"\t" '{{print $1"\t"$3"\t"$4"\t"$6"\t"$5"\t"$8}}' | sed '1d' | gzip >> ./PreDATA/gwas_summary/${{trait}}.txt.gz
        }} &
        done
        wait
        '''

# 7. Run QTLEnrich for each trait
rule run_qtlenrich:
    input:
        ntrait=268
    output:
        outdir="output_QTLEnrich_eQTL.gwasp0.05"
    shell:
        '''
        mkdir -p {output.outdir}
        files=(`ls ./PreDATA/gwas_summary/`)
        echo ${{#files[@]}} 
        for i in `seq 0 $[ntrait-1]` 
        do
        {{
            trait=${{files[i]/.txt.gz/}}
            echo ${{i}}
            echo ${{trait}}
            file="./PreDATA/gwas_summary/${{trait}}.txt.gz"
            QTL_Directory='../QTLEnrich/data/eQTL_top/'
            File_Extension='.cis_qtl_fdr0.05.txt.gz' # suffix of eGenes file in GTEx v8
            QTL_Type='best_eqtl'
            Confounders_Table='./PreDATA/eQTL_output_confounder.txt.gz'
            Null_Table='./PreDATA/eQTL_null_variants_table.txt.gz'
            gencode_file='./PreDATA/Sus_scrofa.Sscrofa11.1.100.gene.gtf'
            exp_label="eQTL_${{trait}}"
            output_directory="{output.outdir}/QTLEnrich_eQTL_${{trait}}/"
            
            mkdir -p $output_directory
            
            #purge and load modules as necessary        
            python3 {src}/QTLEnrichV2.py --gwas_file $file \
                             --qtl_directory $QTL_Directory \
                             --file_name $File_Extension \
                             --qtl_type $QTL_Type \
                             --confounders_table $Confounders_Table \
                             --null_table $Null_Table \
                             --gencode_file $gencode_file \
                             --gwas_p_value 0.05 \
                             --exp_label $exp_label \
                             --output_directory ${{output_directory}} \
                             --GeneEnrich_input
            echo ">>> Done!"
        }}
        done
        wait
        '''
