# wget http://ftp.ensembl.org/pub/release-100/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.100.gtf.gz
sample = ""  # sample id

# 1. Specify gene biotypes to filter from the GTF file
rule prepare_reference:
    shell:
        '''        
        cellranger mkgtf Sus_scrofa.Sscrofa11.1.100.gtf Sus_scrofa.Sscrofa11.1.100.filtered.gtf \
            --attribute=gene_biotype:protein_coding \
            --attribute=gene_biotype:lincRNA \
            --attribute=gene_biotype:antisense \
            --attribute=gene_biotype:IG_LV_gene \
            --attribute=gene_biotype:IG_V_gene \
            --attribute=gene_biotype:IG_V_pseudogene \
            --attribute=gene_biotype:IG_D_gene \
            --attribute=gene_biotype:IG_J_gene \
            --attribute=gene_biotype:IG_J_pseudogene \
            --attribute=gene_biotype:IG_C_gene \
            --attribute=gene_biotype:IG_C_pseudogene \
            --attribute=gene_biotype:TR_V_gene \
            --attribute=gene_biotype:TR_V_pseudogene \
            --attribute=gene_biotype:TR_D_gene \
            --attribute=gene_biotype:TR_J_gene \
            --attribute=gene_biotype:TR_J_pseudogene \
            --attribute=gene_biotype:TR_C_gene
        '''

# 2. Create a Reference Package
rule create_ref:
    shell:
        '''
        cellranger mkref --genome=Sus_scrofa.Sscrofa11.1.100 \
            --fasta=Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
            --genes=Sus_scrofa.Sscrofa11.1.100.filtered.gtf \
            --ref-version=11.1.100
        '''

# 3. Mapping, filtering, barcode counting, and UMI counting
rule mapping:
    shell:
        '''
        cellranger count --id={sample} \
            --fastqs={sample} \
            --sample={sample} \
            --transcriptome=Sus_scrofa.Sscrofa11.1.100.cr3
        '''

# 4. Downstream analysis using Seurat
rule run_haplotypecaller:
    input:
        ""
    output:
        ""
    script:
        "Downstream_analysis.R"
