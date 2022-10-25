# 1. Genotype imputation with beagle
rule run_imputation:
    input:
        beagle="./beagle.18May20.d20.jar",  # beagle software
        CHR=1  # Chromosome 1~18
    output:
        out="./imputed_genotype/PigGTEx.chr{CHR}.imputed"
    threads: 20
    params:
        seed=9823,
        ne=1000
    shell:
        '''        
        java -Xmx10g -jar $beagle \
            gt=./unimputed_genotype/PigGTEx.chr{CHR}.filtered.vcf.gz \
            ref=./ref_WGS_panel_chr/PigGTEx.ReferencePanel.chr{CHR}.vcf.gz \ # Reference panel
            out={output.out} \
            impute=true \
            ap=true \
            gp=true \
            seed={params.seed} \
            ne={params.ne} \
            nthreads={threads}
        bcftools index -f -t {output.out}.vcf.gz  # Build index
        '''

# 2. Quality control
rule qc:
    input:
        CHR=1  # Chromosome 1~18
    output:
        out="./filtered_genotype/PigGTEx.chr{CHR}.filtered.vcf.gz"
    shell:
        '''
        # Keep snps with DR2>=0.85 and minor allele frequency>=0.05
        bcftools view -q 0.05:minor -i "DR2>=0.85" ./imputed_genotype/PigGTEx.chr{CHR}.imputed.vcf.gz -Oz > {output.out}
        '''
