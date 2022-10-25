# wget http://ftp.ensembl.org/pub/release-100/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.100.gtf.gz
sample = ""  # sample id

# 1. Add read groups, sort, mark duplicates, and create index
rule task1:
    input:
        TMPDIR="",
        in_bam="{sample}-STARAligned.sortedByCoord.out.bam"
    output:
        "rg_added_{sample}-STARAligned.sortedByCoord.out.bam",
        "dedupped_{sample}-STARAligned.sortedByCoord.out.bam"
    threads: 10
    shell:
        '''        
        java -jar -Xmx15g -Djava.io.tmpdir={input.TMPDIR} picard.jar AddOrReplaceReadGroups I={input.in_bam} O=rg_added_{sample}-STARAligned.sortedByCoord.out.bam \
        RGID=4 RGLB=lib1 RGPL=illumina RGPU=run RGSM=20 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate
        
        java -jar -Xmx15g -Djava.io.tmpdir={input.TMPDIR} picard.jar MarkDuplicates I=rg_added_{sample}-STARAligned.sortedByCoord.out.bam O=dedupped_{sample}-STARAligned.sortedByCoord.out.bam \
        CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=marked_dup_metrics.txt
        '''

# 2. Split N trim and reassign mapping qualities
rule task2:
    input:
        reference="Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
    output:
        "{sample}.featureCounts_exon.txt"
    threads: 10
    shell:
        '''
        gatk --java-options "-Djava.io.tmpdir=/tmp" SplitNCigarReads --spark-runner LOCAL -R {input.reference} -L chr.list -I dedupped_{sample}-STARAligned.sortedByCoord.out.bam -O SplitNCigarReads_{sample}-STARAligned.sortedByCoord.out.bam
        '''

# 3. Base recalibration (BQSR)
rule task3:
    input:
        gtf="Sus_scrofa.Sscrofa11.1.100.gtf",
        dbSNP="sus_scrofa_BQSR.vcf.gz"
    output:
        expr_file="{sample}.enhancer.exp.txt"
    shell:
        '''
        gatk --java-options "-Djava.io.tmpdir=/tmp" BaseRecalibrator --spark-runner LOCAL -I SplitNCigarReads_{sample}-STARAligned.sortedByCoord.out.bam -R {input.reference} --known-sites {input.dbSNP} -O {sample}-recal.table
        gatk ApplyBQSR --spark-runner LOCAL -R {input.reference} -I SplitNCigarReads_{sample}-STARAligned.sortedByCoord.out.bam --bqsr-recal-file {sample}-recal.table -O BQSR_{sample}-STARAligned.sortedByCoord.out.bam
        '''

# 4. Run the haplotypecaller
rule run_haplotypecaller:
    input:
        interval="sus_scrofa_BQSR.bed"
    output:
        ""
    threads: 10
    shell:
        '''
        gatk --java-options "-Djava.io.tmpdir=/tmp" HaplotypeCaller --spark-runner LOCAL -R {input.reference} -I BQSR_{sample}-STARAligned.sortedByCoord.out.bam -O {sample}.vcf.gz -dbsnp {input.dbSNP} -L {input.interval} --dont-use-soft-clipped-bases --output-mode EMIT_ALL_CONFIDENT_SITES -stand-call-conf 0
        '''

# 5. Filter the snps to get the high quality ones
rule filter_snps:
    input:
        reference="Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
    output:
        out_filter="{sample}.filtered.vcf.gz"
    shell:
        '''
        gatk --java-options "-Djava.io.tmpdir=/tmp" VariantFiltration -R {input.reference} -V {sample}.vcf.gz -O {output.out_filter} \
         -window 35 -cluster 3 \
         --filter-name one \
         --filter-expression "FS>30.0" \
         --filter-name two \
         --filter-expression "QD<2.0"
        '''

# 6. Selected variants
rule select_variants:
    input:
        reference="Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
    output:
        out_filter="selected.{sample}.filtered.vcf.gz"
    shell:
        '''
        gatk SelectVariants -R {input.reference} -V {sample}.filtered.vcf.gz --select-type-to-include SNP --exclude-filtered -O {output.out_filter}
        '''
