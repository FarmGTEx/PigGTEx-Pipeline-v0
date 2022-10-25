# 1. clean data (trimmomatic v0.39)
rule clean_data:
    input:
        i="", # What is this? Please provide comment.
        id="sampleid"
    output:
        ""
    threads: 4
    shell:
        "trimmomatic PE -threads {threads} {i} {input.id}_R1.fastq.gz {input.id}_R2.fastq.gz {input.id}_R1.paired.fastq.gz {input.id}_R1.unpaired.fastq.gz {input.id}_R2.paired.fastq.gz {input.id}_R2.unpaired.fastq.gz MINLEN:50 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20"

# 2. bwa align and sorted (bwa v0.7.5a-r405, samtools v1.9)
rule bwa_align_and_sorted:
    input:
        id="sampleid"
    output:
        "{id}.sort.bam",
        "{id}.unsort.bam"
    threads: 4
    shell:
        '''
        bwa mem -t {threads} -M -R "@RG\tID:{input.id}\tSM:{input.id}" {input.ref} {input.id}_R1.paired.fastq.gz {input.id}_R2.paired.fastq.gz | samtools view -Sb -o {input.id}.unsort.bam -
        samtools sort -@ 4 -O bam -o {output[0]} {output[1]}
        '''

# 3. MarkDuplicates and generate index (GATK v4.1.4.1)
rule mark_index:
    input:
        id="sampleid"
    output:
        ""
    shell:
        '''
        gatk MarkDuplicates -I {input.id}.sort.bam -M {input.id}.marked_dup_metrics.txt -O {input.id}.rmdup.bam
        samtools index {input.id}.rmdup.bam
        '''

# 4. generate gvcf file for indi in chr region
rule generate_gvcf:
    input:
        fa="~/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa",
        id="sampleid",
        chr_region="1"  # 1..70 ?
    output:
        "{id}.{chr_region}.g.vcf.gz"
    shell:
        '''
        gatk HaplotypeCaller --emit-ref-confidence GVCF -R {input.fa} -I {input.id}.rmdup.bam -L {input.chr_region}.intervals -O {output}
        '''

# 5. CombineGVCFs for all samples and GenotypeGVCFs in single region
rule combineGVCFs:
    input:
        ref="ref",
        chr_region="chr_region"
    output:
        outdir="output"
    shell:
        '''
        gatk GenomicsDBImport -R {input.ref} --sample-name-map sample_map --genomicsdb-workspace-path database_{input.chr_region} --intervals {input.chr_region}.intervals --tmp-dir=/vol3/agis/likui_group/yinhongwei/gte_pig/gatk_tmp --reader-threads 20 
        # Optional, if your sample size > 1000.
        gatk GenotypeGVCFs -R {input.ref} -V gendb://database_{input.chr_region} -O {output.outdir}/{input.outnames}.{input.chr_region}.vcf.gz
        '''

# 6. combine genotype in single chr and SelectVariants
rule combine_chrs_SelectVariants:
    input:
        ref="ref",
        outnames="",# ??
        chr_region="1"  # 1..70 ?
    output:
        outdir="output"
    shell:
        '''
        gatk SelectVariants -R {input.ref} -V {output.outdir}/{input.outnames}.{input.chr_region}.vcf.gz -select-type SNP -O {output.outdir}/{input.outnames}.{input.chr_region}.raw.SNP.vcf.gz

        '''

# 7. mark the hard filter for snp
rule mark_hard_filter_snp:
    input:
        ref="ref",
        chr_region="chr_region"
    output:
        outdir="output"
    shell:
        '''
        gatk VariantFiltration -R {input.ref} -V {output.outdir}/503_chr${i}.raw.SNP.vcf.gz -filter \"QD < 2.0\" --filter-name \"QD\" -filter \"MQ < 40.0\" --filter-name \"MQ\" -filter \"FS > 60.0\" --filter-name \"FS\" -filter \"SOR > 3.0\" --filter-name \"SQR\" -filter \"MQRankSum < -12.5\" --filter-name \"MQRS\" -filter \"ReadPosRankSum < -8.0\" --filter-name \"RPRS\" -O {output.outdir}/{input.outnames}.{input.chr_region}.SNP.mark.vcf.gz
        '''

# 8. select snps which passed the hard filter and was the binary snp for autosome (bcftools v***)
rule filter_snps:
    input:
        chr="1",# 1..18
        outnames="",# ??
    output:
        outdir="output"
    shell:
        '''
        gatk MergeVcfs -I /vol3/agis/likui_group/yinhongwei/fst_5/chr_list/chr${i}.list -O {output.outdir}/{input.outnames}.{input.chr}.SNP.mark.vcf.gz
        bcftools view {output.outdir}/{input.outnames}.{input.chr}.SNP.mark.vcf.gz --min-af 0.01:minor -e 'F_MISSING>0.9' -m 2 -M 2 -v snps -f PASS -O z -o {output.outdir}/{input.outnames}.{input.chr}.filter.maf0.01.missing.binary.vcf.gz
        '''

# 9. phase (beagle v***)
rule phase:
    input:
        beagle="~/software/beagle/beagle.18May20.d20.jar",
        chr="1",# 1..18
        outnames=""  # ??
    output:
        "{outnames}.{chr}.filter.maf0.01.missing.binary.phased"
    params: seed=9823
    threads: 20
    shell:
        '''
        java -Xmx100g -jar {input.beagle} gt={input.outnames}.{input.chr}.filter.maf0.01.missing.binary.vcf.gz out={output} seed={seed} nthreads={threads}
        '''
