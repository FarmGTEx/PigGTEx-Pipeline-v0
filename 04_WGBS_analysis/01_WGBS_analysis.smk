# wget http://ftp.ensembl.org/pub/release-100/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.100.gtf.gz
sample = ""  # sample id

# 1. Prepare reference
rule prepare_reference:
    shell:
        '''        
        bismark_genome_preparation --bowtie2 pig_genome/
        '''

# 2. Quality control
rule qc:
    output:
        "{sample}_1.fastq.gz",
        "{sample}_2.fastq.gz"
    threads: 10
    shell:
        '''
        fastqc --noextract -t {threads} -f fastq -o ./fastqc *_1.fastq.gz
        fastqc --noextract -t {threads} -f fastq -o ./fastqc *_2.fastq.gz
        trim_galore --paired --fastqc --max_n 15 -o ./trim_reads {sample}_1.fastq.gz {sample}_2.fastq.gz
        '''

# 3. Mapping
rule mapping:
    threads: 2
    shell:
        '''
        bismark --multicore {threads} --bowtie2 --gzip -p 4 -N 0 -o bamfile pig_genome -1 {sample}_1_val_1.fq.gz -2 {sample}_2_val_2.fq.gz
        '''

# 4. Deduplicate the alignment BAM file
rule run_haplotypecaller:
    input:
        "{sample}_1_val_1_bismark_bt2_pe.bam"
    output:
        ""
    shell:
        '''
        deduplicate_bismark -p --bam {sample}_1_val_1_bismark_bt2_pe.bam
        '''

# 5. Extract context-dependent (CpG/CHG/CHH) methylation
rule filter_snps:
    input:
        ""
    output:
        ""
    threads: 8
    shell:
        '''
        bismark_methylation_extractor -p --ignore_r2 6 --multicore {threads} --bedgraph -o ../methylation --cytosine_report --genome_folder pig_genome *.deduplicated.bam
        '''

# 6. Prepare methylation levels
rule pre_meth:
    input:
        bio_id="" # What is this?
    output:
        ""
    shell:
        '''
        awk -F "[_.\t]" '{{print $i"_"$2}}' {sample}_1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.CpG_report.txt > {sample}_1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.CpG_report1.txt
        awk 'OFS="\t"{{if (($4+$5)>0) print $1,$2,$3,"CpG",$4/($4+$5),$4+$5}}' {input.bio_id}_1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.CpG_report1.txt | sort -k1,1 -k2,2n > {sample}.meth
        symmetric-cpgs {sample}.meth > {sample}.cpg.meth
        awk 'OFS="\t"{{if ($6>=5) print $1,$2,$2+1,$5}}' {sample}.cpg.meth > {sample}.meth.bed
        sed -i 's/^/chr/g' {sample}.meth.bed
        '''

# 7. Tissue-specific MR
rule tissue_spcific_mr:
    shell:
        '''
        awk -v OFS="\t" '{{print $1":"$2":"$3,$0}}' sample.meth.bed | sed '/^#/d' | awk '{{print $1"\t"$NF}}' | sed 's/:/\t/g' > sample_new.meth
        
        cp Sus_scrofa.chrAll.CpG.positions2.txt 0.matrix
        
        c=0
        for i in `ls *_new.meth`
        do
        b=$(($c+1))
        awk 'FNR==NR{{a[$1"\t"$2"\t"$3]=$4; next}}{{k=$1"\t"$2; print $0,a[$1"\t"$2"\t"$3]?a[$1"\t"$2"\t"$3]:"-"}}' OFS="\t" $i ${{c}}.matrix > ${{b}}.matrix
        rm ${c}.matrix
        c=$(($c+1))
        done
        
        SMART final.matrix -t DeNovoDMR -n tissue_specific_grouped -c case_control_matrix.txt -o . -MR 0.5 -AG 1.0 -MS 0.5 -ED 0.2 -SM 0.6 -CD 500 -CN 5 -SL 20 -PD 0.05 -PM 0.05 -AD 0.3
        '''
# 8. Allele-specific methylation
rule allele_specific_methylation:
    shell:
        '''
        # Step8.1 Convert bam file to mr (mapped fread) file
        to-mr -o sample.mr -L 10000 sample.deduplicated.bam -m bismark
        sort -k 1,1 -k 2,2n -k 3,3n sample.mr > sample.sort.mr
        
        # Step8.2 Epiread Format
        methstates -c pig_genome/ -o sample.epiread sample.sort.mr
        
        # Step8.3 Single-site ASM scoring
        allelicmeth -c pig_genome/ -o sample.allelic sample.epiread
        
        # Step8.4 Allelically methylated regions (AMRs)
        /amrfinder -o sample.amr -c pig_genome/ sample.epiread
        '''

# 9. SNP calling for Muscle
rule snp_calling:
    shell:
        '''
        # Step 9.1 Sort bam files
        samtools sort -o sample.deduplicated.sort.bam -m 3G -@ 4 sample.deduplicated.bam
        
        # Step 9.2 Index bam files
        samtools index -@ 4 sample.deduplicated.sort.bam sample.deduplicated.sort.bam.bai
        java -Xmx10g -Djava.io.tmpdir=$TMPDIR -jar picard.jar AddOrReplaceReadGroups I=samplededuplicated.sort.bam O=sample.deduplicated.sort.RG.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=run RGSM=20 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate
        
        # Step 9.3 Base quality recalibration SNP file must indexed by IndexFeatureFile gatk
        gatk BaseRecalibrator -R Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -I sample.deduplicated.sort.RG.bam --known-sites sus_scrofa_BQSR.vcf -O sample_recalFile_before.csv 
        
        # Step 9.4 Write recalibrated base quality score into BAM file
        gatk ApplyBQSR -R Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -I sample.deduplicated.sort.RG.bam --bqsr-recal-file sample_recalFile_before.csv -O sample.deduplicated.sort.RG.recal.bam
        
        # Step 9.5 Genotyping SNP file must indexed by IndexFeatureFile gatk
        java -Xmx10g -Djava.io.tmpdir=$TMPDIR -jar BisSNP-1.0.1.jar -R Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -T BisulfiteGenotyper -I sample.deduplicated.sort.RG.recal.bam -D sus_scrofa_BQSR7.vcf -vfn1 sample.dbsnp.raw.vcf -nt 4 -stand_call_conf 10 -mmq 30 -mbq 17 -out_modes EMIT_ALL_SITES -L sus_scrofa_BQSR.bed
        
        # Step 9.6 Sort the dbsnp.row.vcf 
        sed -i '/^$/d' sample.dbsnp.raw.vcf
        vcf-sort -c -p 4 sample.dbsnp.raw.vcf > sample.dbsnp.raw.vcf.sorted
        
        # Step 9.7 Filter fake SNPs--flter out SNPs with quality score less than 20, reads coverage more than 120, strand bias more than -0.02, quality score by depth less than 1.0, mapping quality zero reads fraction more than 0.1 and 2 SNPs within the same 20 bp window
        java -Xmx10g -Djava.io.tmpdir=$TMPDIR -jar BisSNP-1.0.1.jar -R Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -T VCFpostprocess -oldVcf sample.dbsnp.raw.vcf.sorted -newVcf sample.dbsnp.filtered.vcf -snpVcf sample.dbsnp.raw.vcf.sorted -qual 20 -windSizeForSNPfilter 10 -o sample.dbsnp.raw.filter.summary.txt
        '''