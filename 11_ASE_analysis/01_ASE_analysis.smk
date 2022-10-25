# 1. Mask reference genome for filtering out mapping bias
rule task1:
    input:
        CHR=1,# 1..18
        fa="Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
    output:
        ""  #
    shell:
        '''        
        bedtools maskfasta -fi {input.fa} -fo Ssc.chr1_18.mask.fa -bed RNA_SNPs.7008Samples.chr{input.CHR}.imputed.vcf.gz
        '''

# 2. Mapping by STAR
rule combine_results:
    input:
        sampleID="sample",
        genomeDir="~/ASE/06.maskfasta/chr",
        gtf="Sus_scrofa.Sscrofa11.104.chr1_18.gtf"
    output:
        ""
    threads: 5
    shell:
        '''
        #After filtering out mapping bias,this codes is for mapping by Star and detecting ASE by phASER-software 
        STAR --runThreadN {threads} \
            --genomeDir {input.genomeDir} \ # the path of mask reference genome
            --sjdbGTFfile {input.gtf} \ 
            --quantMode TranscriptomeSAM \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --readFilesCommand zcat \
            --outFilterMismatchNmax 3 \
            --readFilesIn {input.sampleID}_1.clean.fq.gz {input.sampleID}_2.clean.fq.gz \ #the path of clean data
            --outFileNamePrefix {input.sampleID}-STAR
        
        samtools index {input.sampleID}-STARAligned.sortedByCoord.out.bam
        samtools view -bq 255 {input.sampleID}-STARAligned.sortedByCoord.out.bam {input.sampleID}-sort.unique.bam
        samtools index {input.sampleID}-sort.unique.bam
        
        gatk AddOrReplaceReadGroups -I {input.sampleID}-sort.unique.bam -O rg_added_{input.sampleID}-STARsorted.unique.bam -RGID 4 -RGLB lib1 -RGPL illumina -RGPU run -RGSM 20 -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT -SORT_ORDER coordinate
        
        gatk MarkDuplicates -I rg_added_{input.sampleID}-STARsorted.unique.bam -O dedupped_{input.sampleID}-STARsorted.unique.bam -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT --READ_NAME_REGEX null -M dedupped_{input.sampleID}-STARsorted.unique.marked_dup_metrics.txt
        '''


# 3. Haplotype phasing
rule phasing:
    input:
        sampleID="sample"
    shell:
        '''
        #Performs haplotype phasing using read alignments in BAM format from RNA based assays.
        python2 ~/phaser-master/phaser/phaser.py \
            --vcf RNA_SNPs.7008Samples.chr1_18.imputed.vcf.gz \
            --bam dedupped_sampleID-STARsorted.unique.bam \
            --paired_end 1 \
            --mapq 255 \
            --baseq 10 \
            --sample {input.sampleID} \
            --threads 5 \
            --temp_dir tmp2 \
            --haplo_count_blacklist Sus_scrofa_chrAll_75_2_0.5.bed \
            --gw_phase_vcf 1 \
            --o sampleID \
            --output_read_ids 1
        '''

# 4. Produce gene level haplotype counts for allelic expression studies
rule gene_level_hap_counts:
    shell:
        '''
        # phASER Gene AE. Uses output from phASER to produce gene level haplotype counts for allelic expression studies. It does this by summing reads from both single variants and phASER haplotype blocks using their phase for each gene. 
        #Developed by [Stephane E. Castel](mailto:scastel@nygenome.org) in the [Lappalainen Lab](http://tllab.org) at the New York Genome Center and Columbia University Department of Systems Biology.

        python2 /phaser-master/phaser_gene_ae/phaser_gene_ae.py \
            --haplotypic_counts sampleID.haplotypic_counts.txt \
            --features Sus_scrofa.Sscrofa11.1.100.gene.phaser.bed \
            --o sampleID_phaser.gene_ae.txt
        '''

# 5. Aggregates gene-level haplotypic expression measurement files
rule aggr_gene_level_expr:
    shell:
        '''
        #phASER-pop
        #Aggregates gene-level haplotypic expression measurement files across samples to produce a single haplotypic expression matrix, where each row is a gene and each column is a sample
        python2 /phaser-master/phaser_pop/phaser_expr_matrix.py \
            --gene_ae_dir /path/all.phase_gene.ae \ # the path contains all sample's phaser.gene_ae.txt
            --features Sus_scrofa.Sscrofa11.1.100.gene.phaser.bed \
            --t 5 \
            --o phaser_expr.7008
        '''

# 6. Measure cis-eQTL effect size by ASE
rule measure_ase_afc:
    shell:
        '''
        tis=${tis_names[tis_i]}
        python phaser_cis_var.py --bed ./expression/haploCount_gw/${tis}.ase_expr_count.gw_phased.bed.gz --vcf ./phaser_vcf/${tis}.phaser.vcf.gz --pair ./data/ASE_aFC.test_pairs.txt --map ./data/ASE.sample_map.txt --o ./output_ASE_aFC/${tis}.aFC.txt --t 23
        '''
