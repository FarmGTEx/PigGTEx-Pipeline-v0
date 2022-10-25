# wget http://ftp.ensembl.org/pub/release-100/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.100.gtf.gz
sample="" # sample id

# 1. Expression quantification for all genes (stringtie v***, featureCounts v***)
rule expr_gene:
    input:
        gtf="Sus_scrofa.Sscrofa11.1.100.gtf",
        STAR_dir=""  # Path to STAR
    output:
        "{sample}-STAR_STARgenome/{sample}.gtf",
        "{sample}.tsv",
        "{sample}.featureCounts.txt" # Please provide the output files names from this step.
    threads: 10
    shell:
        '''        
        # FPKM quantification
        stringtie -p 5 -e -B -G {input.gtf} -o ./{sample}-STAR_STARgenome/{sample}.gtf -A {sample}.tsv {sample}-STARAligned.sortedByCoord.out.bam
        # Count quantification (genes)
        featureCounts -T {threads} -p -t exon -g gene_id -a {input.gtf} -o {sample}.featureCounts.txt {sample}-STARAligned.sortedByCoord.out.bam
        '''

# 2. Expression quantification for exons (featureCounts v***)
rule expr_exon:
    input:
        STAR_dir="", # Path to STAR
        trimmomatic_dir="", # Path to trimmomatic
        gtf="Sus_scrofa.Sscrofa11.1.100.gtf",
        bio="" # What is this?
    output:
        "{sample}.featureCounts_exon.txt"
    threads: 10
    shell:
        '''
        # Count quantification (exons)
        featureCounts -T {threads} -p -t exon -g gene_id -f -a {input.gtf} -o {sample}.featureCounts_exon.txt {sample}-STARAligned.sortedByCoord.out.bam
        '''

# 3. Expression quantification for enhancers
rule expr_enhancer:
    input:
        gtf="Sus_scrofa.Sscrofa11.1.100.gtf"
    output:
        out_saf="enhancer.saf",
        expr_file="{sample}.enhancer.exp.txt"
    shell:
        '''
        # Step 1. Prepare saf file
        grep -v '#' {input.gtf} | awk -v OFS="\t" '{{print $1,$4,$5,$7}}' | sort -k1,1 -k2,2n | uniq > Sus_scrofa.Sscrofa11.1.100.sort.bed
        bedtools intersect -a enhancer.bed -b Sus_scrofa.Sscrofa11.1.100.sort.bed -wa -wb | awk -v OFS="\t" '{{print $1,$2,$3,$NF}}' | uniq | awk '{{print "enhancer&"$1":"$2":"$3"\t"$1"\t"$2"\t"$3"\t"$4}}' | sed '1i GeneID\tChr\tStart\tEnd\tStrand' > {output.out_saf}
        
        # Step 2. Count quantification
        featureCounts -T 10 -p -F SAF -a {output.out_saf} -o {sample}.enhancer.exp.txt {sample}-STARAligned.sortedByCoord.out.bam
        '''

# 4. Expression quantification for lncRNAs
rule expr_lncrna:
    input:
        ""
    output:
        out_saf="novel_lncrna.saf"
    threads: 10
    shell:
        '''
        # Step 1. Detection of novel lncRNAs (Elisabetta)
        # novel.feelnc_biotype.ok.gff
        
        # Step 2. Prepare the saf file
        grep 'feelnc_biotype "lncRNA"' novel.feelnc_biotype.ok.gff | grep -w 'transcript'> novel_lncrna.gff
        agat_convert_sp_gff2gtf.pl --gff novel_lncrna.gff -o novel_lncrna.gtf
        awk -F '\t' '{{print $1,$4,$5,$7,$NF}}' novel_lncrna.gff | awk -F ';' '{{print $1,$2}}' | sed 's/gene_id //g' | sed 's/ transcript_id //g' | sed 's/\"//g' | awk '{{print $5,"&",$6,$1,$2,$3,$4}}' | sed 's/\ & /\&/g' | sed 's/ /\t/g' | sed '1i GeneID\tChr\tStart\tEnd\tStrand'> {out_saf}
        
        # Step 3. Count quantification
        featureCounts -T {threads} -p -F SAF -a {out_saf} -o {sample}.novel_lnc.exp.txt {sample}-STARAligned.sortedByCoord.out.bam
        '''

# 5. Quantification for splicing events
rule expr_splicing:
    input:
        bam2junc="leafcutter/scripts/bam2junc.sh"
    output:
        ""
    shell:
        '''
        # Step 1. BAM to junc
        sh {bam2junc} {sample}-STARAligned.sortedByCoord.out.bam {sample}.junc
        
        # Step 2. Generating intron excision ratios with LeafCutter
        python3 cluster_prepare_fastqtl.py \
            juncfile.txt \
            Sus_scrofa.Sscrofa11.1.100_exon.bed \
            Sus_scrofa.Sscrofa11.1.100.gtf \
            sample_tissue \
            --min_clu_reads 30 --min_clu_ratio 0.001 --max_intron_len 500000 --num_pcs 15 \
            --leafcutter_dir leafcutter/ \
            --output_dir sample_tissue
        '''
