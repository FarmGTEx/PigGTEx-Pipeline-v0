
sample="" # sample id

# 1. Build reference (STAR v***)
rule build_reference:
    input:
        fa="Sus_scrofa.Sscrofa11.1.dna.toplevel.fa", # What is this? Please provide comment.
        gtf="Sus_scrofa.Sscrofa11.1.100.gtf",
        STAR_dir=""  # Path to STAR
    output:
        "" # Please provide the output files names from this step.
    threads: 20
    shell:
        '''
        {input.STAR_dir}/STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir star-geneome \
        --genomeFastaFiles {input.fa} \
        --sjdbGTFfile {input.gtf}
        '''

# 2. Quality control & Mapping to genome (trimmomatic v0.39)
rule qc_and_map:
    input:
        STAR_dir="", # Path to STAR
        trimmomatic_dir="", # Path to trimmomatic
        gtf="Sus_scrofa.Sscrofa11.1.100.gtf",
        bio="" # What is this?
    output:
        "{sample}-STARAligned.sortedByCoord.out.bam"
    threads: 5
    shell:
        '''
        cd {sample}/
        
        if [ -e *_2.fastq.gz ]  # Paired-end ???
        then
            java -jar {input.trimmomatic_dir}/trimmomatic-0.39.jar PE -phred33 *_1.fastq.gz *_2.fastq.gz {sample}_1.clean.fq.gz {sample}_1_unpaired.fastq.gz {sample}_2.clean.fq.gz {sample}_2_unpaired.fastq.gz -threads {threads} \
            ILLUMINACLIP:{input.trimmomatic_dir}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        
            ${input.STAR_dir}/STAR --runThreadN {threads} \
                --genomeDir star-genome  \
                --sjdbGTFfile {input.gtf} \
                --quantMode TranscriptomeSAM \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped Within \
                --readFilesCommand zcat \
                --outFilterMismatchNmax 3 \
                --readFilesIn {sample}_1.clean.fq.gz {sample}_2.clean.fq.gz \
                --outFileNamePrefix {sample}-STAR
        
        else    # Single-ended sequencing ???
            java -jar {input.trimmomatic_dir}/trimmomatic-0.39.jar SE -phred33 *.fastq.gz {input.bio}.clean.fq.gz -threads {threads} \
            ILLUMINACLIP:{input.trimmomatic_dir}/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        
            ${STAR_dir}/STAR --runThreadN {threads} \
                --genomeDir  star-genome  \
                --sjdbGTFfile {input.gtf} \
                --quantMode TranscriptomeSAM \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMunmapped Within \
                --readFilesCommand zcat \
                --outFilterMismatchNmax 3 \
                --readFilesIn {sample}.clean.fq.gz \
                --outFileNamePrefix {sample}-STAR
        fi
        '''
