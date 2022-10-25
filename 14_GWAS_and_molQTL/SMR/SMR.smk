wkdir="SMR_analysis"
# 1. Prepare eQTL annotation
rule prepare_annot:
    input:
        dir_eqtl_allpairs="~/eqtl_mapping/allpairs_pcg/"
    output:
        ""
    shell:
        '''
        cd {wkdir}
        mkdir PreDATA
        mkdir ./PreDATA/FastQTL_format_pcg/
        mkdir ./PreDATA/BESD_file_pcg/
        
        tis_names=("Adipose" "Large_intestine" "Small_intestine" "Blood" "Fetal_thymus" "Lymph_node" "Milk" "Spleen" "Cartilage" "Synovial_membrane" "Brain" "Embryo" "Heart" "Kidney" "Liver" "Lung" "Testis" "Muscle" "Oocyte" "Artery" "Ovary" "Uterus" "Placenta" "Colon" "Duodenum" "Ileum" "Jejunum" "Macrophage" "Frontal_cortex" "Hypothalamus" "Pituitary" "Blastocyst" "Blastomere" "Morula")
        
        for tis_i in {{0..33}}
        do
        {{
            tissue=${{tis_names[tis_i]}}
            echo $tis_i
            echo $tissue
            zcat {input.dir_eqtl_allpairs}/${{tissue}}.eqtl_allpairs.txt.gz | awk -F"\t" 'FNR>1{{print $1"\t"$2"\t"$3"\t"$7"\t"$8}}' | pigz > ./PreDATA/FastQTL_format_pcg/${{tissue}}.eqtl_allpairs.txt.gz
        }}
        done
        wait
        
        for tis_i in {{0..33}}
        do
        {{
            tissue=${{tis_names[tis_i]}}
            echo $tis_i
            echo $tissue
            smr --eqtl-summary ./PreDATA/FastQTL_format_pcg/${{tissue}}.eqtl_allpairs.txt.gz --fastqtl-nominal-format --make-besd --out ./PreDATA/BESD_file_pcg/${{tissue}}.eqtl_allpairs
            rm ./PreDATA/FastQTL_format_pcg/${{tissue}}.eqtl_allpairs.txt.gz
        }}
        done
        wait
        
        rm -r /BIGDATA2/scau_hzhang_1/USER/tengjy/Pig_GTEx_GWAS_and_eQTLs/SMR/PreDATA/FastQTL_format_pcg
        '''

# 2. Update_BESD_file
rule update_besd_file:
    input:
        file_tss="Sus_scrofa.Sscrofa11.1.100.tss.gz",
        file_gtf="./PreDATA/Sus_scrofa.Sscrofa11.1.100.gene.gtf",
        dir_frq="~/snp_frq/", # allele frq from Plink
        dir_BESD_file="./PreDATA/BESD_file_pcg/"
    output:
        outdir="./PreDATA/expression/"
    shell:
        '''
        for f in `ls {input.dir_BESD_file}/*.eqtl_allpairs.* | grep -P "epi|esi"`
        do
            mv $f ${{f}}.old
        done
        Rscript Update_BESD_file.R {input.file_tss} {input.file_gtf} {input.dir_frq} {input.dir_BESD_file} 
        '''

# 3. Prepare GWAS summary
rule prepare_gwas_summary:
    input:
        dir_gwas_summary="~/GWAS_Summary_Release/"
    output:
        outdir="./PreDATA/gwas_summary/"
    shell:
        '''
        mkdir -p {output.outdir}
        files=(`ls {input.dir_gwas_summary}`)
        ntrait=${{#files[@]}} 
        for i in `seq 0 $[ntrait-1]`
        do
        {{
            trait=${{files[i]/.txt.gz/}}
            echo ${{i}}
            echo ${{trait}}
            zcat {input.dir_gwas_summary}/${{trait}}.txt.gz | awk -F"\t" '{{print $1" "$5" "$6" "$7" "$10" "$11" "$8" "$12}}' | \
            sed 's/variant_id effect_allele non_effect_allele frequency effect_size standard_error pvalue sample_size/SNP A1 A2 freq b se p N/g' \
            > {output.outdir}/${{trait}}.ma
        }} &
        done
        wait
        '''

# 4. Run SMR
rule run_smr:
    input:
        genobfile="" # specify the reference genotype file
    output:
        outdir="./output_SMR_pcg/"
    threads: 23
    shell:
        '''
        mkdir {output.outdir}
           
        tis_names=("Adipose" "Large_intestine" "Small_intestine" "Blood" "Fetal_thymus" "Lymph_node" "Milk" "Spleen" "Cartilage" "Synovial_membrane" "Brain" "Embryo" "Heart" "Kidney" "Liver" "Lung" "Testis" "Muscle" "Oocyte" "Artery" "Ovary" "Uterus" "Placenta" "Colon" "Duodenum" "Ileum" "Jejunum" "Macrophage" "Frontal_cortex" "Hypothalamus" "Pituitary" "Blastocyst" "Blastomere" "Morula")
        gwasfiles=(`ls ./PreDATA/gwas_summary/*.ma`)
        
        for tis_i in {{0..33}}
        do
            tissue=${{tis_names[tis_i]}}
            echo $tis_i
            echo $tissue
            mkdir {output.outdir}/${{tissue}}
            
            for gf_i in {{0..267}}
            do
            {{
                gf=${{gwasfiles[gf_i]}}
                trait=${{gf/"./PreDATA/gwas_summary/"/}}
                trait=${{trait/.ma/}}
                echo $gf_i
                echo $trait
                
                besdprefix=./PreDATA/BESD_file_pcg/${{tissue}}.eqtl_allpairs
                smr --bfile {input.genobfile} --gwas-summary ${{gf}} --beqtl-summary ${{besdprefix}} --peqtl-smr 1e-5 --diff-freq 0.9 \
                --diff-freq-prop 0.1 --extract-probe ./eGene/PCG/${tissue}.SMR.egenelist \
                --out {output.outdir}/${{tissue}}/${{tissue}}.pcg.${{trait}} --thread-num {threads}
                pigz {output.outdir}/${{tissue}}/${{tissue}}.pcg.${{trait}}.smr
            }} &
            done
            wait
        done
        wait
        '''
