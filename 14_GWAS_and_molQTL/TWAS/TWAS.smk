# 1. Prepare covariates file
rule prepare_covar:
    input:
        tissue="tissue_name"
    output:
        outdir="./PreDATA/covar/"
    shell:
        '''
        mkdir -p {output.outdir}
        covar={input.tissue}.covariates.txt # covariates file used in eQTL mapping
        awk -F"\t" '{{i=1;while(i <= NF){{col[i]=col[i] $i " ";i=i+1}}}} END {{i=1;while(i<=NF){{print col[i];i=i+1}}}}' ${{covar}} | \ 
        awk -F" " '{{print "0 "$1" "$0}}' | sed 's/0 pc1  /FID IID IID /g'| sed 's/[ ]*$//g' | awk -F" " '{{$3=null;print $0}}'| \
        sed 's/  / /g' > {output.outdir}/${{tissue}}.cov
        '''

# 2. Prepare expression matrix
rule prepare_expr:
    input:
        tissue="tissue_name"
    output:
        outdir="./PreDATA/expression/"
    shell:
        '''
        mkdir -p {output.outdir}
        bed={input.tissue}.bed.gz # bed file used in eQTL mapping
        zcat ${bed} | awk -F"\t" '{{$2=$1;$1=$4;$4=null;print $0}}' | sed 's/  / /g' | sed 's/gene_id #Chr end/GENE CHR GENE_COORD/g' | pigz > {output.outdir}/${tissue}.bed.gz
        '''

# 3. Prepare GWAS summary
rule prepare_gwas_summary:
    input:
        trait="ADG100" # specify the trait name used for analysis
    output:
        outdir="./PreDATA/gwas_summary/"
    shell:
        '''
        mkdir -p {output.outdir}
        gwasfile="{input.trait}.gwas.results" # specify the file name of gwas results
        echo "SNP N Z A1 A2" > {output.outdir}/{input.trait}.sumstats.gz
        zcat ${gwasfile} | awk -F"\t" '{{print $1" "$12" "$9" "$5" "$6}}' | pigz >> {output.outdir}/{input.trait}.sumstats.gz # you should specify the column index of "SNP N Z A1 A2" in your gwas result file
        '''

# 4. Prepare expression score
rule prepare_expr_score:
    input:
        expbfile="", # specify the genotype file used in eQTL mapping
        genobfile="", # specify the reference genotype file
        snplist="", # specify the file that includes the all SNP ID used for analysis (extracted from the reference genotype file)
        expr_dir="./PreDATA/expression/",
        cov_dir="./PreDATA/covar/",
        tissue="tissue_name"
    output:
        outdir="./output_MESC/"
    shell:
        '''
        for CHR in {{1..18}}; do
        python run_mesc.py --compute-expscore-indiv --plink-path plink --expression-matrix {output.expr_dir}/{input.tissue}.bed.gz \
         --exp-bfile {input.expbfile} --geno-bfile {input.genobfile} --chr ${{CHR}} --covariates {input.cov_dir}/{input.tissue}.cov \
         --keep {input.snplist} --gcta-path gcta64 --tmp ./tmp_{input.tissue} --out ./output_MESC/{input.tissue}
        done
        '''

# 5. Prepare LD Score
rule prepare_ldsc:
    input:
        genobfile = "",  # specify the reference genotype file
    output:
        outdir="./l2.ldscore/"
    shell:
        '''
        plink_prefix="" # specify the reference genotype file
        for CHR in {{1..18}}; do
        python ldsc.py --bfile {input.genobfile} --l2 --ld-wind-kb 500 --out {output.outdir}/RefPanel.chr${{CHR}}
        done
        '''

# 6. Run h2med
rule run_h2med:
    input:
        dir_gwas_summary="./PreDATA/gwas_summary/",
        dir_expr_score="./output_MESC/",
        dir_l2="./l2.ldscore/",
        dir_frq="./PreDATA/frq_RefPanel/",
        trait="trait_name",
        tissue="tissue_name"
    output:
        outdir="./output_h2med/"
    shell:
        '''
        mkdir -p {output.outdir}
        python run_mesc.py --h2med {input.dir_gwas_summary}/{input.trait}.sumstats.gz --exp-chr {input.dir_expr_score}/{input.tissue} \
        --ref-ld-chr {input.dir_l2}/RefPanel.chr@ --w-ld-chr {input.dir_l2}/RefPanel.chr@ --frqfile-chr {input.dir_frq}/RefPanel.chr@ \
        --out {output.outdir}/{input.tissue}.{input.trait}
        '''

