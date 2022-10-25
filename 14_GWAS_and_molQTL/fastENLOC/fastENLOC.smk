# 1. Prepare eQTL annotation
rule prepare_cis_qtl_annotation:
    input:
        dapg_results="dapg_tissue_dir",
        geno_vcf="snp_vcf_file",
        tissue="tissue_name"
    output:
        "fastenloc.eqtl.annotation.vcf.gz"
    shell:
        "summarize_dap2enloc.pl -dir {input.dapg_results} -vcf {input.geno_vcf} -tissue {input.tissue} | gzip - > {output}"

# 2. Match LD block
rule match_ld_block:
    input:
        GWAS_dir="~/MetaGWAS_data",
        trait="trait_name",
        LD_map="~/LD_map/HB_Pig.blocks.total.txt"
    output:
        zval_dir="~/colocalization/Pig_GTEx_GWAS/zval",
        zvalgz="{trait}.zval.gz"
    shell:
        "Rscript ./match_ld_block.r {input.GWAS_dir} {input.trait} {input.LD_map} {output.zval_dir}"
        "cd {zval_dir}"
        "cat {input.trait}*.txt > {input.trait}.zval"
        "rm -f {input.trait}*.txt"
        "gzip -f {input.trait}.zval"

# 3. Prepare GWAS summary
rule prepare_gwas_summary:
    input:
        pip_dir="~/colocalization/Pig_GTEx_GWAS/pip",
        zval_dir="~/colocalization/Pig_GTEx_GWAS/zval",
        trait="trait"
    output:
        "{trait}.pip"
    shell:
        "torus -d {input.zval_dir}/{input.trait}.zval.gz --load_zval -dump_pip {input.pip_dir}/{output}"
        "gzip {input.pip_dir}/{output}"

# 4. Run fastENLOC
rule run_fastenloc:
    input:
        eqtl_annot="fastenloc.eqtl.annotation.vcf.gz",
        gwas_pip="gwas.pip.gz",
        tissue="tissue_name"
    output:
        prefix="prefix_name"
    params:
        shrinkage=1
    threads:
        10
    shell:
        "fastenloc -eqtl {input.eqtl_annot} -gwas {input.gwas_pip} -t {input.tissue} -thread {threads} -prefix {output.prefix} -s {params.shrinkage}"
