#**
# @description: cis-heritability estimation (unpublished)
# @author: Jinyan Teng
# @email: kingyan312@live.cn
# @date: 2021/7/13
#**

using Base.Threads
using Statistics
using HypothesisTests
using Random
using DataFrames
using StatsBase
using CSV
using Dates
using Distributed
using LinearAlgebra
using Mmap
using GZip
#==============================================================================#
#==============================================================================#
CHR=ARGS[1]
TIS=ARGS[2]
# Parameters setting ----------------------------------------------------------#
home_dir = "~/USER/tengjy/Pig_GTEx_eQTL/heritability-cis"
cd(home_dir)
# tissue=["Kidney"][1]  # name of tissue/study
tissue=TIS
expression_bed_file = string("../eqtl_mapping/bed/",tissue,".expr_tmm_inv.bed.gz")
geno_file_prefix = string("../eqtl_mapping/geno/",tissue,".filtered_maf0.05_mac6")
cov_file = string("../eqtl_mapping/covar_5or10pc10peer/",tissue,".covariates.txt")
snp_annot_file = string("../genotypes/filtered/snp_annot_maf0.05_mac6/",tissue,".snp.map")
# gene_annot_file = "../eqtl_mapping/Sus_scrofa.Sscrofa11.1.100.tss.gz"
out_dir=string("./output_covar_5or10pc10peer/", tissue, "/")
# window=1e6
# seed=2021
# -----------------------------------------------------------------------------#
# Load data from disk
io = GZip.open(expression_bed_file)
expression = CSV.File(io,header=true) |> DataFrame
close(io)
n_samples = ncol(expression) - 4
gene_annot_all = expression[:,1:4]
snp_annot_all = CSV.File(snp_annot_file,header=false) |> DataFrame
covar = CSV.File(cov_file,transpose=true) |> DataFrame
rename!(covar,:Column1=>:SampleID)
#
mkpath(out_dir)
# Pull covariates
outcov = hcat(repeat([0],n_samples),covar)
tmpdir=string(out_dir,"_temp/")
mkpath(tmpdir)
# for chrom in 1:18
for chrom in parse.(Int,CHR)
    # chrom = 1
    outcov_file = string(tmpdir,tissue,".chr",chrom,".qcovar")
    CSV.write(outcov_file,outcov,header=false,delim="\t")
    #
    gene_annot = gene_annot_all[gene_annot_all."#Chr" .== chrom,:]
    snp_annot = snp_annot_all[snp_annot_all.Column1 .== chrom,:]
    # exp_genes = names(expression)
    exp_genes = gene_annot.gene_id
    sampleid = names(expression)[5:end]
    n_genes = nrow(gene_annot)
    # seed = ifelse(isnan(seed), sample(1:2021,1)[1],seed)
    # Random.seed!(seed)
    h2_results = DataFrame(geneid=[],nidv=[],nsnps=[],window=[],typeofh2=[],h2=[],se=[],Vg=[],Ve=[],genename=[])
    h2_file = string(out_dir,tissue,".chr",chrom,".cis_h2.txt")
    CSV.write(h2_file,h2_results,delim="\t")
    #
    coveff_results = "geneid\twindow\tintercept\t"
    for x in names(outcov)[3:end]
        coveff_results = string(coveff_results,"\t",x)
    end
    coveff_results = string(coveff_results,"\n")
    coveff_file = string(out_dir,tissue,".chr",chrom,".cis_coveff.txt")
    write(coveff_file,coveff_results)
    #
    coveffse_results = "geneid\twindow\tintercept\t"
    for x in names(outcov)[3:end]
        coveffse_results = string(coveffse_results,"\t",x)
    end
    coveffse_results = string(coveffse_results,"\n")
    coveffse_file = string(out_dir,tissue,".chr",chrom,".cis_coveffse.txt")
    write(coveffse_file,coveffse_results)
    #
    # windows = Int.([1e6, 2e6, 3e6])
    windows = Int.([1e6])
    for window in windows
        # window = windows[1]
        for i in 1:n_genes
            # i = 1
            println("GENE: ",i,"/",n_genes)
            gene = exp_genes[i]
            # Reduce genotype data to only include SNPs within specified window of gene.
            geneinfo = gene_annot[gene_annot.gene_id .== gene,:]
            genename = ""
            start_pos = geneinfo.start[1] - window
            start_pos = ifelse(start_pos<0, 0, start_pos)
            end_pos = geneinfo.end[1] + window
            # Pull cis-SNP info
            cissnps = snp_annot[(snp_annot.Column2 .>= start_pos) .& (snp_annot.Column2 .<= end_pos),:]
            if nrow(cissnps) < 1
                # Need 1 or more cis-snps
                @warn string("Zero cis-snps for ", gene)
                h2_DF = DataFrame([:geneid=>gene,:nidv=>n_samples,:nsnps=>nrow(cissnps),:window=>window,:typeofh2=>"cis",
                    :h2=>NaN,:se=>NaN,:Vg=>NaN,:Ve=>NaN,:genename=>genename])
                CSV.write(h2_file,h2_DF,append=true,delim="\t")
            else
                tmpdir = string(out_dir,"_temp/", tissue,".",gene,".Chr",chrom)
                mkpath(tmpdir)
                # Pull cis-SNP genotypes
                # extract snp from plink
                range_cissnps = string(chrom,"\t",start_pos,"\t",end_pos,"\t",gene)
                rangefile = string(tmpdir,"/",gene,".Chr",chrom,".cissnps.range.txt")
                write(rangefile,range_cissnps)
                # CSV.write(rangefile,cissnps,header=false,delim=" ")
                runplink = `plink --bfile $geno_file_prefix --extract range $rangefile --make-bed --out $tmpdir/$gene.Chr$chrom`
                run(runplink)
                # Pull expression data for gene
                exppheno = Matrix(expression[expression.gene_id .== gene,5:end])[1:end]
                outphe = DataFrame([:FID=>0,:IID=>sampleid,:EXP=>exppheno])
                phefile = string(tmpdir,"/",gene,".Chr",chrom,".phen")
                CSV.write(phefile,outphe,header=false,delim=" ")
                # grm
                rungcta = `gcta64 --bfile $tmpdir/$gene.Chr$chrom --make-grm --out $tmpdir/$gene.Chr$chrom`
                run(rungcta)
                try
                    # reml
                    rungcta = `gcta64 --grm $tmpdir/$gene.Chr$chrom --pheno $phefile --qcovar $outcov_file --reml-est-fix --reml --out $tmpdir/$gene.Chr$chrom`
                    run(rungcta)
                catch e
                    h2_DF = DataFrame([:geneid=>gene,:nidv=>n_samples,:nsnps=>nrow(cissnps),:window=>window,:typeofh2=>"cis",
                    :h2=>NaN,:se=>NaN,:Vg=>NaN,:Ve=>NaN,:genename=>genename])
                    CSV.write(h2_file,h2_DF,append=true,delim="\t")
                    rm(tmpdir, recursive=true, force=true)
                    continue
                end
                abc = CSV.read(string(tmpdir,"/",gene,".Chr",chrom,".hsq"),DataFrame,silencewarnings=true)
                cis_h2 = Dict(
                    :ΣG => abc[1,2],
                    :ΣG_se => abc[1,3],
                    :Σe => abc[2,2],
                    :Σe_se => abc[2,3],
                    :Σp => abc[3,2],
                    :Σp_se => abc[3,3],
                    :ΣG╱Σp => abc[4,2],
                    :ΣG╱Σp_se => abc[4,3],
                    :logL => abc[5,2],
                    :logL0 => abc[6,2],
                    :n => abc[10,2],
                    :Fix_eff => abc[12:end,1],
                    :Fix_eff_se => abc[12:end,2],
                    )
                #
                rm(tmpdir, recursive=true, force=true)
                #
                h2_DF = DataFrame([:geneid=>gene,:nidv=>n_samples,:nsnps=>nrow(cissnps),:window=>window,:typeofh2=>"cis",
                    :h2=>cis_h2[:ΣG╱Σp],:se=>cis_h2[:ΣG╱Σp_se],:Vg=>cis_h2[:ΣG],:Ve=>cis_h2[:Σe],:genename=>genename])
                CSV.write(h2_file,h2_DF,append=true,delim="\t")
                #
                coveff_out = string(gene,"\t",window,"\t")
                for x in cis_h2[:Fix_eff]
                    coveff_out = string(coveff_out,"\t",x)
                end
                coveff_out = string(coveff_out,"\n")
                coveff_file = string(out_dir,tissue,".chr",chrom,".cis_coveff.txt")
                open(coveff_file,"a+") do io
                    write(io,coveff_out)
                end
                #
                coveffse_out = string(gene,"\t",window,"\t")
                for x in cis_h2[:Fix_eff_se]
                    coveffse_out = string(coveffse_out,"\t",x)
                end
                coveffse_out = string(coveffse_out,"\n")
                coveffse_file = string(out_dir,tissue,".chr",chrom,".cis_coveffse.txt")
                open(coveffse_file,"a+") do io
                    write(io,coveffse_out)
                end
            end
        end
    end
end
