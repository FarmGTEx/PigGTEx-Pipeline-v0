using Base.Threads
using Statistics
# using HypothesisTests
using Random
using DataFrames
# using StatsBase
using CSV
using Dates
using Mmap
using GZip
#==============================================================================#
#==============================================================================#
CHR=ARGS[1]
TIS=ARGS[2]
# Parameters setting ----------------------------------------------------------#
home_dir = "/BIGDATA1/scau_hzhang_1/USER/tengjy/Pig_GTEx_eQTL/eqtl_mapping_mlm"
cd(home_dir)
# tissue=["Kidney"][1]  # name of tissue/study
tissue=TIS
expression_bed_file = string("../eqtl_mapping/bed/",tissue,".expr_tmm_inv.bed.gz")
out_dir=string("./output_fastGWA_10peer/", tissue, "/")
# window=1e6
# seed=2021
# -----------------------------------------------------------------------------#
# Load data from disk
io = GZip.open(expression_bed_file)
expression = CSV.File(io,header=true) |> DataFrame
close(io)
n_samples = ncol(expression) - 4
gene_annot_all = expression[:,1:4]
#
mkpath(out_dir)
# for chrom in 1:18
for chrom in parse.(Int,CHR)
    # chrom = 1
    path_eqtls = string(out_dir,"Chr",chrom,"/")
    gene_annot = gene_annot_all[gene_annot_all."#Chr" .== chrom,:]
    exp_genes = gene_annot.gene_id
    n_genes = nrow(gene_annot)
    window = Int.(1e6)
    chr_assoc = DataFrame(:phenotype_id=>Vector{String}[],:variant_id=>Vector{String}[],:tss_distance=>Vector{Int64}[],:af=>Vector{Float64}[],:pval_mlm=>Vector{Float64}[],:slope=>Vector{Float64}[],:slope_se=>Vector{Float64}[])
    top_assoc = DataFrame(:phenotype_id=>Vector{String}[],:variant_id=>Vector{String}[],:tss_distance=>Vector{Int64}[],:af=>Vector{Float64}[],:pval_mlm=>Vector{Float64}[],:slope=>Vector{Float64}[],:slope_se=>Vector{Float64}[])
    for i in 1:n_genes
        # i = 1
        println("GENE: ",i,"/",n_genes)
        gene = exp_genes[i]
        geneinfo = gene_annot[gene_annot.gene_id .== gene,:]
        tss_pos = geneinfo.end[1]
        assoc_fastGWA_file = string(path_eqtls,tissue,".Chr",chrom,".",gene,".assoc.fastGWA")
        assoc_mlma_file = string(path_eqtls,tissue,".Chr",chrom,".",gene,".assoc.mlma")
        
        if !isfile(assoc_fastGWA_file) & !isfile(assoc_mlma_file)
            @warn string("Cannot found assoc file: ", gene, ".")
            continue
        elseif isfile(assoc_fastGWA_file)
            assoc_fastGWA = CSV.File(assoc_fastGWA_file) |> DataFrame
            tss_distance = assoc_fastGWA.POS .- tss_pos
            gene_assoc = DataFrame(:phenotype_id=>repeat([gene],nrow(assoc_fastGWA)),:variant_id=>assoc_fastGWA.SNP,:tss_distance=>tss_distance,:af=>assoc_fastGWA.AF1,:pval_mlm=>assoc_fastGWA.P,:slope=>assoc_fastGWA.BETA,:slope_se=>assoc_fastGWA.SE)
        elseif isfile(assoc_mlma_file)
            assoc_mlma = CSV.File(assoc_mlma_file) |> DataFrame
            tss_distance = assoc_mlma.bp .- tss_pos
            gene_assoc = DataFrame(:phenotype_id=>repeat([gene],nrow(assoc_mlma)),:variant_id=>assoc_mlma.SNP,:tss_distance=>tss_distance,:af=>assoc_mlma.Freq,:pval_mlm=>assoc_mlma.p,:slope=>assoc_mlma.b,:slope_se=>assoc_mlma.se)
        end
        minidx = findmin(gene_assoc.pval_mlm)[2]
        if nrow(chr_assoc) == 0
            chr_assoc = gene_assoc
            top_assoc = gene_assoc[minidx:minidx,:]
        else
            chr_assoc = vcat(chr_assoc,gene_assoc, cols=:orderequal)
            top_assoc = vcat(top_assoc,gene_assoc[minidx:minidx,:], cols=:orderequal)
        end
    end
    CSV.write(string(out_dir,tissue,".cis_qtl_pairs.",chrom,".assoc.fastGWA"),chr_assoc,delim="\t")
    run(`gzip $out_dir$tissue.cis_qtl_pairs.$chrom.assoc.fastGWA`)
    CSV.write(string(out_dir,tissue,".cis_qtl_pairs.",chrom,".best.fastGWA"),top_assoc,delim="\t")
end




