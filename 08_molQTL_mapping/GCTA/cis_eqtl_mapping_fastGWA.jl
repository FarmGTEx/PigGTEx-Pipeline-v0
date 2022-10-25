# using GLMNet
using Base.Threads
using Statistics
using HypothesisTests
using Random
using DataFrames
using StatsBase
using CSV
# using SnpArrays
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
home_dir = "/BIGDATA1/scau_hzhang_1/USER/tengjy/Pig_GTEx_eQTL/eqtl_mapping_mlm"
cd(home_dir)
# tissue=["Kidney"][1]  # name of tissue/study
tissue=TIS
expression_bed_file = string("../eqtl_mapping/bed/",tissue,".expr_tmm_inv.bed.gz")
geno_file_prefix = string("../eqtl_mapping/geno/",tissue,".filtered_maf0.05_mac6")
cov_file = string("../eqtl_mapping/covar_5or10pc10peer/",tissue,".covariates.txt")
snp_annot_file = string("../genotypes/filtered/snp_annot_maf0.05_mac6/",tissue,".snp.map")
# gene_annot_file = "../eqtl_mapping/Sus_scrofa.Sscrofa11.1.100.tss.gz"
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
snp_annot_all = CSV.File(snp_annot_file,header=false) |> DataFrame
covar = CSV.File(cov_file,transpose=true) |> DataFrame
rename!(covar,:Column1=>:SampleID)
#
mkpath(out_dir)
# Pull covariates
outcov = hcat(repeat([0],n_samples),covar)
outcov = outcov[:,["x1","SampleID","peer1","peer2","peer3","peer4","peer5","peer6","peer7","peer8","peer9","peer10"]]
tmpdir=string(out_dir,"_temp/")
mkpath(tmpdir)
# for chrom in 1:18
for chrom in parse.(Int,CHR)
    # chrom = 1
    path_Chr = string(out_dir,"Chr",chrom,"/")
    mkpath(path_Chr)
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
            cissnps = snp_annot[(snp_annot.Column2 .> start_pos) .& (snp_annot.Column2 .<= end_pos),:]
            if nrow(cissnps) < 1
                # Need 1 or more cis-snps
                @warn string("Zero cis-snps for ", gene)
            else
                tmpdir = string(out_dir,"_temp/", tissue,".",gene,".Chr",chrom)
                mkpath(tmpdir)
                # Pull cis-SNP genotypes
                # extract snp from plink
                range_cissnps = string(chrom,"\t",start_pos,"\t",end_pos,"\t",gene)
                rangefile = string(tmpdir,"/",gene,".Chr",chrom,".cissnps.range.txt")
                write(rangefile,range_cissnps)
                runplink = `plink --bfile $geno_file_prefix --extract range $rangefile --make-bed --out $tmpdir/$gene.Chr$chrom`
                run(runplink)
                # Pull expression data for gene
                exppheno = Matrix(expression[expression.gene_id .== gene,5:end])[1:end]
                outphe = DataFrame([:FID=>0,:IID=>sampleid,:EXP=>exppheno])
                phefile = string(tmpdir,"/",gene,".Chr",chrom,".phen")
                CSV.write(phefile,outphe,header=false,delim=" ")
                #
                if nrow(cissnps) < 100
                    # mlma
                    @warn string("Less than 100 cis-snps for ", gene,", '--fastGWA-mlm' was replaced with '--mlma'.")
                    rungcta = `gcta64 --bfile $tmpdir/$gene.Chr$chrom --grm ./GRM/$tissue --mlma --pheno $phefile --qcovar $outcov_file --thread-num 23 --out $path_Chr$tissue.Chr$chrom.$gene.assoc`
                else
                    # fastGWA-mlm
                    rungcta = `gcta64 --bfile $tmpdir/$gene.Chr$chrom --grm-sparse ./GRM/$tissue --fastGWA-mlm --pheno $phefile --qcovar $outcov_file --threads 23 --nofilter --out $path_Chr$tissue.Chr$chrom.$gene.assoc`
                end                
                run(rungcta)
                
                ### Check whether the task have done.
                assoc_fastGWA_file = string(path_Chr,tissue,".Chr",chrom,".",gene,".assoc.fastGWA")
                assoc_mlma_file = string(path_Chr,tissue,".Chr",chrom,".",gene,".assoc.mlma")
                if !isfile(assoc_fastGWA_file) & !isfile(assoc_mlma_file)
                    @warn string("Perform '--mlma' for ", gene,".")
                    # mlma
                    rungcta = `gcta64 --bfile $tmpdir/$gene.Chr$chrom --grm ./GRM/$tissue --mlma --pheno $phefile --qcovar $outcov_file --thread-num 23 --out $path_Chr$tissue.Chr$chrom.$gene.assoc`
                    try
                        run(rungcta)
                        rm(tmpdir, recursive=true, force=true)
                    catch e
                        @warn string("Error for ", gene,".")
                        continue
                    end
                end                
                rm(tmpdir, recursive=true, force=true)                
            end
        end
    end
end
