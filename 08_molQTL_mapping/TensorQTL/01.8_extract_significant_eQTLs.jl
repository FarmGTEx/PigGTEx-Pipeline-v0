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

file_cis_qtl_fdr = ARGS[1]  # "./output_pcg/",tis,"/",tis,".cis_qtl_fdr0.05.txt.gz"
nominal_dir = ARGS[2] # "./output_pcg/",tis,"/"
output_dir = ARGS[3]
study_name = ARGS[4]

function gogogo(file_cis_qtl_fdr, nominal_dir, output_dir, study_name)
    mkpath(output_dir)
    io = GZip.open(file_cis_qtl_fdr)
    pval_threshold = CSV.File(io) |> DataFrame    
    close(io)
    dataSig = DataFrame(:phenotype_id=>[],:variant_id=>[],:tss_distance=>[],:af=>[],:ma_samples=>[],:ma_count=>[],:pval_nominal=>[],:slope=>[],:slope_se=>[])
    for CHR in 1:18
        println(CHR)
        io = GZip.open(string(nominal_dir,"/",study_name,".cis_qtl_pairs.",CHR,".txt.gz"))
        eqtls_lr = CSV.File(io) |> DataFrame
        close(io)
        
        eqtls_lr = eqtls_lr[.~isnothing.(indexin(eqtls_lr.phenotype_id,pval_threshold.phenotype_id)),:]
        idx = .!ismissing.(eqtls_lr.pval_nominal)
        eqtls_lr = eqtls_lr[idx,:]
        
        genes = unique(eqtls_lr.phenotype_id)
        n_genes = length(genes)
        dataSig_tmp = DataFrame(:phenotype_id=>[],:variant_id=>[],:tss_distance=>[],:af=>[],:ma_samples=>[],:ma_count=>[],:pval_nominal=>[],:slope=>[],:slope_se=>[])
        for i in 1:n_genes
            idx_gene = eqtls_lr.phenotype_id .== genes[i]
            eqtls_gene = eqtls_lr[idx_gene,:]
            infogene = pval_threshold[pval_threshold.phenotype_id .== genes[i],:]
            if infogene.is_eGene[1] == true
                p_cutoff = infogene.pval_nominal_threshold[1]
                idx_eSNP = eqtls_gene.pval_nominal .< p_cutoff
                if sum(idx_eSNP) > 0
                    if nrow(dataSig_tmp) == 0
                        dataSig_tmp = eqtls_gene[idx_eSNP,:]
                    else
                        dataSig_tmp = vcat(dataSig_tmp,eqtls_gene[idx_eSNP,:])
                    end
                end
            end
        end
        
        if nrow(dataSig) == 0
            dataSig = dataSig_tmp
        else
            dataSig = vcat(dataSig,dataSig_tmp)
        end
    end
    CSV.write(string(output_dir,"/",study_name,".cis_qtl_pairs.significant.txt"),dataSig,delim="\t")
    run(`gzip $output_dir/$study_name.cis_qtl_pairs.significant.txt`)
end
gogogo(file_cis_qtl_fdr, nominal_dir, output_dir, study_name)


