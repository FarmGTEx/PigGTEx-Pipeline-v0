#!/usr/bin/env Rscript
rm(list=ls())
setwd("Pig/WGCNA/Script")

#----Parameters -------
library("SummarizedExperiment")
library(dplyr)
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
library(data.table)
library(R.utils)

#----Working----------
args<-commandArgs(T)
i=print(args[1])
length=nchar(i)
name<-substr(i, 1, length-23)
exp <- read.table(paste("Pig/Rawdata/New_tmm_removecov/",i,sep=""),header=T, sep="\t", check.names=F, stringsAsFactors = F)


std = apply(exp, 1, sd)
exp = exp[std!=0,]

net = blockwiseModules(t(exp),
                       TOMType ="unsigned",
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = F, corType = "pearson",
                       saveTOMFileBase = "ExprTOM",
                       verbose = 3)

resfile = paste0("WGCNA_", name, ".RData")
saveRDS(net, file=file.path('Pig/WGCNA/1_RData', resfile))

moduleColors = labels2colors(net$colors)
Sta = data.frame(table(moduleColors))
Sta$Tissue=group
Sta$Module.num=length(Sta$moduleColors)
write.csv(Sta, file.path('Pig/WGCNA/2_separate_gene.num', paste("WGCNA", name, "separate_gene.Num.csv", sep="_")), quote = F, row.names = F)


res_mat = data.frame(net$MEs)
sample_order = colnames(exp)
row.names(res_mat) = sample_order

res_mat = res_mat[,colnames(res_mat)!="MEgrey"]
write.csv(res_mat, file.path('Pig/WGCNA/3_separate_eigen', paste("WGCNA", name, "separate_eigen.csv", sep="_")), quote = F)

gene.name = data.frame(gene_id=row.names(exp))
moduleLabels = data.frame(net.colors=paste0("ME",net$colors))
moduleColors = data.frame(labels2colors(net$colors))

geneInfo0<-cbind.data.frame("gene"=rownames(moduleLabels),gene.name,moduleLabels,moduleColors)
names(geneInfo0)<-c("gene","gene_id","module_MEnumber","module_Color")
geneInfo0 <- geneInfo0[which(geneInfo0$module_Color!="grey"),]
write.table(geneInfo0, file.path('Pig/WGCNA/4_separate_clusters', paste("WGCNA", name, "separate_geneInfo.tsv", sep="_")), quote = F, sep = "\t",row.names = F)

print("WGCNA is done")


