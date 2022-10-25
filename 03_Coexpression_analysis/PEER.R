#!/usr/bin/env Rscript
rm(list=ls())
setwd("Pig/PEER/Script")
suppressPackageStartupMessages(library("argparse"))

#----------Working-----------
library(dplyr)
library(PLIER)
library(peer)
library(reshape2)

args<-commandArgs(T)
i=print(args[1])
length=nchar(i)
name<-substr(i, 1, length-23)

exp <- read.table(paste("Pig/Rawdata/New_tmm_removecov/",i,sep=""),header=T, sep="\t", check.names=F, stringsAsFactors = F)

std = apply(exp, 1, sd)
exp = exp[std!=0,]

exp = as.matrix(exp)
mat_scale = rowNorm(exp)
nr.pc = num.pc(mat_scale)

mat_scale = as.matrix(t(mat_scale))

# Create the model object
print("Creating the model object")
model = PEER()

# set the observed data
PEER_setPhenoMean(model, mat_scale)
print(dim(PEER_getPhenoMean(model)))

# we want to infer K hidden confounders
PEER_setNk(model, nr.pc)

# include covariate to account for the mean expression
print("Account for mean expression")
PEER_setAdd_mean(model, TRUE)

# perform the inference
print("Performing the inference")
PEER_update(model)

# observing output
print("Get factors")
factors = PEER_getX(model)

print("Get weights")
weights = PEER_getW(model)

print("Get precision")
precision = PEER_getAlpha(model)

print("Get residuals")
residuals = PEER_getResiduals(model)

print("Merge the results")
res = list()
res[["factors"]] = factors
res[["weights"]] = weights
res[["precision"]] = precision
res[["residuals"]] = residuals

print("PEER is done")

#-------Export--------
group = name
resfile = paste0("PEER_", group, ".rds")
saveRDS(res, file=file.path('Pig/PEER/1_RDS', resfile))

res_mat = data.frame(res$factors)
res_mat = res_mat[,-1]

colnames(res_mat) = paste0("F", seq(1, ncol(res_mat)))

sample_order = colnames(exp)
row.names(res_mat) = sample_order

write.csv(res_mat, file.path('Pig/PEER/1_separate_eigen', paste("PEER", name, "separate_eigen.csv", sep="_")), quote = F)

weights <- res$weights[,-1]
genenames = row.names(exp)
row.names(weights) = genenames
colnames(weights) = paste0("F", seq(1, dim(weights)[2]))

sd.outlier <- function(x) abs(x-mean(x))/sd(x) > 2
sds = apply(weights, 2, sd.outlier)
cldf = melt(sds)
cldf = cldf[cldf$value,c("Var1", "Var2")]
colnames(cldf) = c("gene_id", "cl_id")
write.csv(cldf, file.path('Pig/PEER/2_separate_clusters', paste("PEER", name, "separate_clusters.csv", sep="_")), row.names = F)

num=data.frame(table(cldf$cl_id))
names(num)=c("Cluster_ID","Number")
write.csv(num, file.path('Pig/PEER/3_separate_gene.num' , paste("PEER", name, "separate_gene.num.csv", sep="_")), row.names = F)

