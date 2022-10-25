#!/usr/bin/env Rscript
rm(list=ls())
setwd("Pig/ICA/Script")

#-------Working----------------------
args<-commandArgs(T)
i=print(args[1])
length=nchar(i)
name<-substr(i, 1, length-23)
exp <- read.table(paste("Pig/Rawdata/New_tmm_removecov/",i,sep=""),header=T, sep="\t", check.names=F, stringsAsFactors = F)

#---------Run ICA-----------------------
icares = runICA(exp,n_runs = 15, n_cores = 10, scale_pheno = T, max_iter = 10, var_cutoff = 70)


#--------Save the components-----------
resfile = paste0("ICA_", name, ".RData") #combine
saveRDS(icares, file=file.path('Pig/ICA/1_RData', resfile))

ic_covar_mx = t(icares$A)
write.csv(ic_covar_mx, file.path('Pig/ICA/2_separate_eigen', paste("ICA", name, "separate_eigen.csv", sep="_")), quote = F)

#--------Save cluster info--------------
cl_ids = colnames(ic_covar_mx)
sds = apply(icares$S, 2, function(x)(abs(x)>2*stats::sd(x)))
sds = sds[,cl_ids]
cldf = melt(sds)
cldf = cldf[cldf$value,c("Var1", "Var2")]
colnames(cldf) = c("gene_id", "cl_id")
write.csv(cldf, file.path('Pig/ICA/3_separate_clusters', paste("ICA", name, "separate_clusters.csv", sep="_")), row.names = F)

#---------Save cluster number-----------
num=data.frame(table(cldf$cl_id))
names(num)=c("Cluster_ID","Number")
write.csv(num, file.path('Pig/ICA/4_separate_gene.num', paste("ICA", name, "separate_gene.num.csv", sep="_")), row.names = F)

print("ICA is done")
