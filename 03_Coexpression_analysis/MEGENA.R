#!/usr/bin/env Rscript
#----Parameters-------
rm(list=ls())
setwd("Pig/MEGENA/Script/Blastula")
library(MEGENA)
# input parameters
n.cores <- 24; # number of cores/threads to call for PCP
doPar <- TRUE; # do we want to parallelize
method = "pearson" # method for correlation.
FDR.cutoff = 0.05 # FDR threshold to define significant correlations upon shuffling samples.
module.pval = 0.05 # module significance p-value. Recommended is 0.05.
hub.pval = 0.05 # connectivity significance p-value based random tetrahedral networks
cor.perm = 10; # number of permutations for calculating FDRs for all correlation pairs.
hub.perm = 100; # number of permutations for calculating connectivity significance p-value.

# annotation to be done on the downstream
annot.table=NULL #sample_annot
id.col = 1
symbol.col= 2

#----Working----------
args<-commandArgs(T)
i=print(args[1])
length=nchar(i)
name<-substr(i, 1, length-23)
datExpr <- read.table(paste("Pig/Rawdata/New_tmm_removecov/",i,sep=""),header=T, sep="\t", check.names=F, stringsAsFactors = F)

#-----calculate correlation-------------------#
ijw <- calculate.correlation(datExpr,doPerm = cor.perm)


#------calculate PFN------
run.par = doPar & (getDoParWorkers() == 1)
if(run.par)
{
  cl <- parallel::makeCluster(n.cores)
  registerDoParallel(cl)
  cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep =""))
}


el <- calculate.PFN(ijw[,1:3],doPar = doPar,num.cores = n.cores)

g <- graph.data.frame(el,directed = FALSE)

#-----perform clustering---------
##### perform MCA clustering.
MEGENA.output <- do.MEGENA(g,
                           mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = TRUE,
                           min.size = 150,max.size = vcount(g)/2,
                           doPar = doPar,num.cores = n.cores,n.perm = hub.perm,
                           save.output = TRUE)

resfile = paste0("MEGENA_", name, ".RData")
saveRDS(MEGENA.output, file=file.path('Pig/MEGENA/1_RData', resfile))

#------summarize results------------
summary.output <- MEGENA.ModuleSummary(MEGENA.output,
                                       mod.pvalue = module.pval,hub.pvalue = hub.pval,
                                       min.size = 150,max.size = vcount(g)/2,
                                       annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                       output.sig = TRUE)

if (!is.null(annot.table))
{
  # update annotation to map to gene symbols
  V(g)$name <- paste(annot.table[[symbol.col]][match(V(g)$name,annot.table[[id.col]])],V(g)$name,sep = "|")
  summary.output <- output[c("mapped.modules","module.table")]
  names(summary.output)[1] <- "modules"
}

print(head(summary.output$modules,2))

print(summary.output$module.table)

summary.output$module.table

#------Get module table---------
module.df = module_convert_to_table(MEGENA.output,mod.pval = module.pval,
                                    hub.pval = hub.pval,min.size = 150,max.size = vcount(g)/2)

write.csv(module.df, file.path('Pig/MEGENA/3_separate_clusters', paste("MEGENA", name, "separate_clusters.csv", sep="_")), row.names = F)

#------Save cluster number--------
num=data.frame(table(module.df$module))
names(num)=c("Cluster_ID","Number")
write.csv(num, file.path('Pig/MEGENA/4_separate_gene.num' , paste("MEGENA", name, "separate_gene.num.csv", sep="_")), row.names = F)

#------Calculate PC---------------
datExpr <- data.frame(datExpr)
library(dplyr)
exp <-tibble::rownames_to_column(datExpr, "id")
dim(exp)

#----Loop Calculate PC---------------
module_name<-as.vector(unique(module.df$module))
modulenumber <- length(module_name)

for(i in 1:modulenumber){
  if(i ==1 ){
    c1_2<- subset(module.df,module==module_name[1])
    c1_2<- data.frame(id=c1_2[,1])

    subExp=left_join(c1_2,exp,by="id")
    rownames(subExp)=subExp$id
    subExp= data.frame(t(subExp[,-1]))

    df_pca <- prcomp(subExp)
    df_pc1<-data.frame(df_pca$x[,1])
    file_all<- df_pc1
  }else{
    c1_2<- subset(module.df,module==module_name[i])
    c1_2<- data.frame(id=c1_2[,1])

    subExp=left_join(c1_2,exp,by="id")
    rownames(subExp)=subExp$id
    subExp= data.frame(t(subExp[,-1]))

    df_pca <- prcomp(subExp)
    df_pc1<-data.frame(df_pca$x[,1])
    file_all<-cbind(file_all,df_pc1)
  }
}
names(file_all)[1:ncol(file_all)] <- as.character(module_name)
write.csv(file_all, file.path('Pig/MEGENA/2_separate_eigen', paste("MEGENA", name, "separate_eigen.csv", sep="_")), quote = F)

print("MEGENA is done")
