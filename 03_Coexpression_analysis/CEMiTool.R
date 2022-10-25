#!/usr/bin/env Rscript
rm(list = ls())
setwd("Pig/CEMitools/Script")

#---------Raw data----------------
args<-commandArgs(T)
i=print(args[1])
length=nchar(i)
name<-substr(i, 1, length-23)
expr0 <- read.table(paste("Pig/Rawdata/New_tmm_removecov/",i,sep=""),
                    header=T, sep="\t", check.names=F, stringsAsFactors = F)

#----style, echo=FALSE, results="asis", message=FALSE---------
knitr::opts_chunk$set(tidy = FALSE,
                      warning = FALSE,
                      message = FALSE,
                      cache=TRUE)

#------------------------------------------------------------
library("CEMiTool")
cem <- cemitool(expr0)

ME <- mod_summary(cem)
ME <- filter(ME, ME$modules!='Not.Correlated')
write.csv(ME, file.path('Pig/CEMitools/1_separate_eigen', paste("CEMiTool", name, "separate_eigen.csv", sep="_")), quote = F, row.names = F)

module <- module_genes(cem)
module <- filter(module,module$modules!='Not.Correlated' )
write.csv(module, file.path('Pig/CEMitools/2_separate_clusters', paste("CEMiTool", name, "separate_clusters.csv", sep="_")), row.names = F)

module <- as.data.frame(module)
Sta = data.frame(table(module$modules))
names(Sta)[1] <- c("module")
Sta <- filter(Sta,Sta$module!='Not.Correlated')
Sta$Tissue=name
Sta$Module.num=length(Sta$module)
write.csv(Sta, file.path('Pig/CEMitools/3_separate_gene.num', paste("CEMiTool", name, "separate_gene.Num.csv", sep="_")), quote = F, row.names = F)

print("CEMitools is done")
