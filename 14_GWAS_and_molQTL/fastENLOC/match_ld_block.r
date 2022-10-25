library(data.table)
library("dplyr")
library("parallel")

argv<-commandArgs(T)
setwd(argv[1])
gwas<-data.table::fread(paste0(argv[2],".txt.gz"),header=T)
gwas=as.data.frame(gwas)
gwas = gwas %>%
select(variant_id,chromosome,position,zscore)
gwas<-as.data.frame(gwas)


block<-data.table::fread(argv[3],header=T)
block<-cbind(block,Loc=paste0("Loc",c(1:dim(block)[1])))
block<-as.data.frame(block)
block[,1]<-as.numeric(block[,1])
block[,2]<-as.numeric(block[,2])
block[,3]<-as.numeric(block[,3])
gwas[,3]<-as.numeric(gwas[,3])

loc<-function(i){
sub_gwas<-gwas[which(gwas[,2]==paste0("chr",i)),]
sub_block<-block[which(block[,1]==i),]
block_seq<-seq(1,dim(sub_block)[1],length.out=501)
block_seq<-round(block_seq)
null<-matrix(0,1,3)
for(j in 1:500){
data_block=sub_block[c(block_seq[block_seq]:(block_seq[block_seq+1]-1)),]
data_gwas=sub_gwas[which(sub_gwas[,3]>=data_block[1,2]&sub_gwas[,3]<=data_block[dim(data_block)[1],3]),]
list_gwas<-matrix(0,dim(data_gwas)[1],3)
b<-matrix(0,1,3)
for(z in 1:dim(data_block)[1]){
x2<-which(data_gwas[,3]>=data_block[z,2]&data_gwas[,3]<=data_block[z,3])
if(length(x2)>0){
b<-cbind(data_gwas[x2,1],data_block[z,4],data_gwas[x2,4])
b<-as.matrix(b)
}
list_gwas[x2,]<-b
}
x3<-which(list_gwas[,1]==0)
x3<-as.numeric(x3)
list_gwas[x3,]<-cbind(data_gwas[x3,1],data_gwas[x3,1],data_gwas[x3,4])
c<-rbind(c,list_gwas,use.names=F)
}
c<-c[-1,]
inter_block<-match(sub_gwas[,1],c[,1])
inter_block<-which(is.na(inter_block))
inter_block<-as.numeric(inter_block)
inter_block<-cbind(sub_gwas[inter_block,1],sub_gwas[inter_block,1],sub_gwas[inter_block,4])
c<-rbind(c,inter_block,use.names=F)
write.table(c,paste0(argv[4],argv[2],".chr",i,".txt"),col.names=F,row.names=F,quote=F)

}

detectCores()
cl = makeCluster(18)
clusterExport(cl, varlist = c("gwas","block","argv"))
clusterEvalQ(cl, "data.table")
parLapply(cl, 1:18,  loc)
stopCluster(cl)


