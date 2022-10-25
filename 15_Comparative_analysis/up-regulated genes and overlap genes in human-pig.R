
rm(list = ls())


##########3
commontissue<-rownames(Results1)
Results<-array(data = NA,dim = c(length(commontissue),3))
tissue_name<-commontissue
rownames(Results)<-tissue_name
colnames(Results)<-c("human_up","pig_up","overlap")


logFC_cutoff<-1.5
for(i in 1:length(commontissue)){
  tissue_human<-read.table(file = paste0("/Users/baizhonghao/Downloads/human-pig GTEx/tissue-specific limma/human/",tissue_name[i],".txt"))
  tissue_pig<-read.table(file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/tissue-specific limma/pig/",tissue_name[i],".txt"))

  tissue_human$change = as.factor(ifelse(tissue_human$adj.P.Val < 0.05 & abs(tissue_human$logFC) > logFC_cutoff,
                                         ifelse(tissue_human$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  
  human_upgene<-tissue_human[tissue_human$change =='UP',]
  human_upgene<-na.omit(human_upgene)
  write.table(human_upgene,file=paste0("~/Downloads/human-pig GTEx/tissue-specific limma/human-up/",tissue_name[i],".txt"))
  Results[i,1]<- length(rownames(human_upgene)) 
  
  tissue_pig$change = as.factor(ifelse(tissue_pig$adj.P.Val < 0.05 & abs(tissue_pig$logFC) > logFC_cutoff,
                                       ifelse(tissue_pig$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  
  pig_upgene<-tissue_pig[tissue_pig$change =='UP',]
  pig_upgene<-na.omit(pig_upgene)
  write.table(pig_upgene,file=paste0("~/Downloads/human-pig GTEx/tissue-specific limma/pig-up/",tissue_name[i],".txt"))
  Results[i,2]<- length(rownames(pig_upgene))    
  
  overlap<-intersect(rownames(human_upgene),rownames(pig_upgene))
  Results[i,3]<-length(overlap)
  write.table(overlap,file=paste0("~/Downloads/human-pig GTEx/tissue-specific limma/overlap-up/",tissue_name[i],".txt"))                                                            
}



