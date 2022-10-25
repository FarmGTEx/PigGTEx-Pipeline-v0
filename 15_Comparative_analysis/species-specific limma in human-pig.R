library(limma)
library(xlsx)
library(readxl)
library(ggplot2)

load("/Users/baizhonghao/Downloads/human-pig GTEx/DataMatrix22761.RData")

humantissue_name<-unique(sort(Human_experiment$Tissue_new))
pigtissue_name<-unique(sort(Pig_experiment$Tissue_new))

human_pig_tissue<-intersect(pigtissue_name,humantissue_name)

#human-pig#

for (i in 1:length(human_pig_tissue)){
  tissue_meta<-experiment[which(experiment$Tissue_new==human_pig_tissue[i]),]
  tissue_exp<-expression[match(tissue_meta$SampleID,rownames(expression)),]
  tissue_exp<-log2(tissue_exp)
  tissue_exp<-t(tissue_exp)
  human_number<-length(which(tissue_meta$Species=='Human'))
  pig_number<-length(which(tissue_meta$Species=='Pig'))
  
  df<-c(rep("Human",human_number),rep("Pig",pig_number))
  table(df)
  design<-model.matrix(~-1+factor(df))
  dim(design)
  colnames(design)<-c("Pig","Human")
  rownames(design)<-colnames(tissue_exp)
  design
  contrastmatrix<-makeContrasts(Human-Pig,levels = design)
  contrastmatrix
  fit<-lmFit(tissue_exp,design)
  fit2<-contrasts.fit(fit,contrastmatrix)
  fit2<-eBayes(fit2)
  fit2
  myresults<-topTable(fit2,coef=1,adjust='BH',number=nrow(tissue_exp))
  write.table(myresults,file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/species-specific limma/origin_",human_pig_tissue[i],".txt"),quote = F)
}

Results<-NULL

human_pig_file<-dir("/Users/baizhonghao/Downloads/human-pig GTEx/species-specific limma",full.names = T)
human_pig_file
human_pig_tissue<-substr(human_pig_file,nchar("/Users/baizhonghao/Downloads/human-pig GTEx/species-specific limma/")+1,nchar(human_pig_file)-nchar(".txt"))
human_pig_tissue


for (i in 1:length(human_pig_tissue)){
  a<-read.table(file = paste0("/Users/baizhonghao/Downloads/human-pig GTEx/species-specific limma/",human_pig_tissue[i],".txt"))
  if(length(which(is.na(a$adj.P.Val)==TRUE))!=0){
    a<-a[-which(is.na(a$adj.P.Val)==TRUE),]
  }
  Summary<-a
  Summary$Species<-'Pig'
  Summary$Tissue<-human_pig_tissue[i]
  for(j in 1:length(rownames(Summary))){
    if(Summary$adj.P.Val[j]<0.01){
      if(Summary$logFC[j]>0){
        Summary$Change[j]<-'higher'
      }else if(Summary$logFC[j]<0){
        Summary$Change[j]<-'lower'
      }else{
        Summary$Change[j]<-"undifferential"
      }}
    else{
      Summary$Change[j]<-'undifferential'
    }
  }
  Results<-rbind(Results,Summary)
}

Results$Category<-paste0(Results$Species,sep='-',Results$Tissue)
category<-unique(Results$Category)
data_summary<-array(NA,dim = c(length(category),5))
colnames(data_summary)<-c('Species','Tissues','higher','undifferential','lower')
rownames(data_summary)<-category
data_summary<-as.data.frame(data_summary)
for (i in 1:length(category)){
  a<-Results[which(Results$Category==category[i]),]
  data_summary[category[i],'Species']<-unique(a$Species)
  data_summary[category[i],'Tissues']<-unique(a$Tissue)
  data_summary[category[i],'higher']<-length(which(a$Change=='higher'))
  data_summary[category[i],'undifferential']<-length(which(a$Change=='undifferential'))
  data_summary[category[i],'lower']<-length(which(a$Change=='lower'))
}
write.xlsx(data_summary,file=paste0('/Users/baizhonghao/Downloads/human-pig GTEx/species-specific limma/data_summary.xlsx'))

library(readxl)
data_summary<-as.data.frame(read_excel('/Users/baizhonghao/Downloads/human-pig GTEx/species-specific limma/data_summary.xlsx',sheet = 'Sheet1',col_names = T))
data_summary$Changes<-factor(data_summary$Changes,levels = c('Higher expression in human','not differentially expressed','Lower expression in human'))
#tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/species_specific.tiff",
#     res = 300, width = 2000, height = 1500,compression = "lzw")
par(mar=c(5.1,6,6,10))
p1<-ggplot(data_summary,mapping = aes(x=Tissues,y=number,fill=Changes))+
  geom_bar(stat='identity',position="fill")+
  guides(fill = guide_legend(reverse = F))+
  labs(x = 'Tissues',y = 'Frequnency') +
  theme(axis.title =element_text(size = 16),axis.text =element_text(size = 6, color = 'black'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=5.2))+
  theme(axis.title.y.left = element_text(vjust = 2))+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(legend.title = element_text(face="bold"))+
  theme(axis.title=element_text(face="bold"))+
  scale_fill_manual(values = c('Higher expression in human'="#BD3A1E",'not differentially expressed'="#999696",'Lower expression in human'="#2C418E"))
#dev.off()
ggsave(p1,file='/Users/baizhonghao/Downloads/species-specific.pdf',dpi=300,width=8,height = 6)
