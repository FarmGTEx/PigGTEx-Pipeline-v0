library(data.table)
library(ggplot2)
library(readxl)
library(clusterProfiler)
library(dplyr)
#####################################################
#Load data#
#Human#
human_exp<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.txt"))
human_meta<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"))
tissuetype<-unique(human_meta$SMTS)

#Pig#
load("/Users/baizhonghao/Downloads/human-pig GTEx/ALL_Sample_Big_Matrix_20210104_tpm_7095.Rdata")
pig_meta <- read_excel("/Users/baizhonghao/Downloads/human-pig GTEx/Metadata_eQTLmapping_20210702_tjy.xlsx",sheet = "Metatable3_rmDup", col_names = T)
pig_exp<- as.data.frame(t(ALL_Sample))
pig_file<-dir("/Users/baizhonghao/Downloads/human-pig GTEx/SampleID_each_tissue",full.names = T)
pig_file
pig_tissue<-substr(pig_file,nchar("/Users/baizhonghao/Downloads/human-pig GTEx/SampleID_each_tissue/SampleID_")+1,nchar(pig_file)-nchar(".txt"))
pig_tissue
#Homologous annotation#
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
pig_anno<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/pig_annotation.txt"))

############################################
#Preparing the tissue data#
#Human#
Results <- array(NA,dim = c(length(human_exp$Name),length(tissuetype)))
colnames(Results)<-tissuetype
rownames(Results)<-human_exp$Name

for(i in 1:length(tissuetype)){
  
  a<-human_meta[human_meta$SMTS==tissuetype[i],]
  sampleID<-a$SAMPID
  sampleID<-intersect(colnames(human_exp),sampleID)
  Result<-select(human_exp, c(Name, sampleID))
  rownames(Result)<-Result$Name
  Result<-Result[,-1]
  Result$median_value<-apply(Result,1,median)
  Result<-select(Result,c(median_value))
  Results[,i]<-Result$median_value
}
rownames(Results)<-substr(rownames(Results),1,15)
annotation$length<-annotation$`Gene end (bp)` - annotation$`Gene start (bp)` + 1

Human_ortho<-NULL
for(i in 1:length(tissuetype)){
  summary <- array(NA,dim = c(length(human_exp$Name),2))
  colnames(summary)<-c('Gene','median')
  summary<-as.data.frame(summary)
  summary$Gene<-rownames(Results)
  summary$median<-Results[,tissuetype[i]]
  summary$Tissue<-tissuetype[i]
  summary$Species<-'Human'
  summary$length<-annotation$length[match(summary$Gene,annotation$`Gene stable ID`)]
  summary$orthology<-annotation$`Pig homology type`[match(summary$Gene,annotation$`Gene stable ID`)]
  summary$orthology[which(summary$orthology=='ortholog_one2many'|summary$orthology=='ortholog_many2many')]<-'Complex'
  summary$orthology[which(summary$orthology!='ortholog_one2one'&summary$orthology!='Complex')]<-'Others'
  summary$reads<-summary$median * summary$length
  summary<-na.omit(summary)
  Human_ortho<-rbind(Human_ortho,summary)
}

#Pig#
Results_pig <- array(NA,dim = c(length(rownames(pig_exp)),length(pig_tissue)))
colnames(Results_pig)<-pig_tissue
rownames(Results_pig)<-rownames(pig_exp)
pig_meta$`Categories-class2`[which(pig_meta$`Categories-class2`=="Large intestine")]<-"Large_intestine"
for(i in 1:length(pig_tissue)){
  a<-pig_meta[pig_meta$`Categories-class2`==pig_tissue[i]|pig_meta$`Categories-class3`==pig_tissue[i],]
  sampleID<-a$BioSample
  sampleID<-intersect(colnames(pig_exp),sampleID)
  Result<-select(pig_exp, c(sampleID))
  Result$median_value<-apply(Result,1,median)
  Result<-select(Result,c('median_value'))
  Results_pig[,i]<-Result$median_value
} 
pig_anno$length<- pig_anno$`Gene end (bp)` - pig_anno$`Gene start (bp)` + 1

Pig_ortho<-NULL
for(i in 1:length(pig_tissue)){
  summary <- array(NA,dim = c(length(rownames(Results_pig)),2))
  colnames(summary)<-c('Gene','median')
  summary<-as.data.frame(summary)
  summary$Gene<-rownames(Results_pig)
  summary$median<-Results_pig[,pig_tissue[i]]
  summary$Tissue<-pig_tissue[i]
  summary$Species<-'Pig'
  summary$length<-pig_anno$length[match(summary$Gene,pig_anno$`Gene stable ID`)]
  summary$orthology<-annotation$`Pig homology type`[match(summary$Gene,annotation$`Pig gene stable ID`)]
  summary$orthology[which(summary$orthology=='ortholog_one2many'|summary$orthology=='ortholog_many2many')]<-'Complex'
  summary$orthology[which(is.na(summary$orthology)==TRUE)]<-"Others"
  summary$orthology[which(summary$orthology!='ortholog_one2one'&summary$orthology!='Complex')]<-'Others'
  summary$reads<-summary$median * summary$length
  summary<-na.omit(summary)
  Pig_ortho<-rbind(Pig_ortho,summary)
}

#Integrate two species#
ortho<-rbind(Human_ortho,Pig_ortho)
ortho$orthology[which(ortho$orthology=='ortholog_one2one')]<-"1-1"
###############################################################################################

#Plot the figures#

##figures according to the reads##

p1<-ggplot(ortho,mapping = aes(x=Tissue,y=reads,fill=orthology))+
  geom_bar(stat='identity',position="fill") +
  facet_wrap(~Species, scales = 'free_x', ncol = 3) +
  guides(fill = guide_legend(reverse = F))+
  labs(x = 'Tissues',y = 'Proportion of reads') +
  theme(axis.title =element_text(size = 10),axis.text =element_text(size = 6, color = 'black'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 6))+
  scale_fill_manual(values = c("1-1"="#DC9F0A","Complex"="#37937F","Others"="#2C418E"))
ggsave(p1,file = "/Users/baizhonghao/Downloads/orthologous(1-1) contribute to transcriptome.pdf",height=2.5,width =7,dpi=300)

tiff(file = "/Users/baizhonghao/Downloads/orthologous(1-1) contribute to transcriptome.tiff",
     res = 300, width = 2500, height = 1000,compression = "lzw")
ggplot(ortho,mapping = aes(x=Tissue,y=reads,fill=orthology))+
  geom_bar(stat='identity',position="fill") +
  facet_wrap(~Species, scales = 'free_x', ncol = 3) +
  guides(fill = guide_legend(reverse = F))+
  labs(x = 'Tissues',y = 'Proportion of reads') +
  theme(axis.title =element_text(size = 12),axis.text =element_text(size = 6, color = 'black'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c("1-1"="#DC9F0A","Complex"="#37937F","Others"="#2C418E"))
dev.off()


############################################
############IN COMMON 17 TISSUES############
############################################
#Preparing the tissue data#
#Human#
Results <- array(NA,dim = c(length(human_exp$Name),length(tissues)))
colnames(Results)<-tissues
rownames(Results)<-human_exp$Name

for(i in 1:length(tissues)){
  
  a<-Human_experiment[Human_experiment$Tissue_new==tissues[i],]
  sampleID<-a$SampleID
  sampleID<-intersect(colnames(human_exp),sampleID)
  Result<-select(human_exp, c(Name, sampleID))
  rownames(Result)<-Result$Name
  Result<-Result[,-1]
  Result$median_value<-apply(Result,1,median)
  Result<-select(Result,c(median_value))
  Results[,i]<-Result$median_value
}
rownames(Results)<-substr(rownames(Results),1,15)
annotation$length<-annotation$`Gene end (bp)` - annotation$`Gene start (bp)` + 1

Human_ortho<-NULL
for(i in 1:length(tissues)){
  summary <- array(NA,dim = c(length(human_exp$Name),2))
  colnames(summary)<-c('Gene','median')
  summary<-as.data.frame(summary)
  summary$Gene<-rownames(Results)
  summary$median<-Results[,tissues[i]]
  summary$Tissue<-tissues[i]
  summary$Species<-'Human'
  summary$length<-annotation$length[match(summary$Gene,annotation$`Gene stable ID`)]
  summary$orthology<-annotation$`Pig homology type`[match(summary$Gene,annotation$`Gene stable ID`)]
  summary$orthology[which(summary$orthology=='ortholog_one2many'|summary$orthology=='ortholog_many2many')]<-'Complex'
  summary$orthology[which(summary$orthology!='ortholog_one2one'&summary$orthology!='Complex')]<-'Others'
  summary$reads<-summary$median * summary$length
  summary<-na.omit(summary)
  Human_ortho<-rbind(Human_ortho,summary)
}

#Pig#
Results_pig <- array(NA,dim = c(length(rownames(pig_exp)),length(tissues)))
colnames(Results_pig)<-tissues
rownames(Results_pig)<-rownames(pig_exp)

for(i in 1:length(tissues)){
  a<-Pig_experiment[Pig_experiment$Tissue_new==tissues[i],]
  sampleID<-a$SampleID
  sampleID<-intersect(colnames(pig_exp),sampleID)
  Result<-select(pig_exp, c(sampleID))
  Result$median_value<-apply(Result,1,median)
  Result<-select(Result,c('median_value'))
  Results_pig[,i]<-Result$median_value
} 
pig_anno$length<- pig_anno$`Gene end (bp)` - pig_anno$`Gene start (bp)` + 1

Pig_ortho<-NULL
for(i in 1:length(tissues)){
  summary <- array(NA,dim = c(length(rownames(Results_pig)),2))
  colnames(summary)<-c('Gene','median')
  summary<-as.data.frame(summary)
  summary$Gene<-rownames(Results_pig)
  summary$median<-Results_pig[,tissues[i]]
  summary$Tissue<-tissues[i]
  summary$Species<-'Pig'
  summary$length<-pig_anno$length[match(summary$Gene,pig_anno$`Gene stable ID`)]
  summary$orthology<-annotation$`Pig homology type`[match(summary$Gene,annotation$`Pig gene stable ID`)]
  summary$orthology[which(summary$orthology=='ortholog_one2many'|summary$orthology=='ortholog_many2many')]<-'Complex'
  summary$orthology[which(is.na(summary$orthology)==TRUE)]<-"Others"
  summary$orthology[which(summary$orthology!='ortholog_one2one'&summary$orthology!='Complex')]<-'Others'
  summary$reads<-summary$median * summary$length
  summary<-na.omit(summary)
  Pig_ortho<-rbind(Pig_ortho,summary)
}

#Integrate two species#
ortho<-rbind(Human_ortho,Pig_ortho)
ortho$orthology[which(ortho$orthology=='ortholog_one2one')]<-"1-1"
###############################################################################################

#Plot the figures#

##figures according to the reads##

human_data<-ortho[which(ortho$Species=='Human'),]
pig_data<-ortho[which(ortho$Species=='Pig'),]
p1<-ggplot(human_data,mapping = aes(x=Tissue,y=reads,fill=orthology))+
  geom_bar(stat='identity',position="fill") +
  #guides(fill = guide_legend(reverse = F))+
  labs(x = 'Human',y = 'Proportion of reads') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=0))+
  theme(axis.text.y=element_text(size=12))+
  theme(axis.title.x = element_text(vjust=-1))+
  theme(axis.title.y=element_text(size=13))+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("1-1"="#DC9F0A","Complex"="#37937F","Others"="#2C418E"))+
  annotate("point",x=1:17,y=-0.1,shape=21,fill=col,col=col,size=1.8)+
  coord_cartesian(ylim=c(0,1),clip='off')+
  theme(plot.margin=unit(c(0.4,0.4,0.5,0.4),"lines"))
p2<-ggplot(pig_data,mapping = aes(x=Tissue,y=reads,fill=orthology))+
  geom_bar(stat='identity',position="fill") +
  #guides(fill = guide_legend(reverse = F))+
  labs(x = 'Pig',y = NULL) +
  theme(axis.title =element_text(size = 12),axis.text =element_text(size = 6, color = 'black'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=0))+
  theme(axis.text.y=element_blank())+
  theme(axis.title.x = element_text(vjust=-1))+
  scale_fill_manual(values = c("1-1"="#DC9F0A","Complex"="#37937F","Others"="#2C418E"))+
  theme(legend.title = element_text(face='bold',size=12))+
  annotate("point",x=1:17,y=-0.1,shape=21,fill=col,col=col,size=1.8)+
  coord_cartesian(ylim=c(0,1),clip='off')+
  theme(plot.margin=unit(c(0.4,0.4,0.5,0.4),"lines"))+
  theme(legend.position = 'top')

p3<-p1+p2
ggsave(p3,file = "/Users/baizhonghao/Downloads/testorthologous(1-1) contribute to transcriptome in 17 common tissues_test.pdf",
     dpi = 300, width = 7, height = 2.2)

p1<-ggplot(human_data,mapping = aes(x=Tissue,y=reads,fill=Tissue))+
  geom_bar(stat='identity',position="fill") +
  #guides(fill = guide_legend(reverse = F))+
  labs(x = 'Human',y = 'Proportion of reads') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=0))+
  theme(axis.text.y=element_text(size=12))+
  theme(axis.title.x = element_text(vjust=-1))+
  theme(axis.title.y=element_text(size=13))+
  theme(legend.position = 'bottom',legend.box='horizontal')+
  theme(legend.text=element_text(size=14))+
  scale_fill_manual(values = col,guide=guide_legend(title.position = 'top',nrow=6))
  
ggsave(p1,file = "/Users/baizhonghao/Downloads/testorthologous(1-1) contribute to transcriptome in 17 common tissues_legend_4row.pdf",
       dpi = 300, width = 16, height = 2.2)

      tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/orthologous(1-1) contribute to transcriptome in 17 common tissues.tiff",
     res = 300, width = 2500, height = 1000,compression = "lzw")
ggplot(ortho,mapping = aes(x=Tissue,y=reads,fill=orthology))+
  geom_bar(stat='identity',position="fill") +
  facet_wrap(~Species, scales = 'free_x', ncol = 3) +
  guides(fill = guide_legend(reverse = F))+
  labs(x = 'Tissues',y = 'Proportion of reads') +
  theme(axis.title =element_text(size = 12),axis.text =element_text(size = 6, color = 'black'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c("1-1"="#DC9F0A","Complex"="#37937F","Others"="#2C418E"))
dev.off()

#correlation#
one2one_ortho<-array(NA,dim=c(17,2))
colnames(one2one_ortho)<-c("Human","Pig")
rownames(one2one_ortho)<-tissues
one2one_ortho<-as.data.frame(one2one_ortho)
for(i in 1:length(tissues)){
  Hum_tissuesum<-sum(Human_ortho$reads[which(Human_ortho$Tissue==tissues[i]&Human_ortho$orthology=="ortholog_one2one")])
  allsum_Hum<-sum(Human_ortho$reads[which(Human_ortho$Tissue==tissues[i])])
  one2one_ortho[i,1]<-Hum_tissuesum/allsum_Hum
  Pig_tissuesum<-sum(Pig_ortho$reads[which(Pig_ortho$Tissue==tissues[i]&Pig_ortho$orthology=="ortholog_one2one")])
  allsum_Pig<-sum(Pig_ortho$reads[which(Pig_ortho$Tissue==tissues[i])])
  one2one_ortho[i,2]<-Pig_tissuesum/allsum_Pig
}
cor<-cor(one2one_ortho$Human,one2one_ortho$Pig,method="pearson")
pval<-cor.test(one2one_ortho$Human,one2one_ortho$Pig,method="pearson")
pval<-pval$p.value
one2one_ortho$Tissues<-rownames(one2one_ortho)
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/1-1ortho_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(one2one_ortho,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("1-1",subtitle =paste0("spearman cor=",cor,"  p-value=",pval$p.value))
dev.off()

cor<-round(cor,4)
pval<-round(pval,4)
p3<-ggplot(one2one_ortho,aes(x=Human,y=Pig))+
  geom_point(size=4,aes(colour=Tissues))+
  stat_smooth(method="lm")+
  theme_classic() +
  labs(x="Proportion of reads in human", 
       y="Proportion of reads in pig") + 
  #theme(axis.title =element_text(size = 12),axis.text =element_text(size = 6, color = 'black'))+
  #theme(axis.text.x = element_text(hjust = 1,size=5.2))+
  theme(legend.position = "none")+
  #theme(axis.title.y.left = element_text(vjust = 2))+
  #theme(axis.text.x = element_text(size = 16,hjust = 0.5))+
  #theme(axis.text.y = element_text(size = 16))+
  #theme(axis.title.x = element_text(size = 18))+
  #theme(axis.title.y = element_text(size = 18))+
  #theme(axis.text=element_text(color="black"))+
  #theme(axis.title=element_text(face = "bold", color="black"))+
  theme(axis.title = element_text(color = "black", size = unit(7, "pt")),axis.text = element_text(color = "black", size = unit(7, "pt")))+
  theme(legend.title = element_text(size = unit(7, "pt")), legend.text = element_text(size = unit(6, "pt")))+
  scale_color_manual(values=col)
  #annotate("text",x=0.83,y=1,label=paste0("R=",cor,", p=",pval),size=6)
#annotate("text",x=12.3,y=15,label=paste0(italic_R,"=",corr,", ",italic_p,"=",p.val$p.value),size=6)
ggsave(p3,file = "/Users/baizhonghao/Downloads/1-1orthologous genes(0.75,5.2x10-4).pdf",##reqiured to change
     dpi= 300, width = 45, height = 45, units = 'mm')


for(i in 1:length(tissues)){
  Hum_tissuesum<-sum(Human_ortho$reads[which(Human_ortho$Tissue==tissues[i]&Human_ortho$orthology=="Complex")])
  allsum_Hum<-sum(Human_ortho$reads[which(Human_ortho$Tissue==tissues[i])])
  one2one_ortho[i,1]<-Hum_tissuesum/allsum_Hum
  Pig_tissuesum<-sum(Pig_ortho$reads[which(Pig_ortho$Tissue==tissues[i]&Pig_ortho$orthology=="Complex")])
  allsum_Pig<-sum(Pig_ortho$reads[which(Pig_ortho$Tissue==tissues[i])])
  one2one_ortho[i,2]<-Pig_tissuesum/allsum_Pig
}

cor<-cor(one2one_ortho$Human,one2one_ortho$Pig,method="spearman")
pval<-cor.test(one2one_ortho$Human,one2one_ortho$Pig,method="spearman")
pval<-pval$p.value
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/complexortho_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(one2one_ortho,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("Complex",subtitle =paste0("spearman cor=",cor,"  p-value=",pval$p.value))
dev.off()