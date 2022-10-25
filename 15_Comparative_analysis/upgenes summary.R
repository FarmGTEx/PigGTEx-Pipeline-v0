
rm(list = ls())

library(ggplot2)
##########

Human_experiment<-experiment[which(experiment$Species=="Human"),]
tissue_name<-unique(sort(Human_experiment$Tissue_new))
tissue_name
Pig_experiment<-experiment[which(experiment$Species=="Pig"),]
pig_tissue<-unique(sort(Pig_experiment$Tissue_new))
pig_tissue
tissues<-intersect(tissue_name,pig_tissue)
tissues

Results<-array(data = NA,dim = c(17,3))
rownames(Results)<-tissues
colnames(Results)<-c("Human","Pig","overlap")
tissue_name<-tissues

logFC_cutoff<-1.5
for(i in 1:17){
  tissue_human<-read.table(file = paste0("/Users/baizhonghao/Downloads/human-pig GTEx/tissue-specific limma/human/",tissue_name[i],".txt"))
  tissue_pig<-read.table(file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/tissue-specific limma/pig/",tissue_name[i],".txt"))

  tissue_human$change = as.factor(ifelse(tissue_human$adj.P.Val < 0.05 & abs(tissue_human$logFC) > logFC_cutoff,
                                         ifelse(tissue_human$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  
  human_upgene<-tissue_human[tissue_human$change =='UP',]
  #write.table(human_upgene,file=paste0("~/Downloads/human-tissue-specific-gene/",tissue_name[i],".txt"))
  Results[i,1]<- length(rownames(human_upgene)) 
  
  tissue_pig$change = as.factor(ifelse(tissue_pig$adj.P.Val < 0.05 & abs(tissue_pig$logFC) > logFC_cutoff,
                                       ifelse(tissue_pig$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
  )
  
  pig_upgene<-tissue_pig[tissue_pig$change =='UP',]
  #write.table(pig_upgene,file=paste0("~/Downloads/pig-tissue-specific-gene/",tissue_name[i],".txt"))
  Results[i,2]<- length(rownames(pig_upgene))    

  overlap<-intersect(rownames(human_upgene),rownames(pig_upgene))
  Results[i,3]<-length(overlap)
  write.table(overlap,file=paste0("~/Downloads/overlap-tissue-specific-gene/",tissue_name[i],".txt"))
  #############################################################################
  
  #change the gene name#
  tsoverlap<- bitr(overlap, fromType="ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb=org.Hs.eg.db, drop = T)
  
  #############################################################################
  
  #GO enrichment#
  
  ##overlap##
  ego_BP_overlap<- enrichGO(gene = tsoverlap$ENTREZID,
                                   OrgDb= org.Hs.eg.db,
                                   ont = "BP",
                                   pAdjustMethod = "BH",
                                   minGSSize = 1,
                                   pvalueCutoff = 1,
                                   qvalueCutoff = 1,
                                   readable = TRUE)
  
  write.xlsx(ego_BP_overlap,file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/tsoverlap_",tissue[i],".xlsx"),showNA=TRUE)
  write.table(overlap,file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/overlapgene_",tissue[i],".txt"),sep=" ")
}
library(data.table)

#calculate the corraltion of different tissues#
summary<-melt(Results)
colnames(summary)<-c('Tissues','Species','number')
summary<-as.data.frame(summary)
p.val1<-cor.test(Results[,1],Results[,3])
p.val2<-cor.test(Results[,2],Results[,3])
p.val1<-'1.76e-08'
p.val2<-'5.08e-10'
#tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/limma_summary.tiff",##reqiured to change
#     res = 300, width = 4000, height = 2000,compression = "lzw")
par(mar=c(8,8,6,1))
p1<-ggplot(summary,mapping = aes(x=Tissues,y=number,fill=Species))+
  theme_classic()+
  geom_bar(stat='identity',position="dodge")+  expand_limits(x = 0, y = 0)+
  labs(x = 'Tissues',y = 'Number of tissue-specific genes') +
  theme(axis.title =element_text(size = 12),axis.text =element_text(size = 6, color = 'black'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=5.2))+
  theme(axis.title.y.left = element_text(vjust = 2))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text=element_text(color="black"))+
  theme(axis.title=element_text(face = "bold", color="black"))+
  scale_fill_manual(values = c('Human'="Blue",'Pig'="Red",'overlap'="lightgrey"))+
  annotate("text",x=13,y=450,label=paste0("Human-overlap, p=",p.val1,"\nPig-overlap, p=",p.val2),size=4)+
  theme(legend.title = element_text(face="bold"))
#dev.off()
ggsave(p1,file='/Users/baizhonghao/Downloads/limma_summary.pdf',dpi=300,width=10,height=6)

corr<-cor(Results[,1],Results[,2])
corr<-round(corr,4)
p.val<-cor.test(Results[,1],Results[,2])
p.val<-'2.83e-06'
Results<-as.data.frame(Results)
Results$Tissues<-rownames(Results)
tiff(file = "/Users/baizhonghao/Downloads/tissue_specific_gene_cor.tiff",##reqiured to change
     res = 300, width = 1300, height = 1300,compression = "lzw")
p1<-ggplot(Results,aes(x=Human,y=Pig))+
  geom_point(size=4,aes(colour=Tissues))+
  stat_smooth(method="lm")+
  theme_classic() +
  labs(x="Number of tissue-specific \n genes in human", 
       y="Number of tissue-specific \n genes in pig") + 
  #theme(axis.title =element_text(size = 12),axis.text =element_text(size = 6, color = 'black'))+
  #theme(axis.text.x = element_text(hjust = 0.5,size=5.2))+
  theme(legend.position = "none")+
  #theme(axis.title.y.left = element_text(vjust = 2))+
  #theme(axis.text.x = element_text(size = 16))+
  #theme(axis.text.y = element_text(size = 16))+
  #theme(axis.title.x = element_text(size = 18))+
  #theme(axis.title.y = element_text(size = 18))+
  #theme(axis.text=element_text(color="black"))+
  #theme(axis.title=element_text(face = "bold", color="black"))+
  theme(axis.title = element_text(color = "black", size = unit(7, "pt")),axis.text = element_text(color = "black", size = unit(7, "pt")))+
  theme(legend.title = element_text(size = unit(7, "pt")), legend.text = element_text(size = unit(6, "pt")))+
  scale_color_manual(values=col)
  #annotate("text",x=200,y=550,label=paste0("R=",corr,", p=",p.val),size=4)
#annotate("text",x=12.3,y=15,label=paste0(italic_R,"=",corr,", ",italic_p,"=",p.val$p.value),size=6)
#dev.off()
ggsave(p1,file='/Users/baizhonghao/Downloads/tissue-specific(0.88,2.8x10-6).pdf',dpi = 300,width = 45,height = 45,units = 'mm')
