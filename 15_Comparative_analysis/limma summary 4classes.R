i=1
#load data#
load("/Users/baizhonghao/Downloads/human-pig GTEx/DataMatrix20217.RData")

Human_experiment<-experiment[which(experiment$Species=="Human"),]
tissue_name<-unique(sort(Human_experiment$Tissue_new))
tissue_name
Pig_experiment<-experiment[which(experiment$Species=="Pig"),]
pig_tissue<-unique(sort(Pig_experiment$Tissue_new))
pig_tissue
tissues<-intersect(tissue_name,pig_tissue)
tissues


summary<-array(NA,dim=c(17,8))
colnames(summary)<-c('human-hum','pig-hum','both-hum','neither-hum','human-pig','pig-pig','both-pig','neither-pig')
rownames(summary)<-tissues
for(i in 1:length(tissues)){
  a_tissue<-read.table(file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/species-specific limma/",tissues[i],".txt"))
  a_human<-read.table(paste0("/Users/baizhonghao/Downloads/human-pig GTEx/egene list/human_one2one_egenes_",tissues[i],".txt"))
  a_pig<-read.table(paste0("/Users/baizhonghao/Downloads/human-pig GTEx/egene list/pig_one2one_egenes_",tissues[i],".txt"))
  a_both<-read.table(paste0("/Users/baizhonghao/Downloads/human-pig GTEx/egene list/overlap_one2one_egenes_",tissues[i],".txt"))
  a_neither<-read.table(paste0("/Users/baizhonghao/Downloads/human-pig GTEx/egene list/neither_one2one_",tissues[i],".txt"))
  humanonly<-setdiff(a_human$V1,intersect(a_human$V1,a_both$V1))
  pigonly<-setdiff(a_pig$V1,intersect(a_pig$V1,a_both$V1))
  a_human<-as.data.frame(a_human[match(humanonly,a_human$V1),])
  a_pig<-as.data.frame(a_pig[match(pigonly,a_pig$V1),])
  colnames(a_human)<-'V1'
  colnames(a_pig)<-'V1'
  
  
  a_tissue$Category[which(a_tissue$adj.P.Val>0.05)]<-"not differencial expressed"
  a_tissue$Category[which(a_tissue$logFC>0 & a_tissue$adj.P.Val<=0.05)]<-"higher expression in human"
  a_tissue$Category[which(a_tissue$logFC<0 & a_tissue$adj.P.Val<=0.05)]<-"higher expression in pig"
  human_only<-a_tissue[match(a_human$V1,rownames(a_tissue)),]
  pig_only<-a_tissue[match(a_pig$V1,rownames(a_tissue)),]
  both<-a_tissue[match(a_both$V1,rownames(a_tissue)),]
  neither<-a_tissue[match(a_neither$V1,rownames(a_tissue)),]
  
  summary[i,1]<-nrow(human_only[which(human_only$Category=="higher expression in human"),])
  summary[i,2]<-nrow(pig_only[which(pig_only$Category=="higher expression in human"),])
  summary[i,3]<-nrow(both[which(both$Category=="higher expression in human"),])
  summary[i,4]<-nrow(neither[which(neither$Category=="higher expression in human"),])
  summary[i,5]<-nrow(human_only[which(human_only$Category=="higher expression in pig"),])
  summary[i,6]<-nrow(pig_only[which(pig_only$Category=="higher expression in pig"),])
  summary[i,7]<-nrow(both[which(both$Category=="higher expression in pig"),])
  summary[i,8]<-nrow(neither[which(neither$Category=="higher expression in pig"),])
  
}


logFC_sum<-NULL
adjp_sum<-NULL
i=1
for(i in 1:length(tissues)){
  a_tissue<-read.table(file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/species-specific limma/origin_",tissues[i],".txt"))
  a_human<-read.table(paste0("/Users/baizhonghao/Downloads/human-pig GTEx/egene list/human_one2one_egenes_",tissues[i],".txt"))
  a_pig<-read.table(paste0("/Users/baizhonghao/Downloads/human-pig GTEx/egene list/pig_one2one_egenes_",tissues[i],".txt"))
  a_both<-read.table(paste0("/Users/baizhonghao/Downloads/human-pig GTEx/egene list/overlap_one2one_egenes_",tissues[i],".txt"))
  a_neither<-read.table(paste0("/Users/baizhonghao/Downloads/human-pig GTEx/egene list/neither_one2one_",tissues[i],".txt"))
  humanonly<-setdiff(a_human$V1,intersect(a_human$V1,a_both$V1))
  pigonly<-setdiff(a_pig$V1,intersect(a_pig$V1,a_both$V1))
  a_human<-as.data.frame(a_human[match(humanonly,a_human$V1),])
  a_pig<-as.data.frame(a_pig[match(pigonly,a_pig$V1),])
  colnames(a_human)<-'V1'
  colnames(a_pig)<-'V1'
  
  human_only<-a_tissue[match(a_human$V1,rownames(a_tissue)),]
  pig_only<-a_tissue[match(a_pig$V1,rownames(a_tissue)),]
  both<-a_tissue[match(a_both$V1,rownames(a_tissue)),]
  neither<-a_tissue[match(a_neither$V1,rownames(a_tissue)),]
  
  a_logfc<-NULL
  tmp_logfc<-as.data.frame(human_only$logFC)
  colnames(tmp_logfc)<-'logFC'
  tmp_logfc$Category<-'human'
  tmp_logfc$Tissues<-tissues[i]
  a_logfc<-rbind(a_logfc,tmp_logfc)
  tmp_logfc<-as.data.frame(pig_only$logFC)
  colnames(tmp_logfc)<-'logFC'
  tmp_logfc$Category<-'pig'
  tmp_logfc$Tissues<-tissues[i]
  a_logfc<-rbind(a_logfc,tmp_logfc)
  tmp_logfc<-as.data.frame(both$logFC)
  colnames(tmp_logfc)<-'logFC'
  tmp_logfc$Category<-'both'
  tmp_logfc$Tissues<-tissues[i]
  a_logfc<-rbind(a_logfc,tmp_logfc)
  tmp_logfc<-as.data.frame(neither$logFC)
  colnames(tmp_logfc)<-'logFC'
  tmp_logfc$Category<-'neither'
  tmp_logfc$Tissues<-tissues[i]
  a_logfc<-rbind(a_logfc,tmp_logfc)
  logFC_sum<-rbind(logFC_sum,a_logfc)
  
  a_adjp<-NULL
  tmp_adjp<-as.data.frame(human_only$adj.P.Val)
  colnames(tmp_adjp)<-'adjp'
  tmp_adjp$Category<-'human'
  tmp_adjp$Tissues<-tissues[i]
  a_adjp<-rbind(a_adjp,tmp_adjp)
  tmp_adjp<-as.data.frame(pig_only$adj.P.Val)
  colnames(tmp_adjp)<-'adjp'
  tmp_adjp$Category<-'pig'
  tmp_adjp$Tissues<-tissues[i]
  a_adjp<-rbind(a_adjp,tmp_adjp)
  tmp_adjp<-as.data.frame(both$adj.P.Val)
  colnames(tmp_adjp)<-'adjp'
  tmp_adjp$Category<-'both'
  tmp_adjp$Tissues<-tissues[i]
  a_adjp<-rbind(a_adjp,tmp_adjp)
  tmp_adjp<-as.data.frame(neither$adj.P.Val)
  colnames(tmp_adjp)<-'adjp'
  tmp_adjp$Category<-'neither'
  tmp_adjp$Tissues<-tissues[i]
  a_adjp<-rbind(a_adjp,tmp_adjp)
  adjp_sum<-rbind(adjp_sum,a_adjp)
}

logFC_sum$Category<-factor(logFC_sum$Category,levels = c('neither','human','pig','both'))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/limma-logfc-origin.tiff",##reqiured to change
     res = 300, width =1800, height = 1800,compression = "lzw")

ggplot(logFC_sum,mapping = aes(x=Category,y=logFC,fill=Category))+
  geom_violin(aes(fill = Category), trim = FALSE) +
  geom_boxplot(width = 0.15)+
  #geom_signif(comparisons = list(c("neither","human"),
  #                               c("neither","pig"),
  #                               c("neither","both"),
  #                               c("human","pig"),
  #                               c("human","both"),
  #                               c("pig","both")),
  #            y_position = c(5.5,6.3,7.1,7.9,9.5,8.7),
  #            map_signif_level = TRUE)+
  theme_bw()+
  labs(x = 'Categories',y = 'logFC') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title =element_text(size = 12),axis.text =element_text(size = 6, color = 'black'))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1,size=5.2))+
  theme(axis.title.y.left = element_text(vjust = 2))+
  theme(axis.text.x = element_text(size=16))+
  theme(axis.text.y = element_text(size=16))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text=element_text(color="black"))+
  theme(axis.title=element_text(face = "bold", color="black"))+
  geom_hline(aes(yintercept=0),linetype='dashed',color='black')
dev.off()
adjp_sum$adjp<- -log10(adjp_sum$adjp+0.001)
adjp_sum$Category<-factor(adjp_sum$Category,levels = c("neither","human","pig","both"))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/limma-adjp-origin.tiff",##reqiured to change
     res = 300, width =1800, height = 1800,compression = "lzw")

ggplot(adjp_sum,mapping = aes(x=Category,y=adjp,fill=Category))+
  geom_violin(aes(fill = Category), trim = FALSE) +
  geom_boxplot(width = 0.15,outlier.colour = NA)+
  #geom_signif(comparisons = list(c("neither","human"),
  #                               c("neither","pig"),
  #                               c("neither","both"),
  ##                               c("human","pig"),
  #                               c("human","both"),
  #                               c("pig","both")),
  #            y_position = c(3.7,4.0,4.3,4.6,5.2,4.9),
  #            map_signif_level = TRUE)+
  theme_bw()+
  labs(x = 'Categories',y = '-log10(adj.pval)') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title =element_text(size = 12),axis.text =element_text(size = 6, color = 'black'))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1,size=5.2))+
  theme(axis.title.y.left = element_text(vjust = 2))+
  theme(axis.text.x = element_text(size=16))+
  theme(axis.text.y = element_text(size=16))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text=element_text(color="black"))+
  theme(axis.title=element_text(face = "bold", color="black"))
dev.off()



