#violin LOEUF#

eqtlfm<-as.data.frame(fread("/Users/baizhonghao/Downloads/Pig.chrAll.gene_coverage_phastCons.txt"))
eqtlfm$PhastCons[which(eqtlfm$PhastCons==".")]<-0
eqtlfm$PhastCons<-as.numeric(eqtlfm$PhastCons)
one2one_pig<-annotation[which(annotation$`Pig homology type`=='ortholog_one2one'),]
eqtlfm<-eqtlfm[match(intersect(eqtlfm$gene_id,one2one_pig$`Pig gene stable ID`),eqtlfm$gene_id),]
eqtlfm$gene_id<-one2one_pig$`Gene stable ID`[match(eqtlfm$gene_id, one2one_pig$`Pig gene stable ID`)]
LOEUF_sum<-NULL
for ( i in 1:length(tissues)){
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
  
  a_human$name<-annotation$`Gene name`[match(a_human$V1,annotation$`Gene stable ID`)]
  a_pig$name<-annotation$`Gene name`[match(a_pig$V1,annotation$`Gene stable ID`)]
  a_both$name<-annotation$`Gene name`[match(a_both$V1,annotation$`Gene stable ID`)]
  a_neither$name<-annotation$`Gene name`[match(a_neither$V1,annotation$`Gene stable ID`)]
  
  a_human$LOEUF<-df$oe_lof_upper[match(a_human$name,df$gene)]
  a_pig$LOEUF<-df$oe_lof_upper[match(a_pig$name,df$gene)]
  a_both$LOEUF<-df$oe_lof_upper[match(a_both$name,df$gene)]
  a_neither$LOEUF<-df$oe_lof_upper[match(a_neither$name,df$gene)]
  a_human$category<-'human-specific'
  a_pig$category<-'pig-specific'
  a_both$category<-'shared'
  a_neither$category<-'neither'
  a_sum<-rbind(a_human,a_pig,a_both,a_neither)
  a_sum$Tissues<-tissues[i]
  LOEUF_sum<-rbind(LOEUF_sum,a_sum)
}

LOEUF_mean<-NULL
for ( i in 1:length(tissues)){
  a_human<-read.table(paste0("/Users/baizhonghao/Downloads/human-pig GTEx/egene list/human_one2one_egenes_",tissues[i],".txt"))
  a_pig<-read.table(paste0("/Users/baizhonghao/Downloads/human-pig GTEx/egene list/pig_one2one_egenes_",tissues[i],".txt"))
  a_both<-read.table(paste0("/Users/baizhonghao/Downloads/human-pig GTEx/egene list/overlap_one2one_egenes_",tissues[i],".txt"))
  a_neither<-read.table(paste0("/Users/baizhonghao/Downloads/human-pig GTEx/egene list/neither_one2one_",tissues[i],".txt"))
  a_mean<-array(NA, dim=c(4,3))
  colnames(a_mean)<-c("LOEUF","category","Tissues")
  a_mean<-as.data.frame(a_mean)
  humanonly<-setdiff(a_human$V1,intersect(a_human$V1,a_both$V1))
  pigonly<-setdiff(a_pig$V1,intersect(a_pig$V1,a_both$V1))
  a_human<-as.data.frame(a_human[match(humanonly,a_human$V1),])
  a_pig<-as.data.frame(a_pig[match(pigonly,a_pig$V1),])
  colnames(a_human)<-'V1'
  colnames(a_pig)<-'V1'
  
  a_human$name<-annotation$`Gene name`[match(a_human$V1,annotation$`Gene stable ID`)]
  a_pig$name<-annotation$`Gene name`[match(a_pig$V1,annotation$`Gene stable ID`)]
  a_both$name<-annotation$`Gene name`[match(a_both$V1,annotation$`Gene stable ID`)]
  a_neither$name<-annotation$`Gene name`[match(a_neither$V1,annotation$`Gene stable ID`)]
  
  a_human$LOEUF<-df$oe_lof_upper[match(a_human$name,df$gene)]
  a_pig$LOEUF<-df$oe_lof_upper[match(a_pig$name,df$gene)]
  a_both$LOEUF<-df$oe_lof_upper[match(a_both$name,df$gene)]
  a_neither$LOEUF<-df$oe_lof_upper[match(a_neither$name,df$gene)]
  
  a_human<-a_human[-which(is.na(a_human$LOEUF)=="TRUE"),]
  a_pig<-a_pig[-which(is.na(a_pig$LOEUF)=="TRUE"),]
  a_both<-a_both[-which(is.na(a_both$LOEUF)=="TRUE"),]
  a_neither<-a_neither[-which(is.na(a_neither$LOEUF)=="TRUE"),]
  
  a_mean$LOEUF[1]<-mean(a_human$LOEUF)
  a_mean$LOEUF[2]<-mean(a_pig$LOEUF)
  a_mean$LOEUF[3]<-mean(a_both$LOEUF)
  a_mean$LOEUF[4]<-mean(a_neither$LOEUF)
  a_mean$category[1]<-'human-specific'
  a_mean$category[2]<-'pig-specific'
  a_mean$category[3]<-'shared'
  a_mean$category[4]<-'neither'
  a_mean$Tissues<-tissues[i]
  LOEUF_mean<-rbind(LOEUF_mean,a_mean)
}

library(ggsignif)
LOEUF_mean$category<-factor(LOEUF_mean$category,levels = c("neither","human-specific","pig-specific","shared"))
#tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/eqtl/Phastconsviolinplot-categoryx.tiff",##reqiured to change
#     res = 300, width =2000, height = 2000,compression = "lzw")

LOEUFvio<-ggplot(LOEUF_mean,mapping = aes(x=category,y=LOEUF,fill='red'))+
  geom_violin(aes(fill = 'red'), trim = FALSE) +
  geom_boxplot(width = 0.15)+
  geom_signif(comparisons = list(c("neither","human-specific"),
                                 c("neither","pig-specific"),
                                 c("neither","shared"),
                                 #c("human_only","pig_only"),
                                 c("human-specific","shared"),
                                 c("pig-specific","shared")),
              y_position = c(1.2,1.24,1.28,1.32,1.36),
              map_signif_level = TRUE)+
  theme_classic()+
  labs(x = NULL,y = 'LOEUF') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title =element_text(size = 12),axis.text =element_text(size = 6, color = 'black'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=5.2))+
  theme(axis.title.y.left = element_text(vjust = 2))+
  theme(axis.text.x = element_text(size=16))+
  theme(axis.text.y = element_text(size=16))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text=element_text(color="black"))+
  theme(axis.title=element_text(face = "bold", color="black"))+
  theme(legend.title=element_text(face='bold',size=18))+
  theme(legend.text = element_text(size=16))+
  theme(legend.position= 'none')
#dev.off()
ggsave(LOEUFvio,file='/Users/baizhonghao/Downloads/LOEUFviolin-test.pdf',dpi=300,width=6,height=6)

t.test(LOEUF_mean$LOEUF[which(LOEUF_mean$category=='pig-specific')],LOEUF_mean$LOEUF[which(LOEUF_mean$category=='shared')])
