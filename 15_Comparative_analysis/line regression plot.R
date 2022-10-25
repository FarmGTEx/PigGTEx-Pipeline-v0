#load data#
load("/Users/baizhonghao/Downloads/human-pig GTEx/DataMatrix20217.RData")

Human_sample<-experiment[which(experiment$Species=="Human"),]
Pig_sample<-experiment[which(experiment$Species=="Pig"),]

expressed_genes_tissues<-array(NA, dim=c(length(tissues),2))
colnames(expressed_genes_tissues)<-c("Human","Pig")
rownames(expressed_genes_tissues)<-tissues
expressed_genes_tissues<-as.data.frame(expressed_genes_tissues)
for(i in 1:length(tissues)){
  Human_exp<-expression[match(Human_sample$SampleID[which(Human_sample$Tissue_new==tissues[i])],rownames(expression)),]
  Pig_exp<-expression[match(Pig_sample$SampleID[which(Pig_sample$Tissue_new==tissues[i])],rownames(expression)),]
  Human_exp<-as.data.frame(t(Human_exp))
  Pig_exp<-as.data.frame(t(Pig_exp))
  Human_exp$median<-apply(Human_exp,1,median)
  Pig_exp$median<-apply(Pig_exp,1,median)
  a_Hum<-Human_exp[which(Human_exp$median>0.1),]
  a_Pig<-Pig_exp[which(Pig_exp$median>0.1),]
  expressed_genes_tissues$Human[i]<-nrow(a_Hum)
  expressed_genes_tissues$Pig[i]<-nrow(a_Pig)
}
#calculate correlation#
expressed_genes_tissues<-expressed_genes_tissues/1000
corr<-cor(expressed_genes_tissues$Human,expressed_genes_tissues$Pig)
p.val<-cor.test(expressed_genes_tissues$Human,expressed_genes_tissues$Pig,method = "spearman")
corr<-round(corr,2)

expressed_genes_tissues$Tissues<-rownames(expressed_genes_tissues)


p1<-ggplot(expressed_genes_tissues,aes(x=Human,y=Pig))+
  geom_point(size=4,aes(colour=Tissues))+
  stat_smooth(method="lm")+
  theme_classic() +
  labs(x="Number of expressed \ngenes in human (x1000)", 
       y="Number of expressed \ngenes in pig (x1000)") + 
  theme(legend.position = "none")+
  theme(axis.title = element_text(color = "black", size = unit(7, "pt")),axis.text = element_text(color = "black", size = unit(7, "pt")))+
  theme(legend.title = element_text(size = unit(7, "pt")), legend.text = element_text(size = unit(6, "pt")))+
  scale_color_manual(values=col)
ggsave(p1,file = "/Users/baizhonghao/Downloads/expressed_gene_cor(0.84,p4.3x10-4).pdf",dpi=300,width=45,height=45,units = 'mm')
