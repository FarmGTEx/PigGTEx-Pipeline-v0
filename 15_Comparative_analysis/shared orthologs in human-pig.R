rm(list = ls()) 
load("/Users/baizhonghao/Downloads/human-pig GTEx/DataMatrix20217.RData")

#tissue info#
Human_experiment<-experiment[which(experiment$Species=="Human"),]
tissue_name<-unique(sort(Human_experiment$Tissue_new))
tissue_name
Pig_experiment<-experiment[which(experiment$Species=="Pig"),]
pig_tissue<-unique(sort(Pig_experiment$Tissue_new))
pig_tissue
tissues<-intersect(tissue_name,pig_tissue)
tissues

Results1 <- array(NA,dim = c(length(tissues),10))
experiment$annotation<-paste(experiment$Tissue_new,experiment$Species,sep="-")
expression_tpm<-expression[match(rownames(expression),rownames(expression)),]

median_tissue_tpm<-aggregate(expression_tpm,list(as.factor(experiment$annotation)),median) ## also calculated mad/ variance with the same script
median_tissue_tpm[1:10,1:10]
rownames(median_tissue_tpm)<-median_tissue_tpm[,1]
median_tissue_tpm<-median_tissue_tpm[,-1]
median_tissue_tpm[,1:10]


rownames(Results1)<-tissues
colnames(Results1)<-c("top10%","10%-20%","20%-30%","30%-40%","40%-50%","50%-60%","60%-70%","70%-80%","80%-90%","bottom10%")

### calculated how many overlapped genes in each 10% window between human and cattle
for (i in 1:length(tissues)){

  human_name<-paste(tissues[i],"Human",sep="-")
  pig_name<-paste(tissues[i],"Pig",sep="-")
  human_tissue<-sort(median_tissue_tpm[human_name,],decreasing = T)
  pig_tissue<-sort(median_tissue_tpm[pig_name,],decreasing = T)
  
  top10_human_tissue<-colnames(human_tissue)[1:1594]
  top10_pig_tissue<-colnames(pig_tissue)[1:1594]
  
  tissue_inter<-intersect(top10_human_tissue,top10_pig_tissue)
  Results1[i,1]<-length(tissue_inter)
  
  top10_20_human_tissue<-colnames(human_tissue)[1595:3188]
  top10_20_pig_tissue<-colnames(pig_tissue)[1595:3188]

  
  tissue_inter_10_20<-intersect(top10_20_human_tissue,top10_20_pig_tissue)
  length(tissue_inter_10_20)
  Results1[i,2]<-length(tissue_inter_10_20)
  
  top20_30_human_tissue<-colnames(human_tissue)[3189:4782]
  top20_30_pig_tissue<-colnames(pig_tissue)[3189:4782]
  
  tissue_inter_20_30<-intersect(top20_30_human_tissue,top20_30_pig_tissue)
  length(tissue_inter_20_30)
  Results1[i,3]<-length(tissue_inter_20_30)
  
  top30_40_human_tissue<-colnames(human_tissue)[4783:6376]
  top30_40_pig_tissue<-colnames(pig_tissue)[4783:6376]
  
  tissue_inter_30_40<-intersect(top20_30_human_tissue,top30_40_pig_tissue)
  length(tissue_inter_30_40)
  Results1[i,4]<-length(tissue_inter_30_40)
  
  top40_50_human_tissue<-colnames(human_tissue)[6377:7970]
  top40_50_pig_tissue<-colnames(pig_tissue)[6377:7970]

  tissue_inter_40_50<-intersect(top40_50_human_tissue,top40_50_pig_tissue)
  length(tissue_inter_40_50)
  Results1[i,5]<-length(tissue_inter_40_50)
  
  top50_60_human_tissue<-colnames(human_tissue)[7971:9564]
  top50_60_pig_tissue<-colnames(pig_tissue)[7971:9564]
  
  tissue_inter_50_60<-intersect(top50_60_human_tissue,top50_60_pig_tissue)
  length(tissue_inter_50_60)
  Results1[i,6]<-length(tissue_inter_50_60)

  top60_70_human_tissue<-colnames(human_tissue)[9565:11158]
  top60_70_pig_tissue<-colnames(pig_tissue)[9565:11158]
  
  tissue_inter_60_70<-intersect(top60_70_human_tissue,top60_70_pig_tissue)
  length(tissue_inter_60_70)
  Results1[i,7]<-length(tissue_inter_60_70)
  
  top70_80_human_tissue<-colnames(human_tissue)[11159:12752]
  top70_80_pig_tissue<-colnames(pig_tissue)[11159:12752]
  
  tissue_inter_70_80<-intersect(top70_80_human_tissue,top70_80_pig_tissue)
  length(tissue_inter_70_80)
  Results1[i,8]<-length(tissue_inter_70_80)
  
  
  top80_90_human_tissue<-colnames(human_tissue)[12753:14346]
  top80_90_pig_tissue<-colnames(pig_tissue)[12753:14346]
  
  tissue_inter_80_90<-intersect(top80_90_human_tissue,top80_90_pig_tissue)
  length(tissue_inter_80_90)
  
  Results1[i,9]<-length(tissue_inter_80_90)
  
  
  last10_human_tissue<-colnames(human_tissue)[14347:15944]
  last10_pig_tissue<-colnames(pig_tissue)[14347:15944]
 
  tissue_inter_last10<-intersect(last10_human_tissue,last10_pig_tissue)
  length(tissue_inter_last10)
  Results1[i,10]<-length(tissue_inter_last10)
  
}
Results1
Results1_percent<-Results1/1594
Results1_percent
col<-c("#CC66FF","#AAAAFF","#FF0000","#8EABD2","#FFD700","#33CCCC","#FDFDBF", "#8EA9DB","#8B0F55",
       "#E2EFDA","#7570B3","#AAEEFF","#99BB88","#FFDD99","#FF6600","#A6CEE3","#A6761D")
names(col)<-c("Adipose","Artery","Blood","Colon","Frontal_cortex","Heart","Hypothalamus","Ileum",
              "Kidney","Liver","Lung","Muscle","Ovary", "Pituitary","Spleen","Testis","Uterus")
pdf(file = "~/Downloads/median_sorted.pdf", width = 18, height = 8)
par(mar=c(8,8,6,1))
my<-barplot(Results1_percent,beside = T,col=col,cex.axis = 1.8,cex.names = 1.8,ylim = c(0,0.7))
title(ylab="Proportion of shared orthologs", xlab = "Windows of genes sorted by median expression",line=4, cex.lab=2.5)
dev.off()




