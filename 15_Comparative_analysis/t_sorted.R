rm(list = ls()) 
library(clusterProfiler)
library(xlsx)
library(org.Hs.eg.db)

#read Limma results#
human_file<-dir("/Users/baizhonghao/Downloads/human-pig GTEx/tissue-specific limma/human",full.names = T)
human_file
humantissue<-substr(human_file,nchar("/Users/baizhonghao/Downloads/human-pig GTEx/tissue-specific limma/human/")+1,nchar(human_file)-nchar(".txt"))
humantissue

pig_file<-dir("/Users/baizhonghao/Downloads/human-pig GTEx/tissue-specific limma/pig",full.names = T)
pig_file
pigtissue<-substr(pig_file,nchar("/Users/baizhonghao/Downloads/human-pig GTEx/tissue-specific limma/pig/")+1,nchar(pig_file)-nchar(".txt"))
pigtissue

tissue <-intersect(humantissue,pigtissue)
Results1 <- array(NA,dim = c(length(tissue),10))

i=1
rownames(Results1)<-tissue
colnames(Results1)<-c("top10%","10%-20%","20%-30%","30%-40%","40%-50%","50%-60%","60%-70%","70%-80%","80%-90%","bottom10%")
Results1
for (i in 1:length(tissue)){
  human1<- read.table(file = paste0("/Users/baizhonghao/Downloads/human-pig GTEx/tissue-specific limma/human/",tissue[i],".txt"),header = T,stringsAsFactors = F)
  pig1<- read.table(file = paste0("/Users/baizhonghao/Downloads/human-pig GTEx/tissue-specific limma/pig/",tissue[i],".txt"),header = T,stringsAsFactors = F)
  #rownames(human1)<-Orthologous_human_cattle$`Cow gene stable ID`[match(rownames(human1),Orthologous_human_cattle$`Gene stable ID`)]
  human1$t<-abs(human1$t)
  pig1$t<-abs(pig1$t)
  human1<-human1[order(human1$t,decreasing = TRUE),]
  pig1<-pig1[order(pig1$t,decreasing = TRUE),]
  
  top10_human_tissue<-rownames(human1)[1:1594]
  top10_pig_tissue<-rownames(pig1)[1:1594]
  
  tissue_inter<-intersect(top10_human_tissue,top10_pig_tissue)
  Results1[i,1]<-length(tissue_inter)
  
  top10_20_human_tissue<-rownames(human1)[1595:3188]
  top10_20_pig_tissue<-rownames(pig1)[1595:3188]
  
  
  tissue_inter_10_20<-intersect(top10_20_human_tissue,top10_20_pig_tissue)
  length(tissue_inter_10_20)
  Results1[i,2]<-length(tissue_inter_10_20)
  
  top20_30_human_tissue<-rownames(human1)[3189:4782]
  top20_30_pig_tissue<-rownames(pig1)[3189:4782]
  
  tissue_inter_20_30<-intersect(top20_30_human_tissue,top20_30_pig_tissue)
  length(tissue_inter_20_30)
  Results1[i,3]<-length(tissue_inter_20_30)
  
  top30_40_human_tissue<-rownames(human1)[4783:6376]
  top30_40_pig_tissue<-rownames(pig1)[4783:6376]
  
  tissue_inter_30_40<-intersect(top20_30_human_tissue,top30_40_pig_tissue)
  length(tissue_inter_30_40)
  Results1[i,4]<-length(tissue_inter_30_40)
  
  top40_50_human_tissue<-rownames(human1)[6377:7970]
  top40_50_pig_tissue<-rownames(pig1)[6377:7970]
  
  tissue_inter_40_50<-intersect(top40_50_human_tissue,top40_50_pig_tissue)
  length(tissue_inter_40_50)
  Results1[i,5]<-length(tissue_inter_40_50)
  
  top50_60_human_tissue<-rownames(human1)[7971:9564]
  top50_60_pig_tissue<-rownames(pig1)[7971:9564]
  
  tissue_inter_50_60<-intersect(top50_60_human_tissue,top50_60_pig_tissue)
  length(tissue_inter_50_60)
  Results1[i,6]<-length(tissue_inter_50_60)
  
  top60_70_human_tissue<-rownames(human1)[9565:11158]
  top60_70_pig_tissue<-rownames(pig1)[9565:11158]
  
  tissue_inter_60_70<-intersect(top60_70_human_tissue,top60_70_pig_tissue)
  length(tissue_inter_60_70)
  Results1[i,7]<-length(tissue_inter_60_70)
  
  top70_80_human_tissue<-rownames(human1)[11159:12752]
  top70_80_pig_tissue<-rownames(pig1)[11159:12752]
  
  tissue_inter_70_80<-intersect(top70_80_human_tissue,top70_80_pig_tissue)
  length(tissue_inter_70_80)
  Results1[i,8]<-length(tissue_inter_70_80)
  
  
  top80_90_human_tissue<-rownames(human1)[12753:14346]
  top80_90_pig_tissue<-rownames(pig1)[12753:14346]
  
  tissue_inter_80_90<-intersect(top80_90_human_tissue,top80_90_pig_tissue)
  length(tissue_inter_80_90)
  
  Results1[i,9]<-length(tissue_inter_80_90)
  
  
  last10_human_tissue<-rownames(human1)[14347:15944]
  last10_pig_tissue<-rownames(pig1)[14347:15944]
  
  tissue_inter_last10<-intersect(last10_human_tissue,last10_pig_tissue)
  length(tissue_inter_last10)
  Results1[i,10]<-length(tissue_inter_last10)
  
  #############################################################################
  
  #change the gene name#
  gene_overlap_top10<- bitr(tissue_inter, fromType="ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb=org.Hs.eg.db, drop = T)
  gene_overlap_last10<-bitr(tissue_inter_last10, fromType="ENSEMBL", toType = c("SYMBOL", "ENTREZID"), OrgDb=org.Hs.eg.db, drop = T)
  
  #############################################################################
  
  #GO enrichment#
  
  ##overlap##
  ego_BP_overlap_top10 <- enrichGO(gene = gene_overlap_top10$ENTREZID,
  OrgDb= org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  minGSSize = 1,
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE)
  ##human##
  ego_BP_overlap_last10 <- enrichGO(gene = gene_overlap_last10$ENTREZID,
  OrgDb= org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  minGSSize = 1,
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE)
  
  #############################################################################
  
  #save the data#
  #write.xlsx(ego_BP_overlap_top10,file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/t-sorted/overlap_top10_",tissue[i],".xlsx"),showNA=TRUE)
  #write.xlsx(ego_BP_overlap_last10,file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/t-sorted/overlap_last10_",tissue[i],".xlsx"),showNA=TRUE)
  
  #write.table(top10_human_tissue,file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/t-sorted/top10_human_",tissue[i],".txt"),sep=" ")
  #write.table(top10_cattle_tissue,file=paste0("/Users/baizhonghao/Downloads/top10_cattle_",cattletissue[i],".txt"),sep=" ")
  #write.table(top10_pig_tissue,file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/t-sorted/top10_pig_",tissue[i],".txt"),sep=" ")
  #write.table(top10_sheep_tissue,file=paste0("/Users/baizhonghao/Downloads/top10_sheep_",sheeptissue[i],".txt"),sep=" ")
  #write.table(top10_goat_tissue,file=paste0("/Users/baizhonghao/Downloads/top10_goat_",goattissue[i],".txt"),sep=" ")
  #write.table(top10_chicken_tissue,file=paste0("/Users/baizhonghao/Downloads/top10_chicken_",chickentissue[i],".txt"),sep=" ")
  #write.table(last10_human_tissue,file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/t-sorted/last10_human_",tissue[i],".txt"),sep=" ")
  #write.table(last10_cattle_tissue,file=paste0("/Users/baizhonghao/Downloads/last10_cattle_",cattletissue[i],".txt"),sep=" ")
  #write.table(last10_pig_tissue,file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/t-sorted/last10_pig_",tissue[i],".txt"),sep=" ")
  #write.table(last10_sheep_tissue,file=paste0("/Users/baizhonghao/Downloads/last10_sheep_",sheeptissue[i],".txt"),sep=" ")
  #write.table(last10_goat_tissue,file=paste0("/Users/baizhonghao/Downloads/last10_goat_",goattissue[i],".txt"),sep=" ")
  #write.table(last10_chicken_tissue,file=paste0("/Users/baizhonghao/Downloads/last10_chicken_",chickentissue[i],".txt"),sep=" ")
  #write.table(tissue_inter,file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/t-sorted/top10_overlap_",tissue[i],".txt"),sep=" ")
  #write.table(tissue_inter_last10,file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/t-sorted/last10_overlap_",tissue[i],".txt"),sep=" ")
}
Results1
Results1_percent<-Results1/1594
Results1_percent
col<-c("#CC66FF","#AAAAFF","#FF0000","#8EABD2","#FFD700","#33CCCC","#FDFDBF", "#8EA9DB","#8B0F55",
       "#E2EFDA","#7570B3","#AAEEFF","#99BB88","#FFDD99","#FF6600","#A6CEE3","#A6761D")
names(col)<-c("Adipose","Artery","Blood","Colon","Frontal_cortex","Heart","Hypothalamus","Ileum",
              "Kidney","Liver","Lung","Muscle","Ovary", "Pituitary","Spleen","Testis","Uterus")
pdf(file = "/Users/baizhonghao/Downloads/t-value_new.pdf",width = 19, height = 8)
par(mar=c(8,8,6,1))
my<-barplot(Results1_percent,beside = T,col=col,cex.axis = 2,cex.names = 2,ylim = c(0,0.8))
title(ylab="Proportion of shared orthologs", xlab = "Windows of genes sorted by tissue specificity",line=4,cex.lab=2.5)
dev.off()
