rm(list = ls()) 
load("/Users/baizhonghao/Downloads/human-pig GTEx/DataMatrix22761.RData")
# BiocManager::install("limma")
library(limma)

###########human
Human_experiment<-experiment[which(experiment$Species=="Human"),]
tissue_name<-unique(sort(Human_experiment$Tissue_new))
tissue_name
Meta_data_human_sort<-Human_experiment[order(Human_experiment$Tissue_new),]
Meta_data_human_sort
expression[1:10,1:10]
Human_expression<-expression[match(Human_experiment$SampleID,rownames(expression)),]
human_expression<-t(Human_expression)
human_expression[1:10,1:10]
human_expression_sort<-human_expression[,match(Meta_data_human_sort$SampleID,colnames(human_expression))]
human_expression_sort[1:10,1:10]
table(colnames(human_expression_sort)==Meta_data_human_sort$SampleID)



dim(Meta_data_human_sort)

###########Limma
for (i in 1:length(tissue_name)){
  tissue_index<-which(Meta_data_human_sort$Tissue_new==tissue_name[i])
  tissue_index
  df<-c(rep("others",tissue_index[1]-1),rep("tissue",length(tissue_index)),rep("others",dim(Meta_data_human_sort)[1]-length(tissue_index)-(tissue_index[1]-1)))
  table(df)
  
  df
  design <- model.matrix(~-1+factor(df))
  dim(design)
  colnames(design) <- c("Others","Tissue")
  table(design)
  rownames(design)=colnames(human_expression_sort)
  design
  design[tissue_index[1],]
  contrastmatrix <- makeContrasts(Tissue-Others, levels=design)
  fit <- lmFit(human_expression_sort, design) ##issue these commands to fit the model and make the contrasts
  fit2 <- contrasts.fit(fit, contrastmatrix)
  fit2 <- eBayes(fit2)
  fit2
  myresults <-topTable(fit2,coef=1, adjust="fdr", number=nrow(human_expression_sort))
  write.table(myresults,file=paste0("~/Downloads/human-pig GTEx/tissue-specific limma/human/",tissue_name[i],".txt"),quote = F)
}

########pig

Pig_experiment<-experiment[which(experiment$Species=="Pig"),]
tissue_name<-unique(sort(Pig_experiment$Tissue_new))
tissue_name
Meta_data_pig_sort<-Pig_experiment[order(Pig_experiment$Tissue_new),]
Meta_data_pig_sort
Pig_expression<-expression[match(Pig_experiment$SampleID,rownames(expression)),]
pig_expression<-t(Pig_expression)
pig_expression[1:10,1:10]
pig_expression_sort<-pig_expression[,match(Meta_data_pig_sort$SampleID,colnames(pig_expression))]
pig_expression_sort[1:10,1:10]
table(colnames(pig_expression_sort)==Meta_data_pig_sort$SampleID)



dim(Meta_data_pig_sort)
i=1
###########Limma
for (i in 1:length(tissue_name)){
  
  tissue_index<-which(Meta_data_pig_sort$Tissue_new==tissue_name[i])
  tissue_index
  df<-c(rep("others",tissue_index[1]-1),rep("tissue",length(tissue_index)),rep("others",dim(Meta_data_pig_sort)[1]-length(tissue_index)-(tissue_index[1]-1)))
  table(df)
  
  df
  design <- model.matrix(~-1+factor(df))
  dim(design)
  colnames(design) <- c("Others","Tissue")
  table(design)
  rownames(design)=colnames(pig_expression_sort)
  design
  design[tissue_index[1],]
  contrastmatrix <- makeContrasts(Tissue-Others, levels=design)
  fit <- lmFit(pig_expression_sort, design) ##issue these commands to fit the model and make the contrasts
  fit2 <- contrasts.fit(fit, contrastmatrix)
  fit2 <- eBayes(fit2)
  fit2
  myresults <-topTable(fit2,coef=1, adjust="fdr", number=nrow(pig_expression_sort))
  write.table(myresults,file=paste0("~/Downloads/human-pig GTEx/tissue-specific limma/pig/",tissue_name[i],".txt"),quote = F)
}

