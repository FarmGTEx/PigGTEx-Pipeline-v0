#!/usr/bin/env Rscript
# Step 2. Detect tissue-specific genes
library(limma)
load("Pig_expression_new.Rdata")
color <- read.xlsx("color.xlsx",sheetIndex = 1)
color <- color[color$Categories.class2.md%in%tissue_name,]
tissue_name <- unique(sort(Pig_meta_data_new$Tissue))
Pig_meta_data_new$Tissue_categories <- color[Pig_meta_data_new$Tissue,"Categories.class2"]
Pig_meta_data_new$Tissue_class <- Pig_meta_data_new$Tissue

for (i in 1:34){
  tissue_cata <- unique(Pig_meta_data_new$Tissue_categories[Pig_meta_data_new$Tissue_class==tissue_name[i]])
  tissue_cata_sample <- Pig_meta_data_new$BioSample[Pig_meta_data_new$Tissue_categories==tissue_cata]
  tissue_class_sample <- Pig_meta_data_new$BioSample[Pig_meta_data_new$Tissue_class==tissue_name[i]]  
  pig_meta_data <- Pig_meta_data_new[!Pig_meta_data_new$BioSample%in%remove_sample,]
  Meta_data_pig_sort <- pig_meta_data[order(pig_meta_data$Tissue_class),]
  pig_expression_sort <- pig_expression[,match(Meta_data_pig_sort$BioSample, colnames(pig_expression))]
  table(colnames(pig_expression_sort)==Meta_data_pig_sort$BioSample)
  pig_expression_sort <- log2(pig_expression_sort+0.25)  
  tissue_index <- which(Meta_data_pig_sort$Tissue_class==tissue_name[i])
  df <- c(rep("others",tissue_index[1]-1),rep("tissue",length(tissue_index)),rep("others",dim(Meta_data_pig_sort)[1]-length(tissue_index)-(tissue_index[1]-1)))
  design <- model.matrix(~-1+factor(df))
  colnames(design) <- c("Others","Tissue")
  rownames(design)=colnames(pig_expression_sort)
  design[tissue_index[1],]
  contrastmatrix <- makeContrasts(Tissue-Others, levels=design)
  fit <- lmFit(pig_expression_sort, design)
  fit2 <- contrasts.fit(fit, contrastmatrix)
  fit2 <- eBayes(fit2)
  myresults <-topTable(fit2,coef=1, adjust="fdr", number=nrow(pig_expression_sort))
  write.table(myresults,file=paste0(tissue_name[i],".txt"),quote = F)
}
