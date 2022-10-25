library(data.table)
library(readxl)
library(xlsx)

#load data#
human_meta<-as.data.frame(fread("/Users/baizhonghao/Downloads/Comparative/human/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"))
pig_meta <- read_excel("/Users/baizhonghao/Downloads/human-pig GTEx/Metadata_eQTLmapping_20210702_tjy.xlsx",sheet = "Metatable3_rmDup", col_names = T)

######################################################################
#integrate pigs' and humans' expression data#
pig_file<-dir("/Users/baizhonghao/Downloads/human-pig GTEx/SampleID_each_tissue",full.names = T)
pig_file
pig_tissue<-substr(pig_file,nchar("/Users/baizhonghao/Downloads/human-pig GTEx/SampleID_each_tissue/SampleID_")+1,nchar(pig_file)-nchar(".txt"))
pig_tissue

human_pig_expression[1:10,1:10]
human_sample<-rownames(human_pig_expression)[1:17382]
pig_sample<-rownames(human_pig_expression)[17383:22882]
experiment<-NULL
experiment$SampleID<-human_sample
experiment$Species<-"Human"
experiment$Tissues<-human_meta$SMTS[match(experiment$SampleID,human_meta$SAMPID)]
experiment$Tissue_type<-human_meta$SMTSD[match(experiment$SampleID,human_meta$SAMPID)]
experiment<-as.data.frame(experiment)

experiment_pig<-NULL
experiment_pig$SampleID<-pig_sample
experiment_pig$Species<-"Pig"
experiment_pig$Tissues<-pig_meta$`Categories-class2`[match(experiment_pig$SampleID,pig_meta$BioSample)]
experiment_pig$Tissue_type<-pig_meta$`Categories-class3`[match(experiment_pig$SampleID,pig_meta$BioSample)]
experiment_pig<-as.data.frame(experiment_pig)

experiment<-rbind(experiment,experiment_pig)
experiment<-na.omit(experiment)
######################################################################
#redraw the tissues#
experiment$Tissue_new<-experiment$Tissues
experiment$Tissue_new[which(experiment$Species=='Human'&experiment$Tissue_type=='Brain - Hypothalamus')]<-'Hypothalamus'
experiment$Tissue_new[which(experiment$Species=='Pig'&experiment$Tissue_type=='Colon')]<-'Colon'
experiment$Tissue_new[which(experiment$Species=='Human'&experiment$Tissue_type=='Brain - Frontal Cortex (BA9)')]<-'Frontal_cortex'
experiment$Tissue_new[which(experiment$Species=='Pig'&experiment$Tissue_type=='Frontal_cortex')]<-'Frontal_cortex'
experiment$Tissue_new[which(experiment$Species=='Pig'&experiment$Tissue_type=='Hypothalamus')]<-'Hypothalamus'
experiment$Tissue_new[which(experiment$Species=='Pig'&experiment$Tissue_type=='Pituitary')]<-'Pituitary'
experiment$Tissue_new[which(experiment$Species=='Pig'&experiment$Tissue_type=='Blastocyst')]<-'Blastocyst'
experiment$Tissue_new[which(experiment$Species=='Pig'&experiment$Tissue_type=='Blastomere')]<-'Blastomere'
experiment$Tissue_new[which(experiment$Species=='Pig'&experiment$Tissue_type=='Macrophage')]<-'Macrophage'
experiment$Tissue_new[which(experiment$Species=='Pig'&experiment$Tissue_type=='Duodenum')]<-'Duodenum'
experiment$Tissue_new[which(experiment$Species=='Pig'&experiment$Tissue_type=='Morula')]<-'Morula'
experiment$Tissue_new[which(experiment$Species=='Human'&experiment$Tissue_type=='Small Intestine - Terminal Ileum')]<-'Ileum'
experiment$Tissue_new[which(experiment$Species=='Pig'&experiment$Tissue_type=='Ileum')]<-'Ileum'
experiment$Tissue_new[which(experiment$Species=='Pig'&experiment$Tissue_type=='Jejunum')]<-'Jejunum'
experiment$Tissue_new[which(experiment$Species=='Human'&experiment$Tissue_type=='Cells - EBV-transformed lymphocytes')]<-'Blood'
experiment$Tissue_new[which(experiment$Species=='Human'&experiment$Tissue_type=='Breast - Mammary Tissue')]<-'Breast'
experiment<-experiment[,-5]
save(experiment,expression,file="/Users/baizhonghao/Downloads/human-pig GTEx/DataMatrix22761.RData")
write.xlsx(experiment,file="/Users/baizhonghao/Downloads/human-pig GTEx/pigGTEx metadata.xlsx")
commontissues<-which(experiment)

expression<-human_pig_expression[match(experiment$SampleID,rownames(human_pig_expression)),]
save(experiment,expression,file="/Users/baizhonghao/Downloads/human-pig GTEx/DataMatrix20217.RData")
