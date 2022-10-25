library(data.table)
library(ggplot2)
library(stringr)
#load data#
load("/Users/baizhonghao/Downloads/human-pig GTEx/DataMatrix20217.RData")
#Adipose#
eqtl_Adipose<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/sig_cis_qtl_pairs_pcg//Adipose.cis_qtl_pairs.significant.txt"))
eqtl_Adipose[1:10,1:9]
length(unique(eqtl_Adipose$phenotype_id)) #4078

eqtl_Adipose_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Adipose_Subcutaneous.v8.egenes.txt"))
eqtl_Adipose_hum2<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Adipose_Visceral_Omentum.v8.egenes.txt"))
eqtl_Adipose_hum<-rbind(eqtl_Adipose_hum1,eqtl_Adipose_hum2)
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))

#overlap egene#
pig_egene_Adipose<-unique(eqtl_Adipose$phenotype_id)
hum_egene_Adipose<-unique(substr(eqtl_Adipose_hum$gene_id[which(eqtl_Adipose_hum$qval<=0.05)],1,15))
pig2hum_Adipose<-annotation$`Gene stable ID`[match(pig_egene_Adipose,annotation$`Pig gene stable ID`)]
overlap_egene_Adipose<-intersect(hum_egene_Adipose,pig2hum_Adipose)#3474
overlap_one2oneE_Adipose<-intersect(overlap_egene_Adipose,one2one)
#non-egene#
one2one<-colnames(expression)#all one2one orthologous genes
non_humegene_Adipose<-setdiff(one2one,intersect(hum_egene_Adipose,one2one))#non-egene in hum
non_pigegene_Adipose<-setdiff(one2one,intersect(pig2hum_Adipose,one2one))#non-egene in pig
neither_egene_Adipose<-intersect(non_humegene_Adipose,non_pigegene_Adipose)#non-egene in hum and pig

#species 1-1 egenes#
one2oneE_hum_Adipose<-intersect(hum_egene_Adipose,one2one)#1-1 orthologous eGene in hum
one2oneE_pig_Adipose<-intersect(pig2hum_Adipose,one2one)#1-1 orthologous eGene in pig

#Artery#
eqtl_Artery<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/sig_cis_qtl_pairs_pcg/Artery.cis_qtl_pairs.significant.txt"))
eqtl_Artery_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Artery_Aorta.v8.egenes.txt"))
eqtl_Artery_hum2<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Artery_Coronary.v8.egenes.txt"))
eqtl_Artery_hum3<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Artery_Tibial.v8.egenes.txt"))
eqtl_Artery_hum<-rbind(eqtl_Artery_hum1,eqtl_Artery_hum2,eqtl_Artery_hum3)
one2one<-colnames(expression)#all one2one orthologous genes

#overlap egene#
pig_egene_Artery<-unique(eqtl_Artery$phenotype_id)
hum_egene_Artery<-unique(substr(eqtl_Artery_hum$gene_id[which(eqtl_Artery_hum$qval<=0.05)],1,15))
pig2hum_Artery<-annotation$`Gene stable ID`[match(pig_egene_Artery,annotation$`Pig gene stable ID`)]
overlap_egene_Artery<-intersect(hum_egene_Artery,pig2hum_Artery)
overlap_one2oneE_Artery<-intersect(overlap_egene_Artery,one2one)
#non-egene#

non_humegene_Artery<-setdiff(one2one,intersect(hum_egene_Artery,one2one))#non-egene in hum
non_pigegene_Artery<-setdiff(one2one,intersect(pig2hum_Artery,one2one))#non-egene in pig
neither_egene_Artery<-intersect(non_humegene_Artery,non_pigegene_Artery)#non-egene in hum and pig

#species 1-1 egenes#
one2oneE_hum_Artery<-intersect(hum_egene_Artery,one2one)#1-1 orthologous eGene in hum
one2oneE_pig_Artery<-intersect(pig2hum_Artery,one2one)#1-1 orthologous eGene in pig

#Blood#
eqtl_Blood<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/sig_cis_qtl_pairs_pcg/Blood.cis_qtl_pairs.significant.txt"))
eqtl_Blood_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt"))
eqtl_Blood_hum2<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Cells_EBV-transformed_lymphocytes.v8.egenes.txt"))
eqtl_Blood_hum<-rbind(eqtl_Blood_hum1,eqtl_Blood_hum2)

#overlap egene#
pig_egene_Blood<-unique(eqtl_Blood$phenotype_id)
hum_egene_Blood<-unique(substr(eqtl_Blood_hum$gene_id[which(eqtl_Blood_hum$qval<=0.05)],1,15))
pig2hum_Blood<-annotation$`Gene stable ID`[match(pig_egene_Blood,annotation$`Pig gene stable ID`)]
overlap_egene_Blood<-intersect(hum_egene_Blood,pig2hum_Blood)#3474
overlap_one2oneE_Blood<-intersect(overlap_egene_Blood,one2one)
#non-egene#

non_humegene_Blood<-setdiff(one2one,intersect(hum_egene_Blood,one2one))#non-egene in hum
non_pigegene_Blood<-setdiff(one2one,intersect(pig2hum_Blood,one2one))#non-egene in pig
neither_egene_Blood<-intersect(non_humegene_Blood,non_pigegene_Blood)#non-egene in hum and pig

#species 1-1 egenes#
one2oneE_hum_Blood<-intersect(hum_egene_Blood,one2one)#1-1 orthologous eGene in hum
one2oneE_pig_Blood<-intersect(pig2hum_Blood,one2one)#1-1 orthologous eGene in pig

#Colon#
eqtl_Colon<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/sig_cis_qtl_pairs_pcg/Colon.cis_qtl_pairs.significant.txt"))
eqtl_Colon_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Colon_Sigmoid.v8.egenes.txt"))
eqtl_Colon_hum2<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Colon_Transverse.v8.egenes.txt"))
eqtl_Colon_hum<-rbind(eqtl_Colon_hum1,eqtl_Colon_hum2)

#overlap egene#
pig_egene_Colon<-unique(eqtl_Colon$phenotype_id)
hum_egene_Colon<-unique(substr(eqtl_Colon_hum$gene_id[which(eqtl_Colon_hum$qval<=0.05)],1,15))
pig2hum_Colon<-annotation$`Gene stable ID`[match(pig_egene_Colon,annotation$`Pig gene stable ID`)]
overlap_egene_Colon<-intersect(hum_egene_Colon,pig2hum_Colon)#3474
overlap_one2oneE_Colon<-intersect(overlap_egene_Colon,one2one)
#non-egene#

non_humegene_Colon<-setdiff(one2one,intersect(hum_egene_Colon,one2one))#non-egene in hum
non_pigegene_Colon<-setdiff(one2one,intersect(pig2hum_Colon,one2one))#non-egene in pig
neither_egene_Colon<-intersect(non_humegene_Colon,non_pigegene_Colon)#non-egene in hum and pig

#species 1-1 egenes#
one2oneE_hum_Colon<-intersect(hum_egene_Colon,one2one)#1-1 orthologous eGene in hum
one2oneE_pig_Colon<-intersect(pig2hum_Colon,one2one)#1-1 orthologous eGene in pig

#Heart#
eqtl_Heart<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/sig_cis_qtl_pairs_pcg/Heart.cis_qtl_pairs.significant.txt"))
eqtl_Heart_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Heart_Atrial_Appendage.v8.egenes.txt"))
eqtl_Heart_hum2<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Heart_Left_Ventricle.v8.egenes.txt"))
eqtl_Heart_hum<-rbind(eqtl_Heart_hum1,eqtl_Heart_hum2)

#overlap egene#
pig_egene_Heart<-unique(eqtl_Heart$phenotype_id)
hum_egene_Heart<-unique(substr(eqtl_Heart_hum$gene_id[which(eqtl_Heart_hum$qval<=0.05)],1,15))
pig2hum_Heart<-annotation$`Gene stable ID`[match(pig_egene_Heart,annotation$`Pig gene stable ID`)]
overlap_egene_Heart<-intersect(hum_egene_Heart,pig2hum_Heart)#3474
overlap_one2oneE_Heart<-intersect(overlap_egene_Heart,one2one)
#non-egene#

non_humegene_Heart<-setdiff(one2one,intersect(hum_egene_Heart,one2one))#non-egene in hum
non_pigegene_Heart<-setdiff(one2one,intersect(pig2hum_Heart,one2one))#non-egene in pig
neither_egene_Heart<-intersect(non_humegene_Heart,non_pigegene_Heart)#non-egene in hum and pig

#species 1-1 egenes#
one2oneE_hum_Heart<-intersect(hum_egene_Heart,one2one)#1-1 orthologous eGene in hum
one2oneE_pig_Heart<-intersect(pig2hum_Heart,one2one)#1-1 orthologous eGene in pig

#Hypothalamus#
eqtl_Hypothalamus<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/sig_cis_qtl_pairs_pcg/Hypothalamus.cis_qtl_pairs.significant.txt"))
eqtl_Hypothalamus_hum<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Brain_Hypothalamus.v8.egenes.txt"))

#overlap egene#
pig_egene_Hypothalamus<-unique(eqtl_Hypothalamus$phenotype_id)
hum_egene_Hypothalamus<-unique(substr(eqtl_Hypothalamus_hum$gene_id[which(eqtl_Hypothalamus_hum$qval<=0.05)],1,15))
pig2hum_Hypothalamus<-annotation$`Gene stable ID`[match(pig_egene_Hypothalamus,annotation$`Pig gene stable ID`)]
overlap_egene_Hypothalamus<-intersect(hum_egene_Hypothalamus,pig2hum_Hypothalamus)#3474
overlap_one2oneE_Hypothalamus<-intersect(overlap_egene_Hypothalamus,one2one)
#non-egene#

non_humegene_Hypothalamus<-setdiff(one2one,intersect(hum_egene_Hypothalamus,one2one))#non-egene in hum
non_pigegene_Hypothalamus<-setdiff(one2one,intersect(pig2hum_Hypothalamus,one2one))#non-egene in pig
neither_egene_Hypothalamus<-intersect(non_humegene_Hypothalamus,non_pigegene_Hypothalamus)#non-egene in hum and pig

#species 1-1 egenes#
one2oneE_hum_Hypothalamus<-intersect(hum_egene_Hypothalamus,one2one)#1-1 orthologous eGene in hum
one2oneE_pig_Hypothalamus<-intersect(pig2hum_Hypothalamus,one2one)#1-1 orthologous eGene in pig

#Frontal cortex#
eqtl_Frontal_cortex<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/sig_cis_qtl_pairs_pcg/Frontal_cortex.cis_qtl_pairs.significant.txt"))
eqtl_Frontal_cortex_hum<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Brain_Frontal_Cortex_BA9.v8.egenes.txt"))

#overlap egene#
pig_egene_Frontal_cortex<-unique(eqtl_Frontal_cortex$phenotype_id)
hum_egene_Frontal_cortex<-unique(substr(eqtl_Frontal_cortex_hum$gene_id[which(eqtl_Frontal_cortex_hum$qval<=0.05)],1,15))
pig2hum_Frontal_cortex<-annotation$`Gene stable ID`[match(pig_egene_Frontal_cortex,annotation$`Pig gene stable ID`)]
overlap_egene_Frontal_cortex<-intersect(hum_egene_Frontal_cortex,pig2hum_Frontal_cortex)#3474
overlap_one2oneE_Frontal_cortex<-intersect(overlap_egene_Frontal_cortex,one2one)
#non-egene#

non_humegene_Frontal_cortex<-setdiff(one2one,intersect(hum_egene_Frontal_cortex,one2one))#non-egene in hum
non_pigegene_Frontal_cortex<-setdiff(one2one,intersect(pig2hum_Frontal_cortex,one2one))#non-egene in pig
neither_egene_Frontal_cortex<-intersect(non_humegene_Frontal_cortex,non_pigegene_Frontal_cortex)#non-egene in hum and pig

#species 1-1 egenes#
one2oneE_hum_Frontal_cortex<-intersect(hum_egene_Frontal_cortex,one2one)#1-1 orthologous eGene in hum
one2oneE_pig_Frontal_cortex<-intersect(pig2hum_Frontal_cortex,one2one)#1-1 orthologous eGene in pig

#Ileum#
eqtl_Ileum<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/sig_cis_qtl_pairs_pcg/Ileum.cis_qtl_pairs.significant.txt"))
eqtl_Ileum_hum<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Small_Intestine_Terminal_Ileum.v8.egenes.txt"))

#overlap egene#
pig_egene_Ileum<-unique(eqtl_Ileum$phenotype_id)
hum_egene_Ileum<-unique(substr(eqtl_Ileum_hum$gene_id[which(eqtl_Ileum_hum$qval<=0.05)],1,15))
pig2hum_Ileum<-annotation$`Gene stable ID`[match(pig_egene_Ileum,annotation$`Pig gene stable ID`)]
overlap_egene_Ileum<-intersect(hum_egene_Ileum,pig2hum_Ileum)#3474
overlap_one2oneE_Ileum<-intersect(overlap_egene_Ileum,one2one)
#non-egene#

non_humegene_Ileum<-setdiff(one2one,intersect(hum_egene_Ileum,one2one))#non-egene in hum
non_pigegene_Ileum<-setdiff(one2one,intersect(pig2hum_Ileum,one2one))#non-egene in pig
neither_egene_Ileum<-intersect(non_humegene_Ileum,non_pigegene_Ileum)#non-egene in hum and pig

#species 1-1 egenes#
one2oneE_hum_Ileum<-intersect(hum_egene_Ileum,one2one)#1-1 orthologous eGene in hum
one2oneE_pig_Ileum<-intersect(pig2hum_Ileum,one2one)#1-1 orthologous eGene in pig

#Kidney#
eqtl_Kidney<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/sig_cis_qtl_pairs_pcg/Kidney.cis_qtl_pairs.significant.txt"))
eqtl_Kidney_hum<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Kidney_Cortex.v8.egenes.txt"))

#overlap egene#
pig_egene_Kidney<-unique(eqtl_Kidney$phenotype_id)
hum_egene_Kidney<-unique(substr(eqtl_Kidney_hum$gene_id[which(eqtl_Kidney_hum$qval<=0.05)],1,15))
pig2hum_Kidney<-annotation$`Gene stable ID`[match(pig_egene_Kidney,annotation$`Pig gene stable ID`)]
overlap_egene_Kidney<-intersect(hum_egene_Kidney,pig2hum_Kidney)#3474
overlap_one2oneE_Kidney<-intersect(overlap_egene_Kidney,one2one)
#non-egene#

non_humegene_Kidney<-setdiff(one2one,intersect(hum_egene_Kidney,one2one))#non-egene in hum
non_pigegene_Kidney<-setdiff(one2one,intersect(pig2hum_Kidney,one2one))#non-egene in pig
neither_egene_Kidney<-intersect(non_humegene_Kidney,non_pigegene_Kidney)#non-egene in hum and pig

#species 1-1 egenes#
one2oneE_hum_Kidney<-intersect(hum_egene_Kidney,one2one)#1-1 orthologous eGene in hum
one2oneE_pig_Kidney<-intersect(pig2hum_Kidney,one2one)#1-1 orthologous eGene in pig

#Liver#
eqtl_Liver<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/sig_cis_qtl_pairs_pcg/Liver.cis_qtl_pairs.significant.txt"))
eqtl_Liver_hum<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Liver.v8.egenes.txt"))

#overlap egene#
pig_egene_Liver<-unique(eqtl_Liver$phenotype_id)
hum_egene_Liver<-unique(substr(eqtl_Liver_hum$gene_id[which(eqtl_Liver_hum$qval<=0.05)],1,15))
pig2hum_Liver<-annotation$`Gene stable ID`[match(pig_egene_Liver,annotation$`Pig gene stable ID`)]
overlap_egene_Liver<-intersect(hum_egene_Liver,pig2hum_Liver)#3474
overlap_one2oneE_Liver<-intersect(overlap_egene_Liver,one2one)
#non-egene#

non_humegene_Liver<-setdiff(one2one,intersect(hum_egene_Liver,one2one))#non-egene in hum
non_pigegene_Liver<-setdiff(one2one,intersect(pig2hum_Liver,one2one))#non-egene in pig
neither_egene_Liver<-intersect(non_humegene_Liver,non_pigegene_Liver)#non-egene in hum and pig

#species 1-1 egenes#
one2oneE_hum_Liver<-intersect(hum_egene_Liver,one2one)#1-1 orthologous eGene in hum
one2oneE_pig_Liver<-intersect(pig2hum_Liver,one2one)#1-1 orthologous eGene in pig

#Lung#
eqtl_Lung<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/sig_cis_qtl_pairs_pcg/Lung.cis_qtl_pairs.significant.txt"))
eqtl_Lung_hum<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Lung.v8.egenes.txt"))

#overlap egene#
pig_egene_Lung<-unique(eqtl_Lung$phenotype_id)
hum_egene_Lung<-unique(substr(eqtl_Lung_hum$gene_id[which(eqtl_Lung_hum$qval<=0.05)],1,15))
pig2hum_Lung<-annotation$`Gene stable ID`[match(pig_egene_Lung,annotation$`Pig gene stable ID`)]
overlap_egene_Lung<-intersect(hum_egene_Lung,pig2hum_Lung)#3474
overlap_one2oneE_Lung<-intersect(overlap_egene_Lung,one2one)
#non-egene#

non_humegene_Lung<-setdiff(one2one,intersect(hum_egene_Lung,one2one))#non-egene in hum
non_pigegene_Lung<-setdiff(one2one,intersect(pig2hum_Lung,one2one))#non-egene in pig
neither_egene_Lung<-intersect(non_humegene_Lung,non_pigegene_Lung)#non-egene in hum and pig

#species 1-1 egenes#
one2oneE_hum_Lung<-intersect(hum_egene_Lung,one2one)#1-1 orthologous eGene in hum
one2oneE_pig_Lung<-intersect(pig2hum_Lung,one2one)#1-1 orthologous eGene in pig

#Muscle#
eqtl_Muscle<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/sig_cis_qtl_pairs_pcg/Muscle.cis_qtl_pairs.significant.txt"))
eqtl_Muscle_hum<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Muscle_Skeletal.v8.egenes.txt"))

#overlap egene#
pig_egene_Muscle<-unique(eqtl_Muscle$phenotype_id)
hum_egene_Muscle<-unique(substr(eqtl_Muscle_hum$gene_id[which(eqtl_Muscle_hum$qval<=0.05)],1,15))
pig2hum_Muscle<-annotation$`Gene stable ID`[match(pig_egene_Muscle,annotation$`Pig gene stable ID`)]
overlap_egene_Muscle<-intersect(hum_egene_Muscle,pig2hum_Muscle)#3474
overlap_one2oneE_Muscle<-intersect(overlap_egene_Muscle,one2one)
#non-egene#

non_humegene_Muscle<-setdiff(one2one,intersect(hum_egene_Muscle,one2one))#non-egene in hum
non_pigegene_Muscle<-setdiff(one2one,intersect(pig2hum_Muscle,one2one))#non-egene in pig
neither_egene_Muscle<-intersect(non_humegene_Muscle,non_pigegene_Muscle)#non-egene in hum and pig

#species 1-1 egenes#
one2oneE_hum_Muscle<-intersect(hum_egene_Muscle,one2one)#1-1 orthologous eGene in hum
one2oneE_pig_Muscle<-intersect(pig2hum_Muscle,one2one)#1-1 orthologous eGene in pig

#Ovary#
eqtl_Ovary<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/sig_cis_qtl_pairs_pcg/Ovary.cis_qtl_pairs.significant.txt"))
eqtl_Ovary_hum<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Ovary.v8.egenes.txt"))

#overlap egene#
pig_egene_Ovary<-unique(eqtl_Ovary$phenotype_id)
hum_egene_Ovary<-unique(substr(eqtl_Ovary_hum$gene_id[which(eqtl_Ovary_hum$qval<=0.05)],1,15))
pig2hum_Ovary<-annotation$`Gene stable ID`[match(pig_egene_Ovary,annotation$`Pig gene stable ID`)]
overlap_egene_Ovary<-intersect(hum_egene_Ovary,pig2hum_Ovary)#3474
overlap_one2oneE_Ovary<-intersect(overlap_egene_Ovary,one2one)
#non-egene#

non_humegene_Ovary<-setdiff(one2one,intersect(hum_egene_Ovary,one2one))#non-egene in hum
non_pigegene_Ovary<-setdiff(one2one,intersect(pig2hum_Ovary,one2one))#non-egene in pig
neither_egene_Ovary<-intersect(non_humegene_Ovary,non_pigegene_Ovary)#non-egene in hum and pig

#species 1-1 egenes#
one2oneE_hum_Ovary<-intersect(hum_egene_Ovary,one2one)#1-1 orthologous eGene in hum
one2oneE_pig_Ovary<-intersect(pig2hum_Ovary,one2one)#1-1 orthologous eGene in pig

#Pituitary#
eqtl_Pituitary<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/sig_cis_qtl_pairs_pcg/Pituitary.cis_qtl_pairs.significant.txt"))
eqtl_Pituitary_hum<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Pituitary.v8.egenes.txt"))

#overlap egene#
pig_egene_Pituitary<-unique(eqtl_Pituitary$phenotype_id)
hum_egene_Pituitary<-unique(substr(eqtl_Pituitary_hum$gene_id[which(eqtl_Pituitary_hum$qval<=0.05)],1,15))
pig2hum_Pituitary<-annotation$`Gene stable ID`[match(pig_egene_Pituitary,annotation$`Pig gene stable ID`)]
overlap_egene_Pituitary<-intersect(hum_egene_Pituitary,pig2hum_Pituitary)#3474
overlap_one2oneE_Pituitary<-intersect(overlap_egene_Pituitary,one2one)
#non-egene#

non_humegene_Pituitary<-setdiff(one2one,intersect(hum_egene_Pituitary,one2one))#non-egene in hum
non_pigegene_Pituitary<-setdiff(one2one,intersect(pig2hum_Pituitary,one2one))#non-egene in pig
neither_egene_Pituitary<-intersect(non_humegene_Pituitary,non_pigegene_Pituitary)#non-egene in hum and pig

#species 1-1 egenes#
one2oneE_hum_Pituitary<-intersect(hum_egene_Pituitary,one2one)#1-1 orthologous eGene in hum
one2oneE_pig_Pituitary<-intersect(pig2hum_Pituitary,one2one)#1-1 orthologous eGene in pig

#Spleen#
eqtl_Spleen<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/sig_cis_qtl_pairs_pcg/Spleen.cis_qtl_pairs.significant.txt"))
eqtl_Spleen_hum<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Spleen.v8.egenes.txt"))

#overlap egene#
pig_egene_Spleen<-unique(eqtl_Spleen$phenotype_id)
hum_egene_Spleen<-unique(substr(eqtl_Spleen_hum$gene_id[which(eqtl_Spleen_hum$qval<=0.05)],1,15))
pig2hum_Spleen<-annotation$`Gene stable ID`[match(pig_egene_Spleen,annotation$`Pig gene stable ID`)]
overlap_egene_Spleen<-intersect(hum_egene_Spleen,pig2hum_Spleen)#3474
overlap_one2oneE_Spleen<-intersect(overlap_egene_Spleen,one2one)
#non-egene#

non_humegene_Spleen<-setdiff(one2one,intersect(hum_egene_Spleen,one2one))#non-egene in hum
non_pigegene_Spleen<-setdiff(one2one,intersect(pig2hum_Spleen,one2one))#non-egene in pig
neither_egene_Spleen<-intersect(non_humegene_Spleen,non_pigegene_Spleen)#non-egene in hum and pig

#species 1-1 egenes#
one2oneE_hum_Spleen<-intersect(hum_egene_Spleen,one2one)#1-1 orthologous eGene in hum
one2oneE_pig_Spleen<-intersect(pig2hum_Spleen,one2one)#1-1 orthologous eGene in pig

#Testis#
eqtl_Testis<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/sig_cis_qtl_pairs_pcg/Testis.cis_qtl_pairs.significant.txt"))
eqtl_Testis_hum<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Testis.v8.egenes.txt"))

#overlap egene#
pig_egene_Testis<-unique(eqtl_Testis$phenotype_id)
hum_egene_Testis<-unique(substr(eqtl_Testis_hum$gene_id[which(eqtl_Testis_hum$qval<=0.05)],1,15))
pig2hum_Testis<-annotation$`Gene stable ID`[match(pig_egene_Testis,annotation$`Pig gene stable ID`)]
overlap_egene_Testis<-intersect(hum_egene_Testis,pig2hum_Testis)#3474
overlap_one2oneE_Testis<-intersect(overlap_egene_Testis,one2one)
#non-egene#

non_humegene_Testis<-setdiff(one2one,intersect(hum_egene_Testis,one2one))#non-egene in hum
non_pigegene_Testis<-setdiff(one2one,intersect(pig2hum_Testis,one2one))#non-egene in pig
neither_egene_Testis<-intersect(non_humegene_Testis,non_pigegene_Testis)#non-egene in hum and pig

#species 1-1 egenes#
one2oneE_hum_Testis<-intersect(hum_egene_Testis,one2one)#1-1 orthologous eGene in hum
one2oneE_pig_Testis<-intersect(pig2hum_Testis,one2one)#1-1 orthologous eGene in pig

#Uterus#
eqtl_Uterus<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/sig_cis_qtl_pairs_pcg/Uterus.cis_qtl_pairs.significant.txt"))
eqtl_Uterus_hum<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Uterus.v8.egenes.txt"))

#overlap egene#
pig_egene_Uterus<-unique(eqtl_Uterus$phenotype_id)
hum_egene_Uterus<-unique(substr(eqtl_Uterus_hum$gene_id[which(eqtl_Uterus_hum$qval<=0.05)],1,15))
pig2hum_Uterus<-annotation$`Gene stable ID`[match(pig_egene_Uterus,annotation$`Pig gene stable ID`)]
overlap_egene_Uterus<-intersect(hum_egene_Uterus,pig2hum_Uterus)#3474
overlap_one2oneE_Uterus<-intersect(overlap_egene_Uterus,one2one)
#non-egene#

non_humegene_Uterus<-setdiff(one2one,intersect(hum_egene_Uterus,one2one))#non-egene in hum
non_pigegene_Uterus<-setdiff(one2one,intersect(pig2hum_Uterus,one2one))#non-egene in pig
neither_egene_Uterus<-intersect(non_humegene_Uterus,non_pigegene_Uterus)#non-egene in hum and pig

#species 1-1 egenes#
one2oneE_hum_Uterus<-intersect(hum_egene_Uterus,one2one)#1-1 orthologous eGene in hum
one2oneE_pig_Uterus<-intersect(pig2hum_Uterus,one2one)#1-1 orthologous eGene in pig

################################################################################
#tissue info#
Human_experiment<-experiment[which(experiment$Species=="Human"),]
tissue_name<-unique(sort(Human_experiment$Tissue_new))
tissue_name
Pig_experiment<-experiment[which(experiment$Species=="Pig"),]
pig_tissue<-unique(sort(Pig_experiment$Tissue_new))
pig_tissue
tissues<-intersect(tissue_name,pig_tissue)
tissues

#summarize 4 gene groups#
all_non_egene<-get(paste0("neither_egene_",tissues[i]))
for(i in 2:length(tissues)){
  all_non_egene<-intersect(all_non_egene,get(paste0("neither_egene_",tissues[i])))
}       #non-egene in all common tissues


summary_hum<-array(NA, dim=c(length(colnames(expression)),length(tissues)))
summary_pig<-array(NA, dim=c(length(colnames(expression)),length(tissues)))
rownames(summary_hum)<-colnames(expression)
colnames(summary_hum)<-tissues
rownames(summary_pig)<-colnames(expression)
colnames(summary_pig)<-tissues

for(i in 1:length(rownames(summary_hum))){
  for(j in 1:length(tissues)){
    if(rownames(summary_hum)[i]%in%get(paste0("one2oneE_hum_",tissues[j]))==TRUE){
      summary_hum[i,j]<-1
    }else{
      summary_hum[i,j]<-0
    }
  }
}

for(i in 1:length(rownames(summary_pig))){
  for(j in 1:length(tissues)){
    if(rownames(summary_pig)[i]%in%get(paste0("one2oneE_pig_",tissues[j]))==TRUE){
      summary_pig[i,j]<-1
    }else{
      summary_pig[i,j]<-0
    }
  }
}

summary_hum<-as.data.frame(summary_hum)
summary_pig<-as.data.frame(summary_pig)
summary_hum$number<-apply(summary_hum,1,sum)
summary_pig$number<-apply(summary_pig,1,sum)
Results<-array(NA,dim=c(length(rownames(summary_hum)),2))
colnames(Results)<-c("Human","Pig")
rownames(Results)<-rownames(summary_hum)
Results<-as.data.frame(Results)
Results$Human<-summary_hum$number
Results$Pig<-summary_pig$number

cor<-cor(Results$Human,Results$Pig,method="spearman")
pval<-t.test(Results$Human,Results$Pig,method="spearman")

#calculate the genes were identified as eGenes in how many tissues#

Result<-NULL
for(i in 0:17){
  for(j in 0:17){
    summary<-NULL
    a<-nrow(Results[Results$Human==i&Results$Pig==j,])
    summary$Human<-i
    summary$Pig<-j
    summary$number<-a
    summary<-as.data.frame(summary)
    Result<-rbind(Result,summary)
  }
}


Result<-as.data.frame(Result)
Result$number[which(Result$number==0)]<-NA
nrow(Results[Results$Human==17&Results$Pig==0,])

dotplot<-ggplot(Result,aes(x=Human,y=Pig))+
  theme_classic()+
  geom_point(aes(size=number))+
  labs(x='Number of tissues eGenes \nare actived in Human',y='Number of tissues eGenes \nare actived in Pig')+
  theme(axis.title = element_text(color = "black", size = unit(12.4, "pt")),axis.text = element_text(color = "black", size = unit(12.4, "pt")))+
  theme(legend.title = element_text(size = unit(7, "pt")), legend.text = element_text(size = unit(6, "pt")))
ggsave(dotplot,file='/Users/baizhonghao/Downloads/eGenetissuenumber.pdf',dpi = 300,width = 80,height=80,units = 'mm')
###############################################################################
#calculate OR value#
ORsum<-array(NA, dim=c(length(tissues),length(tissues)))
colnames(ORsum)<-tissues
rownames(ORsum)<-tissues
ORsum<-as.data.frame(ORsum)
Vennsum<-NULL

for(i in 1:length(tissues)){
  for(j in 1:length(tissues)){
    a<-array(NA,dim=c(1,8))
    colnames(a)<-c("human_tissue","pig_tissue","both","human-only","pig-only","neither","OR","p-value")
    a<-as.data.frame(a)
    a$human_tissue<-tissues[i]
    a$pig_tissue<-tissues[j]
    overlap<-intersect(get(paste0("hum_egene_",tissues[i])),get(paste0("pig2hum_",tissues[j])))
    a$both<-length(intersect(overlap,one2one))#intersect the same egenes between two tissues
    a$`human-only`<-length(get(paste0("one2oneE_hum_",tissues[i]))) - a$both
    a$`pig-only`<-length(get(paste0("one2oneE_pig_",tissues[j]))) - a$both
    a$neither<-length(intersect(get(paste0("non_humegene_",tissues[i])),get(paste0("non_pigegene_",tissues[j]))))
    dat <- matrix(unlist(c(a[3:6])), nrow = 2, ncol = 2)
    test<-fisher.test(dat)
    a$OR<-as.numeric(test$estimate)
    a$`p-value`<-test$p.value + 1e-201
    Venn<-list(group1=get(paste0("one2oneE_hum_",tissues[i])),group2=get(paste0("one2oneE_pig_",tissues[j])))
    venn.diagram(Venn,filename =paste0("/Users/baizhonghao/Downloads/human-pig GTEx/Venn/",tissues[i],"_",tissues[j],".png"),
                 imagetype = "png",
                 fill = c('red', 'blue'), alpha = 0.50, cat.col = rep('black', 2), 
                 col = 'black', cex = 1.5, fontfamily = 'serif', 
                 cat.cex = 1.5, cat.fontfamily = 'serif')
    Vennsum<-rbind(Vennsum,a)
  }
}

ORheatmap<-Vennsum[,c(1,2,7)]
ORheatmap<-dcast(ORheatmap,human_tissue~pig_tissue)
rownames(ORheatmap)<-ORheatmap$human_tissue
ORheatmap<-ORheatmap[,-1]


tissuegroup<-c("Artery","Muscle","Heart","Frontal_cortex","Hypothalamus","Pituitary","Blood","Spleen","Colon","Ileum","Liver","Lung","Kidney","Adipose","Ovary","Testis","Uterus")
ORheatmap_group<-ORheatmap[match(tissuegroup,rownames(ORheatmap)),]
ORheatmap_group<-ORheatmap_group[,match(tissuegroup,colnames(ORheatmap))]
library(pheatmap)
library(RColorBrewer)
cc = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))
bk<-c(seq(min(ORheatmap),1,length.out = 50),seq(1.01,max(ORheatmap),length.out = 50))
pheatmap(ORheatmap_group,file="/Users/baizhonghao/Downloads/hc-OR_clustered.pdf",
         cluster_rows = F,cluster_cols = F,
         fontsize_row = 13,fontsize_col = 13,
         cellwidth = 13,cellheight = 13,
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         breaks = bk)

library(stats)
Vennsum$`p.adj`<-p.adjust(Vennsum$`p-value`,method="BH")
Vennsum$`-log10 p.adj`<--log10(Vennsum$p.adj)
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/eqtl/OR-plot.tiff",##reqiured to change
     res = 300, width =2000, height = 2000,compression = "lzw")
ggplot(Vennsum,aes(x=human_tissue,y=pig_tissue))+
  geom_point(aes(size=OR,color=`-log10 p.adj`))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low="lightgrey",high="red")+
  labs(x='Human',y='Pig')+
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  theme(axis.title=element_text(face = "bold", color="black"))
dev.off()



OR_matrix<-melt(ORheatmap_group)
OR_matrix$Human<-rep(c("Artery","Muscle","Heart","Frontal_cortex","Hypothalamus","Pituitary","Blood","Spleen","Colon","Ileum","Liver","Lung","Kidney","Adipose","Ovary","Testis","Uterus"),17)
colnames(OR_matrix)[1]<-'Pig'
colnames(OR_matrix)[2]<-'OR'
OR_matrix$Category<-'unmatched tissue pairs'
OR_matrix$Category[which(OR_matrix$Pig==OR_matrix$Human)]<-'matched tissue pairs'
OR_matrix<-OR_matrix[order(OR_matrix$Category,decreasing=T),]
OR_matrix$Pig<-factor(OR_matrix$Pig,levels=c("Artery","Muscle","Heart","Frontal_cortex","Hypothalamus","Pituitary","Blood","Spleen","Colon","Ileum","Liver","Lung","Kidney","Adipose","Ovary","Testis","Uterus"))
OR_matrix$Category<-factor(OR_matrix$Category, levels=c('unmatched tissue pairs','matched tissue pairs'))

orviolin<-ggplot(OR_matrix,aes(x=Pig,y=OR))+
  geom_violin(fill='#94E1ED',colour='#94E1ED', trim = FALSE) +
  geom_jitter(aes(x = Pig, y =OR ,color=Category, size=Category), width=0.2)+
  scale_color_manual(values=c('unmatched tissue pairs'='grey','matched tissue pairs'='#C11111'))+
  scale_size_manual(values=c('unmatched tissue pairs'=2, 'matched tissue pairs'=3))+
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45,size = 15))+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.title.y=element_text(size=18, face="bold"))+
  theme(legend.text = element_text(size=15))+
  theme(legend.title = element_text(size=18,face="bold"))+
  theme(legend.position=c(1,1),legend.justification = c(1,1))+
  theme(legend.background = element_rect(color="black"))+
  scale_y_continuous(limits = c(0,4.2))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x='Tissues',y='Odds Ratio')
ggsave(orviolin,file='/Users/baizhonghao/Downloads/orviolin.pdf',dpi = 300,width=15,height = 6)
##################################################
#SNP#

Adipose_SNP<-strsplit(FM_Adipose$variant_id,split = "_")
Adipose_SNP<-as.data.frame(Adipose_SNP)
Adipose_SNP<-as.data.frame(t(Adipose_SNP))
write.table(Adipose_SNP[,c(1,2)],file="/Users/baizhonghao/Downloads/Adipose")
chr1<-read.csv("/Users/baizhonghao/Downloads/liftover/chr1.csv")
chr2<-read.csv("/Users/baizhonghao/Downloads/liftover/chr2.csv")
chr3<-read.csv("/Users/baizhonghao/Downloads/liftover/chr3.csv")
chr4<-read.csv("/Users/baizhonghao/Downloads/liftover/chr4.csv")
chr5<-read.csv("/Users/baizhonghao/Downloads/liftover/chr5.csv")
chr6<-read.csv("/Users/baizhonghao/Downloads/liftover/chr6.csv")
chr7<-read.csv("/Users/baizhonghao/Downloads/liftover/chr7.csv")
chr8<-read.csv("/Users/baizhonghao/Downloads/liftover/chr8.csv")
chr9<-read.csv("/Users/baizhonghao/Downloads/liftover/chr9.csv")
chr10<-read.csv("/Users/baizhonghao/Downloads/liftover/chr10.csv")
chr11<-read.csv("/Users/baizhonghao/Downloads/liftover/chr11.csv")
chr12<-read.csv("/Users/baizhonghao/Downloads/liftover/chr12.csv")
chr13<-read.csv("/Users/baizhonghao/Downloads/liftover/chr13.csv")
chr14<-read.csv("/Users/baizhonghao/Downloads/liftover/chr14.csv")
chr15<-read.csv("/Users/baizhonghao/Downloads/liftover/chr15.csv")
chr16<-read.csv("/Users/baizhonghao/Downloads/liftover/chr16.csv")
chr17<-read.csv("/Users/baizhonghao/Downloads/liftover/chr17.csv")
chr18<-read.csv("/Users/baizhonghao/Downloads/liftover/chr18.csv")

#Adipose#
FM_Adipose<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/fine_mapping_results/common tissues/Adipose.cis_independent_qtl.txt"))
FM_Adipose_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Adipose_Subcutaneous.v8.egenes.txt"))
FM_Adipose_hum2<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Adipose_Visceral_Omentum.v8.egenes.txt"))
FM_Adipose_hum<-rbind(FM_Adipose_hum1,FM_Adipose_hum2)
FM_Adipose_hum$gene_id<-substr(FM_Adipose_hum$gene_id,1,15)
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
FMgene_pig_Adipose<-annotation$`Gene stable ID`[match(unique(FM_Adipose$phenotype_id),annotation$`Pig gene stable ID`)]
FMgene_one2one_Adipose<-intersect(FMgene_pig_Adipose,colnames(expression))
overlap_one2oneF_Adipose<-intersect(FMgene_one2one_Adipose,unique(FM_Adipose_hum$gene_id))
overlap_one2oneF_Adipose_pig<-annotation$`Pig gene stable ID`[match(overlap_one2oneF_Adipose,annotation$`Gene stable ID`)]

SNP_pig_Adipose<-NULL
for(i in 1:length(overlap_one2oneF_Adipose)){
  a<-FM_Adipose[which(FM_Adipose$phenotype_id==overlap_one2oneF_Adipose_pig[i]),]
  a$human_id<-overlap_one2oneF_Adipose[i]
  SNP_pig_Adipose<-rbind(SNP_pig_Adipose,a)
}

SNP_hum_Adipose<-NULL
for(i in 1:length(overlap_one2oneF_Adipose)){
  a<-FM_Adipose_hum[which(FM_Adipose_hum$gene_id==overlap_one2oneF_Adipose[i]),]
  SNP_hum_Adipose<-rbind(SNP_hum_Adipose,a)
}

variant_Adipose<-strsplit(SNP_pig_Adipose$variant_id,split = "_")
variant_Adipose<-as.data.frame(variant_Adipose)
variant_Adipose<-as.data.frame(t(variant_Adipose))
SNP_pig_Adipose$chr <- variant_Adipose[,1]
SNP_pig_Adipose$SNP<-variant_Adipose[,2]
SNP_pig_Adipose$ref <- variant_Adipose[,3]
SNP_pig_Adipose$alt <- variant_Adipose[,4]

SNP_Adipose_pig<-NULL
SNP_Adipose_hum<-NULL
for(i in 1:length(overlap_one2oneF_Adipose)){
  SNP_hum<-NULL
  SNP_pig<-NULL
  pig_variant<-SNP_pig_Adipose[which(SNP_pig_Adipose$human_id==overlap_one2oneF_Adipose[i]),]
  hum_variant<-SNP_hum_Adipose[which(SNP_hum_Adipose$gene_id==overlap_one2oneF_Adipose[i]),]
  for(j in 1:length(nrow(pig_variant))){
    pair_pig<-paste0(pig_variant$ref[j],pig_variant$alt[j])
    pair_pig_trans<-paste0(pig_variant$alt[j],pig_variant$ref[j])
    for(k in 1:length(nrow(hum_variant))){
      pair_hum<-paste0(hum_variant$ref[k],hum_variant$alt[k])
      if(pair_pig==pair_hum|pair_hum==pair_pig_trans){
        SNP_hum<-rbind(SNP_hum,hum_variant[k,])
        SNP_pig<-rbind(SNP_pig,pig_variant[j,])
      }
    }
    SNP_Adipose_pig<-rbind(SNP_Adipose_pig,SNP_pig)
    SNP_Adipose_hum<-rbind(SNP_Adipose_hum,SNP_hum)
  }
}

#SNP_Adipose_hum<-SNP_Adipose_hum[!duplicated(SNP_Adipose_hum, fromLast=TRUE), ]

#Artery#
FM_Artery<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/fine_mapping_results/common tissues/Artery.cis_independent_qtl.txt"))
FM_Artery_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Artery_Aorta.v8.egenes.txt"))
FM_Artery_hum2<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Artery_Coronary.v8.egenes.txt"))
FM_Artery_hum3<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Artery_Tibial.v8.egenes.txt"))
FM_Artery_hum<-rbind(FM_Artery_hum1,FM_Artery_hum2,FM_Artery_hum3)
FM_Artery_hum$gene_id<-substr(FM_Artery_hum$gene_id,1,15)
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
FMgene_pig_Artery<-annotation$`Gene stable ID`[match(unique(FM_Artery$phenotype_id),annotation$`Pig gene stable ID`)]
FMgene_one2one_Artery<-intersect(FMgene_pig_Artery,colnames(expression))
overlap_one2oneF_Artery<-intersect(FMgene_one2one_Artery,unique(FM_Artery_hum$gene_id))
overlap_one2oneF_Artery_pig<-annotation$`Pig gene stable ID`[match(overlap_one2oneF_Artery,annotation$`Gene stable ID`)]

SNP_pig_Artery<-NULL
for(i in 1:length(overlap_one2oneF_Artery)){
  a<-FM_Artery[which(FM_Artery$phenotype_id==overlap_one2oneF_Artery_pig[i]),]
  a$human_id<-overlap_one2oneF_Artery[i]
  SNP_pig_Artery<-rbind(SNP_pig_Artery,a)
}

SNP_hum_Artery<-NULL
for(i in 1:length(overlap_one2oneF_Artery)){
  a<-FM_Artery_hum[which(FM_Artery_hum$gene_id==overlap_one2oneF_Artery[i]),]
  SNP_hum_Artery<-rbind(SNP_hum_Artery,a)
}

variant_Artery<-strsplit(SNP_pig_Artery$variant_id,split = "_")
variant_Artery<-as.data.frame(variant_Artery)
variant_Artery<-as.data.frame(t(variant_Artery))
SNP_pig_Artery$chr <- variant_Artery[,1]
SNP_pig_Artery$SNP<-variant_Artery[,2]
SNP_pig_Artery$ref <- variant_Artery[,3]
SNP_pig_Artery$alt <- variant_Artery[,4]

SNP_Artery_pig<-NULL
SNP_Artery_hum<-NULL
for(i in 1:length(overlap_one2oneF_Artery)){
  SNP_hum<-NULL
  SNP_pig<-NULL
  pig_variant<-SNP_pig_Artery[which(SNP_pig_Artery$human_id==overlap_one2oneF_Artery[i]),]
  hum_variant<-SNP_hum_Artery[which(SNP_hum_Artery$gene_id==overlap_one2oneF_Artery[i]),]
  for(j in 1:length(nrow(pig_variant))){
    pair_pig<-paste0(pig_variant$ref[j],pig_variant$alt[j])
    pair_pig_trans<-paste0(pig_variant$alt[j],pig_variant$ref[j])
    for(k in 1:length(nrow(hum_variant))){
      pair_hum<-paste0(hum_variant$ref[k],hum_variant$alt[k])
      if(pair_pig==pair_hum|pair_hum==pair_pig_trans){
        SNP_hum<-rbind(SNP_hum,hum_variant[k,])
        SNP_pig<-rbind(SNP_pig,pig_variant[j,])
      }
    }
    SNP_Artery_pig<-rbind(SNP_Artery_pig,SNP_pig)
    SNP_Artery_hum<-rbind(SNP_Artery_hum,SNP_hum)
  }
}

#Blood#
FM_Blood<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/fine_mapping_results/common tissues/Blood.cis_independent_qtl.txt"))
FM_Blood_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt"))
FM_Blood_hum2<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Cells_EBV-transformed_lymphocytes.v8.egenes.txt"))
FM_Blood_hum<-rbind(FM_Blood_hum1,FM_Blood_hum2)
FM_Blood_hum$gene_id<-substr(FM_Blood_hum$gene_id,1,15)
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
FMgene_pig_Blood<-annotation$`Gene stable ID`[match(unique(FM_Blood$phenotype_id),annotation$`Pig gene stable ID`)]
FMgene_one2one_Blood<-intersect(FMgene_pig_Blood,colnames(expression))
overlap_one2oneF_Blood<-intersect(FMgene_one2one_Blood,unique(FM_Blood_hum$gene_id))
overlap_one2oneF_Blood_pig<-annotation$`Pig gene stable ID`[match(overlap_one2oneF_Blood,annotation$`Gene stable ID`)]

SNP_pig_Blood<-NULL
for(i in 1:length(overlap_one2oneF_Blood)){
  a<-FM_Blood[which(FM_Blood$phenotype_id==overlap_one2oneF_Blood_pig[i]),]
  a$human_id<-overlap_one2oneF_Blood[i]
  SNP_pig_Blood<-rbind(SNP_pig_Blood,a)
}

SNP_hum_Blood<-NULL
for(i in 1:length(overlap_one2oneF_Blood)){
  a<-FM_Blood_hum[which(FM_Blood_hum$gene_id==overlap_one2oneF_Blood[i]),]
  SNP_hum_Blood<-rbind(SNP_hum_Blood,a)
}

variant_Blood<-strsplit(SNP_pig_Blood$variant_id,split = "_")
variant_Blood<-as.data.frame(variant_Blood)
variant_Blood<-as.data.frame(t(variant_Blood))
SNP_pig_Blood$chr <- variant_Blood[,1]
SNP_pig_Blood$SNP<-variant_Blood[,2]
SNP_pig_Blood$ref <- variant_Blood[,3]
SNP_pig_Blood$alt <- variant_Blood[,4]

SNP_Blood_pig<-NULL
SNP_Blood_hum<-NULL
for(i in 1:length(overlap_one2oneF_Blood)){
  SNP_hum<-NULL
  SNP_pig<-NULL
  pig_variant<-SNP_pig_Blood[which(SNP_pig_Blood$human_id==overlap_one2oneF_Blood[i]),]
  hum_variant<-SNP_hum_Blood[which(SNP_hum_Blood$gene_id==overlap_one2oneF_Blood[i]),]
  for(j in 1:length(nrow(pig_variant))){
    pair_pig<-paste0(pig_variant$ref[j],pig_variant$alt[j])
    pair_pig_trans<-paste0(pig_variant$alt[j],pig_variant$ref[j])
    for(k in 1:length(nrow(hum_variant))){
      pair_hum<-paste0(hum_variant$ref[k],hum_variant$alt[k])
      if(pair_pig==pair_hum|pair_hum==pair_pig_trans){
        SNP_hum<-rbind(SNP_hum,hum_variant[k,])
        SNP_pig<-rbind(SNP_pig,pig_variant[j,])
      }
    }
    SNP_Blood_pig<-rbind(SNP_Blood_pig,SNP_pig)
    SNP_Blood_hum<-rbind(SNP_Blood_hum,SNP_hum)
  }
}

#Colon#
FM_Colon<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/fine_mapping_results/common tissues/Colon.cis_independent_qtl.txt"))
FM_Colon_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Colon_Sigmoid.v8.egenes.txt"))
FM_Colon_hum2<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Colon_Transverse.v8.egenes.txt"))
FM_Colon_hum<-rbind(FM_Colon_hum1,FM_Colon_hum2)
FM_Colon_hum$gene_id<-substr(FM_Colon_hum$gene_id,1,15)
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
FMgene_pig_Colon<-annotation$`Gene stable ID`[match(unique(FM_Colon$phenotype_id),annotation$`Pig gene stable ID`)]
FMgene_one2one_Colon<-intersect(FMgene_pig_Colon,colnames(expression))
overlap_one2oneF_Colon<-intersect(FMgene_one2one_Colon,unique(FM_Colon_hum$gene_id))
overlap_one2oneF_Colon_pig<-annotation$`Pig gene stable ID`[match(overlap_one2oneF_Colon,annotation$`Gene stable ID`)]

SNP_pig_Colon<-NULL
for(i in 1:length(overlap_one2oneF_Colon)){
  a<-FM_Colon[which(FM_Colon$phenotype_id==overlap_one2oneF_Colon_pig[i]),]
  a$human_id<-overlap_one2oneF_Colon[i]
  SNP_pig_Colon<-rbind(SNP_pig_Colon,a)
}

SNP_hum_Colon<-NULL
for(i in 1:length(overlap_one2oneF_Colon)){
  a<-FM_Colon_hum[which(FM_Colon_hum$gene_id==overlap_one2oneF_Colon[i]),]
  SNP_hum_Colon<-rbind(SNP_hum_Colon,a)
}

variant_Colon<-strsplit(SNP_pig_Colon$variant_id,split = "_")
variant_Colon<-as.data.frame(variant_Colon)
variant_Colon<-as.data.frame(t(variant_Colon))
SNP_pig_Colon$chr <- variant_Colon[,1]
SNP_pig_Colon$SNP<-variant_Colon[,2]
SNP_pig_Colon$ref <- variant_Colon[,3]
SNP_pig_Colon$alt <- variant_Colon[,4]

SNP_Colon_pig<-NULL
SNP_Colon_hum<-NULL
for(i in 1:length(overlap_one2oneF_Colon)){
  SNP_hum<-NULL
  SNP_pig<-NULL
  pig_variant<-SNP_pig_Colon[which(SNP_pig_Colon$human_id==overlap_one2oneF_Colon[i]),]
  hum_variant<-SNP_hum_Colon[which(SNP_hum_Colon$gene_id==overlap_one2oneF_Colon[i]),]
  for(j in 1:length(nrow(pig_variant))){
    pair_pig<-paste0(pig_variant$ref[j],pig_variant$alt[j])
    pair_pig_trans<-paste0(pig_variant$alt[j],pig_variant$ref[j])
    for(k in 1:length(nrow(hum_variant))){
      pair_hum<-paste0(hum_variant$ref[k],hum_variant$alt[k])
      if(pair_pig==pair_hum|pair_hum==pair_pig_trans){
        SNP_hum<-rbind(SNP_hum,hum_variant[k,])
        SNP_pig<-rbind(SNP_pig,pig_variant[j,])
      }
    }
    SNP_Colon_pig<-rbind(SNP_Colon_pig,SNP_pig)
    SNP_Colon_hum<-rbind(SNP_Colon_hum,SNP_hum)
  }
}

#Heart#
FM_Heart<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/fine_mapping_results/common tissues/Heart.cis_independent_qtl.txt"))
FM_Heart_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Heart_Atrial_Appendage.v8.egenes.txt"))
FM_Heart_hum2<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Heart_Left_Ventricle.v8.egenes.txt"))
FM_Heart_hum<-rbind(FM_Heart_hum1,FM_Heart_hum2)
FM_Heart_hum$gene_id<-substr(FM_Heart_hum$gene_id,1,15)
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
FMgene_pig_Heart<-annotation$`Gene stable ID`[match(unique(FM_Heart$phenotype_id),annotation$`Pig gene stable ID`)]
FMgene_one2one_Heart<-intersect(FMgene_pig_Heart,colnames(expression))
overlap_one2oneF_Heart<-intersect(FMgene_one2one_Heart,unique(FM_Heart_hum$gene_id))
overlap_one2oneF_Heart_pig<-annotation$`Pig gene stable ID`[match(overlap_one2oneF_Heart,annotation$`Gene stable ID`)]

SNP_pig_Heart<-NULL
for(i in 1:length(overlap_one2oneF_Heart)){
  a<-FM_Heart[which(FM_Heart$phenotype_id==overlap_one2oneF_Heart_pig[i]),]
  a$human_id<-overlap_one2oneF_Heart[i]
  SNP_pig_Heart<-rbind(SNP_pig_Heart,a)
}

SNP_hum_Heart<-NULL
for(i in 1:length(overlap_one2oneF_Heart)){
  a<-FM_Heart_hum[which(FM_Heart_hum$gene_id==overlap_one2oneF_Heart[i]),]
  SNP_hum_Heart<-rbind(SNP_hum_Heart,a)
}

variant_Heart<-strsplit(SNP_pig_Heart$variant_id,split = "_")
variant_Heart<-as.data.frame(variant_Heart)
variant_Heart<-as.data.frame(t(variant_Heart))
SNP_pig_Heart$chr <- variant_Heart[,1]
SNP_pig_Heart$SNP<-variant_Heart[,2]
SNP_pig_Heart$ref <- variant_Heart[,3]
SNP_pig_Heart$alt <- variant_Heart[,4]

SNP_Heart_pig<-NULL
SNP_Heart_hum<-NULL
for(i in 1:length(overlap_one2oneF_Heart)){
  SNP_hum<-NULL
  SNP_pig<-NULL
  pig_variant<-SNP_pig_Heart[which(SNP_pig_Heart$human_id==overlap_one2oneF_Heart[i]),]
  hum_variant<-SNP_hum_Heart[which(SNP_hum_Heart$gene_id==overlap_one2oneF_Heart[i]),]
  for(j in 1:length(nrow(pig_variant))){
    pair_pig<-paste0(pig_variant$ref[j],pig_variant$alt[j])
    pair_pig_trans<-paste0(pig_variant$alt[j],pig_variant$ref[j])
    for(k in 1:length(nrow(hum_variant))){
      pair_hum<-paste0(hum_variant$ref[k],hum_variant$alt[k])
      if(pair_pig==pair_hum|pair_hum==pair_pig_trans){
        SNP_hum<-rbind(SNP_hum,hum_variant[k,])
        SNP_pig<-rbind(SNP_pig,pig_variant[j,])
      }
    }
    SNP_Heart_pig<-rbind(SNP_Heart_pig,SNP_pig)
    SNP_Heart_hum<-rbind(SNP_Heart_hum,SNP_hum)
  }
}

#Frontal_cortex#
FM_Frontal_cortex<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/fine_mapping_results/common tissues/Frontal_cortex.cis_independent_qtl.txt"))
FM_Frontal_cortex_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Brain_Frontal_Cortex_BA9.v8.egenes.txt"))
FM_Frontal_cortex_hum<-rbind(FM_Frontal_cortex_hum1)
FM_Frontal_cortex_hum$gene_id<-substr(FM_Frontal_cortex_hum$gene_id,1,15)
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
FMgene_pig_Frontal_cortex<-annotation$`Gene stable ID`[match(unique(FM_Frontal_cortex$phenotype_id),annotation$`Pig gene stable ID`)]
FMgene_one2one_Frontal_cortex<-intersect(FMgene_pig_Frontal_cortex,colnames(expression))
overlap_one2oneF_Frontal_cortex<-intersect(FMgene_one2one_Frontal_cortex,unique(FM_Frontal_cortex_hum$gene_id))
overlap_one2oneF_Frontal_cortex_pig<-annotation$`Pig gene stable ID`[match(overlap_one2oneF_Frontal_cortex,annotation$`Gene stable ID`)]

SNP_pig_Frontal_cortex<-NULL
for(i in 1:length(overlap_one2oneF_Frontal_cortex)){
  a<-FM_Frontal_cortex[which(FM_Frontal_cortex$phenotype_id==overlap_one2oneF_Frontal_cortex_pig[i]),]
  a$human_id<-overlap_one2oneF_Frontal_cortex[i]
  SNP_pig_Frontal_cortex<-rbind(SNP_pig_Frontal_cortex,a)
}

SNP_hum_Frontal_cortex<-NULL
for(i in 1:length(overlap_one2oneF_Frontal_cortex)){
  a<-FM_Frontal_cortex_hum[which(FM_Frontal_cortex_hum$gene_id==overlap_one2oneF_Frontal_cortex[i]),]
  SNP_hum_Frontal_cortex<-rbind(SNP_hum_Frontal_cortex,a)
}

variant_Frontal_cortex<-strsplit(SNP_pig_Frontal_cortex$variant_id,split = "_")
variant_Frontal_cortex<-as.data.frame(variant_Frontal_cortex)
variant_Frontal_cortex<-as.data.frame(t(variant_Frontal_cortex))
SNP_pig_Frontal_cortex$chr <- variant_Frontal_cortex[,1]
SNP_pig_Frontal_cortex$SNP<-variant_Frontal_cortex[,2]
SNP_pig_Frontal_cortex$ref <- variant_Frontal_cortex[,3]
SNP_pig_Frontal_cortex$alt <- variant_Frontal_cortex[,4]

SNP_Frontal_cortex_pig<-NULL
SNP_Frontal_cortex_hum<-NULL
for(i in 1:length(overlap_one2oneF_Frontal_cortex)){
  SNP_hum<-NULL
  SNP_pig<-NULL
  pig_variant<-SNP_pig_Frontal_cortex[which(SNP_pig_Frontal_cortex$human_id==overlap_one2oneF_Frontal_cortex[i]),]
  hum_variant<-SNP_hum_Frontal_cortex[which(SNP_hum_Frontal_cortex$gene_id==overlap_one2oneF_Frontal_cortex[i]),]
  for(j in 1:length(nrow(pig_variant))){
    pair_pig<-paste0(pig_variant$ref[j],pig_variant$alt[j])
    pair_pig_trans<-paste0(pig_variant$alt[j],pig_variant$ref[j])
    for(k in 1:length(nrow(hum_variant))){
      pair_hum<-paste0(hum_variant$ref[k],hum_variant$alt[k])
      if(pair_pig==pair_hum|pair_hum==pair_pig_trans){
        SNP_hum<-rbind(SNP_hum,hum_variant[k,])
        SNP_pig<-rbind(SNP_pig,pig_variant[j,])
      }
    }
    SNP_Frontal_cortex_pig<-rbind(SNP_Frontal_cortex_pig,SNP_pig)
    SNP_Frontal_cortex_hum<-rbind(SNP_Frontal_cortex_hum,SNP_hum)
  }
}

#Liver#
FM_Liver<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/fine_mapping_results/common tissues/Liver.cis_independent_qtl.txt"))
FM_Liver_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Liver.v8.egenes.txt"))
FM_Liver_hum<-rbind(FM_Liver_hum1)
FM_Liver_hum$gene_id<-substr(FM_Liver_hum$gene_id,1,15)
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
FMgene_pig_Liver<-annotation$`Gene stable ID`[match(unique(FM_Liver$phenotype_id),annotation$`Pig gene stable ID`)]
FMgene_one2one_Liver<-intersect(FMgene_pig_Liver,colnames(expression))
overlap_one2oneF_Liver<-intersect(FMgene_one2one_Liver,unique(FM_Liver_hum$gene_id))
overlap_one2oneF_Liver_pig<-annotation$`Pig gene stable ID`[match(overlap_one2oneF_Liver,annotation$`Gene stable ID`)]

SNP_pig_Liver<-NULL
for(i in 1:length(overlap_one2oneF_Liver)){
  a<-FM_Liver[which(FM_Liver$phenotype_id==overlap_one2oneF_Liver_pig[i]),]
  a$human_id<-overlap_one2oneF_Liver[i]
  SNP_pig_Liver<-rbind(SNP_pig_Liver,a)
}

SNP_hum_Liver<-NULL
for(i in 1:length(overlap_one2oneF_Liver)){
  a<-FM_Liver_hum[which(FM_Liver_hum$gene_id==overlap_one2oneF_Liver[i]),]
  SNP_hum_Liver<-rbind(SNP_hum_Liver,a)
}

variant_Liver<-strsplit(SNP_pig_Liver$variant_id,split = "_")
variant_Liver<-as.data.frame(variant_Liver)
variant_Liver<-as.data.frame(t(variant_Liver))
SNP_pig_Liver$chr <- variant_Liver[,1]
SNP_pig_Liver$SNP<-variant_Liver[,2]
SNP_pig_Liver$ref <- variant_Liver[,3]
SNP_pig_Liver$alt <- variant_Liver[,4]

SNP_Liver_pig<-NULL
SNP_Liver_hum<-NULL
for(i in 1:length(overlap_one2oneF_Liver)){
  SNP_hum<-NULL
  SNP_pig<-NULL
  pig_variant<-SNP_pig_Liver[which(SNP_pig_Liver$human_id==overlap_one2oneF_Liver[i]),]
  hum_variant<-SNP_hum_Liver[which(SNP_hum_Liver$gene_id==overlap_one2oneF_Liver[i]),]
  for(j in 1:length(nrow(pig_variant))){
    pair_pig<-paste0(pig_variant$ref[j],pig_variant$alt[j])
    pair_pig_trans<-paste0(pig_variant$alt[j],pig_variant$ref[j])
    for(k in 1:length(nrow(hum_variant))){
      pair_hum<-paste0(hum_variant$ref[k],hum_variant$alt[k])
      if(pair_pig==pair_hum|pair_hum==pair_pig_trans){
        SNP_hum<-rbind(SNP_hum,hum_variant[k,])
        SNP_pig<-rbind(SNP_pig,pig_variant[j,])
      }
    }
    SNP_Liver_pig<-rbind(SNP_Liver_pig,SNP_pig)
    SNP_Liver_hum<-rbind(SNP_Liver_hum,SNP_hum)
  }
}

#Hypothalamus#
FM_Hypothalamus<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/fine_mapping_results/common tissues/Hypothalamus.cis_independent_qtl.txt"))
FM_Hypothalamus_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Brain_Hypothalamus.v8.egenes.txt"))
FM_Hypothalamus_hum<-rbind(FM_Hypothalamus_hum1)
FM_Hypothalamus_hum$gene_id<-substr(FM_Hypothalamus_hum$gene_id,1,15)
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
FMgene_pig_Hypothalamus<-annotation$`Gene stable ID`[match(unique(FM_Hypothalamus$phenotype_id),annotation$`Pig gene stable ID`)]
FMgene_one2one_Hypothalamus<-intersect(FMgene_pig_Hypothalamus,colnames(expression))
overlap_one2oneF_Hypothalamus<-intersect(FMgene_one2one_Hypothalamus,unique(FM_Hypothalamus_hum$gene_id))
overlap_one2oneF_Hypothalamus_pig<-annotation$`Pig gene stable ID`[match(overlap_one2oneF_Hypothalamus,annotation$`Gene stable ID`)]

SNP_pig_Hypothalamus<-NULL
for(i in 1:length(overlap_one2oneF_Hypothalamus)){
  a<-FM_Hypothalamus[which(FM_Hypothalamus$phenotype_id==overlap_one2oneF_Hypothalamus_pig[i]),]
  a$human_id<-overlap_one2oneF_Hypothalamus[i]
  SNP_pig_Hypothalamus<-rbind(SNP_pig_Hypothalamus,a)
}

SNP_hum_Hypothalamus<-NULL
for(i in 1:length(overlap_one2oneF_Hypothalamus)){
  a<-FM_Hypothalamus_hum[which(FM_Hypothalamus_hum$gene_id==overlap_one2oneF_Hypothalamus[i]),]
  SNP_hum_Hypothalamus<-rbind(SNP_hum_Hypothalamus,a)
}

variant_Hypothalamus<-strsplit(SNP_pig_Hypothalamus$variant_id,split = "_")
variant_Hypothalamus<-as.data.frame(variant_Hypothalamus)
variant_Hypothalamus<-as.data.frame(t(variant_Hypothalamus))
SNP_pig_Hypothalamus$chr <- variant_Hypothalamus[,1]
SNP_pig_Hypothalamus$SNP<-variant_Hypothalamus[,2]
SNP_pig_Hypothalamus$ref <- variant_Hypothalamus[,3]
SNP_pig_Hypothalamus$alt <- variant_Hypothalamus[,4]

SNP_Hypothalamus_pig<-NULL
SNP_Hypothalamus_hum<-NULL
for(i in 1:length(overlap_one2oneF_Hypothalamus)){
  SNP_hum<-NULL
  SNP_pig<-NULL
  pig_variant<-SNP_pig_Hypothalamus[which(SNP_pig_Hypothalamus$human_id==overlap_one2oneF_Hypothalamus[i]),]
  hum_variant<-SNP_hum_Hypothalamus[which(SNP_hum_Hypothalamus$gene_id==overlap_one2oneF_Hypothalamus[i]),]
  for(j in 1:length(nrow(pig_variant))){
    pair_pig<-paste0(pig_variant$ref[j],pig_variant$alt[j])
    pair_pig_trans<-paste0(pig_variant$alt[j],pig_variant$ref[j])
    for(k in 1:length(nrow(hum_variant))){
      pair_hum<-paste0(hum_variant$ref[k],hum_variant$alt[k])
      if(pair_pig==pair_hum|pair_hum==pair_pig_trans){
        SNP_hum<-rbind(SNP_hum,hum_variant[k,])
        SNP_pig<-rbind(SNP_pig,pig_variant[j,])
      }
    }
    SNP_Hypothalamus_pig<-rbind(SNP_Hypothalamus_pig,SNP_pig)
    SNP_Hypothalamus_hum<-rbind(SNP_Hypothalamus_hum,SNP_hum)
  }
}

#Lung#
FM_Lung<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/fine_mapping_results/common tissues/Lung.cis_independent_qtl.txt"))
FM_Lung_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Lung.v8.egenes.txt"))
FM_Lung_hum<-rbind(FM_Lung_hum1)
FM_Lung_hum$gene_id<-substr(FM_Lung_hum$gene_id,1,15)
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
FMgene_pig_Lung<-annotation$`Gene stable ID`[match(unique(FM_Lung$phenotype_id),annotation$`Pig gene stable ID`)]
FMgene_one2one_Lung<-intersect(FMgene_pig_Lung,colnames(expression))
overlap_one2oneF_Lung<-intersect(FMgene_one2one_Lung,unique(FM_Lung_hum$gene_id))
overlap_one2oneF_Lung_pig<-annotation$`Pig gene stable ID`[match(overlap_one2oneF_Lung,annotation$`Gene stable ID`)]

SNP_pig_Lung<-NULL
for(i in 1:length(overlap_one2oneF_Lung)){
  a<-FM_Lung[which(FM_Lung$phenotype_id==overlap_one2oneF_Lung_pig[i]),]
  a$human_id<-overlap_one2oneF_Lung[i]
  SNP_pig_Lung<-rbind(SNP_pig_Lung,a)
}

SNP_hum_Lung<-NULL
for(i in 1:length(overlap_one2oneF_Lung)){
  a<-FM_Lung_hum[which(FM_Lung_hum$gene_id==overlap_one2oneF_Lung[i]),]
  SNP_hum_Lung<-rbind(SNP_hum_Lung,a)
}

variant_Lung<-strsplit(SNP_pig_Lung$variant_id,split = "_")
variant_Lung<-as.data.frame(variant_Lung)
variant_Lung<-as.data.frame(t(variant_Lung))
SNP_pig_Lung$chr <- variant_Lung[,1]
SNP_pig_Lung$SNP<-variant_Lung[,2]
SNP_pig_Lung$ref <- variant_Lung[,3]
SNP_pig_Lung$alt <- variant_Lung[,4]

SNP_Lung_pig<-NULL
SNP_Lung_hum<-NULL
for(i in 1:length(overlap_one2oneF_Lung)){
  SNP_hum<-NULL
  SNP_pig<-NULL
  pig_variant<-SNP_pig_Lung[which(SNP_pig_Lung$human_id==overlap_one2oneF_Lung[i]),]
  hum_variant<-SNP_hum_Lung[which(SNP_hum_Lung$gene_id==overlap_one2oneF_Lung[i]),]
  for(j in 1:length(nrow(pig_variant))){
    pair_pig<-paste0(pig_variant$ref[j],pig_variant$alt[j])
    pair_pig_trans<-paste0(pig_variant$alt[j],pig_variant$ref[j])
    for(k in 1:length(nrow(hum_variant))){
      pair_hum<-paste0(hum_variant$ref[k],hum_variant$alt[k])
      if(pair_pig==pair_hum|pair_hum==pair_pig_trans){
        SNP_hum<-rbind(SNP_hum,hum_variant[k,])
        SNP_pig<-rbind(SNP_pig,pig_variant[j,])
      }
    }
    SNP_Lung_pig<-rbind(SNP_Lung_pig,SNP_pig)
    SNP_Lung_hum<-rbind(SNP_Lung_hum,SNP_hum)
  }
}

#Ileum#
FM_Ileum<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/fine_mapping_results/common tissues/Ileum.cis_independent_qtl.txt"))
FM_Ileum_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Small_Intestine_Terminal_Ileum.v8.egenes.txt"))
FM_Ileum_hum<-rbind(FM_Ileum_hum1)
FM_Ileum_hum$gene_id<-substr(FM_Ileum_hum$gene_id,1,15)
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
FMgene_pig_Ileum<-annotation$`Gene stable ID`[match(unique(FM_Ileum$phenotype_id),annotation$`Pig gene stable ID`)]
FMgene_one2one_Ileum<-intersect(FMgene_pig_Ileum,colnames(expression))
overlap_one2oneF_Ileum<-intersect(FMgene_one2one_Ileum,unique(FM_Ileum_hum$gene_id))
overlap_one2oneF_Ileum_pig<-annotation$`Pig gene stable ID`[match(overlap_one2oneF_Ileum,annotation$`Gene stable ID`)]

SNP_pig_Ileum<-NULL
for(i in 1:length(overlap_one2oneF_Ileum)){
  a<-FM_Ileum[which(FM_Ileum$phenotype_id==overlap_one2oneF_Ileum_pig[i]),]
  a$human_id<-overlap_one2oneF_Ileum[i]
  SNP_pig_Ileum<-rbind(SNP_pig_Ileum,a)
}

SNP_hum_Ileum<-NULL
for(i in 1:length(overlap_one2oneF_Ileum)){
  a<-FM_Ileum_hum[which(FM_Ileum_hum$gene_id==overlap_one2oneF_Ileum[i]),]
  SNP_hum_Ileum<-rbind(SNP_hum_Ileum,a)
}

variant_Ileum<-strsplit(SNP_pig_Ileum$variant_id,split = "_")
variant_Ileum<-as.data.frame(variant_Ileum)
variant_Ileum<-as.data.frame(t(variant_Ileum))
SNP_pig_Ileum$chr <- variant_Ileum[,1]
SNP_pig_Ileum$SNP<-variant_Ileum[,2]
SNP_pig_Ileum$ref <- variant_Ileum[,3]
SNP_pig_Ileum$alt <- variant_Ileum[,4]

SNP_Ileum_pig<-NULL
SNP_Ileum_hum<-NULL
for(i in 1:length(overlap_one2oneF_Ileum)){
  SNP_hum<-NULL
  SNP_pig<-NULL
  pig_variant<-SNP_pig_Ileum[which(SNP_pig_Ileum$human_id==overlap_one2oneF_Ileum[i]),]
  hum_variant<-SNP_hum_Ileum[which(SNP_hum_Ileum$gene_id==overlap_one2oneF_Ileum[i]),]
  for(j in 1:length(nrow(pig_variant))){
    pair_pig<-paste0(pig_variant$ref[j],pig_variant$alt[j])
    pair_pig_trans<-paste0(pig_variant$alt[j],pig_variant$ref[j])
    for(k in 1:length(nrow(hum_variant))){
      pair_hum<-paste0(hum_variant$ref[k],hum_variant$alt[k])
      if(pair_pig==pair_hum|pair_hum==pair_pig_trans){
        SNP_hum<-rbind(SNP_hum,hum_variant[k,])
        SNP_pig<-rbind(SNP_pig,pig_variant[j,])
      }
    }
    SNP_Ileum_pig<-rbind(SNP_Ileum_pig,SNP_pig)
    SNP_Ileum_hum<-rbind(SNP_Ileum_hum,SNP_hum)
  }
}

#Kidney#
FM_Kidney<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/fine_mapping_results/common tissues/Kidney.cis_independent_qtl.txt"))
FM_Kidney_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Kidney_Cortex.v8.egenes.txt"))
FM_Kidney_hum<-rbind(FM_Kidney_hum1)
FM_Kidney_hum$gene_id<-substr(FM_Kidney_hum$gene_id,1,15)
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
FMgene_pig_Kidney<-annotation$`Gene stable ID`[match(unique(FM_Kidney$phenotype_id),annotation$`Pig gene stable ID`)]
FMgene_one2one_Kidney<-intersect(FMgene_pig_Kidney,colnames(expression))
overlap_one2oneF_Kidney<-intersect(FMgene_one2one_Kidney,unique(FM_Kidney_hum$gene_id))
overlap_one2oneF_Kidney_pig<-annotation$`Pig gene stable ID`[match(overlap_one2oneF_Kidney,annotation$`Gene stable ID`)]

SNP_pig_Kidney<-NULL
for(i in 1:length(overlap_one2oneF_Kidney)){
  a<-FM_Kidney[which(FM_Kidney$phenotype_id==overlap_one2oneF_Kidney_pig[i]),]
  a$human_id<-overlap_one2oneF_Kidney[i]
  SNP_pig_Kidney<-rbind(SNP_pig_Kidney,a)
}

SNP_hum_Kidney<-NULL
for(i in 1:length(overlap_one2oneF_Kidney)){
  a<-FM_Kidney_hum[which(FM_Kidney_hum$gene_id==overlap_one2oneF_Kidney[i]),]
  SNP_hum_Kidney<-rbind(SNP_hum_Kidney,a)
}

variant_Kidney<-strsplit(SNP_pig_Kidney$variant_id,split = "_")
variant_Kidney<-as.data.frame(variant_Kidney)
variant_Kidney<-as.data.frame(t(variant_Kidney))
SNP_pig_Kidney$chr <- variant_Kidney[,1]
SNP_pig_Kidney$SNP<-variant_Kidney[,2]
SNP_pig_Kidney$ref <- variant_Kidney[,3]
SNP_pig_Kidney$alt <- variant_Kidney[,4]

SNP_Kidney_pig<-NULL
SNP_Kidney_hum<-NULL
for(i in 1:length(overlap_one2oneF_Kidney)){
  SNP_hum<-NULL
  SNP_pig<-NULL
  pig_variant<-SNP_pig_Kidney[which(SNP_pig_Kidney$human_id==overlap_one2oneF_Kidney[i]),]
  hum_variant<-SNP_hum_Kidney[which(SNP_hum_Kidney$gene_id==overlap_one2oneF_Kidney[i]),]
  for(j in 1:length(nrow(pig_variant))){
    pair_pig<-paste0(pig_variant$ref[j],pig_variant$alt[j])
    pair_pig_trans<-paste0(pig_variant$alt[j],pig_variant$ref[j])
    for(k in 1:length(nrow(hum_variant))){
      pair_hum<-paste0(hum_variant$ref[k],hum_variant$alt[k])
      if(pair_pig==pair_hum|pair_hum==pair_pig_trans){
        SNP_hum<-rbind(SNP_hum,hum_variant[k,])
        SNP_pig<-rbind(SNP_pig,pig_variant[j,])
      }
    }
    SNP_Kidney_pig<-rbind(SNP_Kidney_pig,SNP_pig)
    SNP_Kidney_hum<-rbind(SNP_Kidney_hum,SNP_hum)
  }
}

#Muscle#
FM_Muscle<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/fine_mapping_results/common tissues/Muscle.cis_independent_qtl.txt"))
FM_Muscle_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Muscle_Skeletal.v8.egenes.txt"))
FM_Muscle_hum<-rbind(FM_Muscle_hum1)
FM_Muscle_hum$gene_id<-substr(FM_Muscle_hum$gene_id,1,15)
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
FMgene_pig_Muscle<-annotation$`Gene stable ID`[match(unique(FM_Muscle$phenotype_id),annotation$`Pig gene stable ID`)]
FMgene_one2one_Muscle<-intersect(FMgene_pig_Muscle,colnames(expression))
overlap_one2oneF_Muscle<-intersect(FMgene_one2one_Muscle,unique(FM_Muscle_hum$gene_id))
overlap_one2oneF_Muscle_pig<-annotation$`Pig gene stable ID`[match(overlap_one2oneF_Muscle,annotation$`Gene stable ID`)]

SNP_pig_Muscle<-NULL
for(i in 1:length(overlap_one2oneF_Muscle)){
  a<-FM_Muscle[which(FM_Muscle$phenotype_id==overlap_one2oneF_Muscle_pig[i]),]
  a$human_id<-overlap_one2oneF_Muscle[i]
  SNP_pig_Muscle<-rbind(SNP_pig_Muscle,a)
}

SNP_hum_Muscle<-NULL
for(i in 1:length(overlap_one2oneF_Muscle)){
  a<-FM_Muscle_hum[which(FM_Muscle_hum$gene_id==overlap_one2oneF_Muscle[i]),]
  SNP_hum_Muscle<-rbind(SNP_hum_Muscle,a)
}

variant_Muscle<-strsplit(SNP_pig_Muscle$variant_id,split = "_")
variant_Muscle<-as.data.frame(variant_Muscle)
variant_Muscle<-as.data.frame(t(variant_Muscle))
SNP_pig_Muscle$chr <- variant_Muscle[,1]
SNP_pig_Muscle$SNP<-variant_Muscle[,2]
SNP_pig_Muscle$ref <- variant_Muscle[,3]
SNP_pig_Muscle$alt <- variant_Muscle[,4]

SNP_Muscle_pig<-NULL
SNP_Muscle_hum<-NULL
for(i in 1:length(overlap_one2oneF_Muscle)){
  SNP_hum<-NULL
  SNP_pig<-NULL
  pig_variant<-SNP_pig_Muscle[which(SNP_pig_Muscle$human_id==overlap_one2oneF_Muscle[i]),]
  hum_variant<-SNP_hum_Muscle[which(SNP_hum_Muscle$gene_id==overlap_one2oneF_Muscle[i]),]
  for(j in 1:length(nrow(pig_variant))){
    pair_pig<-paste0(pig_variant$ref[j],pig_variant$alt[j])
    pair_pig_trans<-paste0(pig_variant$alt[j],pig_variant$ref[j])
    for(k in 1:length(nrow(hum_variant))){
      pair_hum<-paste0(hum_variant$ref[k],hum_variant$alt[k])
      if(pair_pig==pair_hum|pair_hum==pair_pig_trans){
        SNP_hum<-rbind(SNP_hum,hum_variant[k,])
        SNP_pig<-rbind(SNP_pig,pig_variant[j,])
      }
    }
    SNP_Muscle_pig<-rbind(SNP_Muscle_pig,SNP_pig)
    SNP_Muscle_hum<-rbind(SNP_Muscle_hum,SNP_hum)
  }
}

#Ovary#
FM_Ovary<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/fine_mapping_results/common tissues/Ovary.cis_independent_qtl.txt"))
FM_Ovary_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Ovary.v8.egenes.txt"))
FM_Ovary_hum<-rbind(FM_Ovary_hum1)
FM_Ovary_hum$gene_id<-substr(FM_Ovary_hum$gene_id,1,15)
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
FMgene_pig_Ovary<-annotation$`Gene stable ID`[match(unique(FM_Ovary$phenotype_id),annotation$`Pig gene stable ID`)]
FMgene_one2one_Ovary<-intersect(FMgene_pig_Ovary,colnames(expression))
overlap_one2oneF_Ovary<-intersect(FMgene_one2one_Ovary,unique(FM_Ovary_hum$gene_id))
overlap_one2oneF_Ovary_pig<-annotation$`Pig gene stable ID`[match(overlap_one2oneF_Ovary,annotation$`Gene stable ID`)]

SNP_pig_Ovary<-NULL
for(i in 1:length(overlap_one2oneF_Ovary)){
  a<-FM_Ovary[which(FM_Ovary$phenotype_id==overlap_one2oneF_Ovary_pig[i]),]
  a$human_id<-overlap_one2oneF_Ovary[i]
  SNP_pig_Ovary<-rbind(SNP_pig_Ovary,a)
}

SNP_hum_Ovary<-NULL
for(i in 1:length(overlap_one2oneF_Ovary)){
  a<-FM_Ovary_hum[which(FM_Ovary_hum$gene_id==overlap_one2oneF_Ovary[i]),]
  SNP_hum_Ovary<-rbind(SNP_hum_Ovary,a)
}

variant_Ovary<-strsplit(SNP_pig_Ovary$variant_id,split = "_")
variant_Ovary<-as.data.frame(variant_Ovary)
variant_Ovary<-as.data.frame(t(variant_Ovary))
SNP_pig_Ovary$chr <- variant_Ovary[,1]
SNP_pig_Ovary$SNP<-variant_Ovary[,2]
SNP_pig_Ovary$ref <- variant_Ovary[,3]
SNP_pig_Ovary$alt <- variant_Ovary[,4]

SNP_Ovary_pig<-NULL
SNP_Ovary_hum<-NULL
for(i in 1:length(overlap_one2oneF_Ovary)){
  SNP_hum<-NULL
  SNP_pig<-NULL
  pig_variant<-SNP_pig_Ovary[which(SNP_pig_Ovary$human_id==overlap_one2oneF_Ovary[i]),]
  hum_variant<-SNP_hum_Ovary[which(SNP_hum_Ovary$gene_id==overlap_one2oneF_Ovary[i]),]
  for(j in 1:length(nrow(pig_variant))){
    pair_pig<-paste0(pig_variant$ref[j],pig_variant$alt[j])
    pair_pig_trans<-paste0(pig_variant$alt[j],pig_variant$ref[j])
    for(k in 1:length(nrow(hum_variant))){
      pair_hum<-paste0(hum_variant$ref[k],hum_variant$alt[k])
      if(pair_pig==pair_hum|pair_hum==pair_pig_trans){
        SNP_hum<-rbind(SNP_hum,hum_variant[k,])
        SNP_pig<-rbind(SNP_pig,pig_variant[j,])
      }
    }
    SNP_Ovary_pig<-rbind(SNP_Ovary_pig,SNP_pig)
    SNP_Ovary_hum<-rbind(SNP_Ovary_hum,SNP_hum)
  }
}

#Pituitary#
FM_Pituitary<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/fine_mapping_results/common tissues/Pituitary.cis_independent_qtl.txt"))
FM_Pituitary_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Pituitary.v8.egenes.txt"))
FM_Pituitary_hum<-rbind(FM_Pituitary_hum1)
FM_Pituitary_hum$gene_id<-substr(FM_Pituitary_hum$gene_id,1,15)
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
FMgene_pig_Pituitary<-annotation$`Gene stable ID`[match(unique(FM_Pituitary$phenotype_id),annotation$`Pig gene stable ID`)]
FMgene_one2one_Pituitary<-intersect(FMgene_pig_Pituitary,colnames(expression))
overlap_one2oneF_Pituitary<-intersect(FMgene_one2one_Pituitary,unique(FM_Pituitary_hum$gene_id))
overlap_one2oneF_Pituitary_pig<-annotation$`Pig gene stable ID`[match(overlap_one2oneF_Pituitary,annotation$`Gene stable ID`)]

SNP_pig_Pituitary<-NULL
for(i in 1:length(overlap_one2oneF_Pituitary)){
  a<-FM_Pituitary[which(FM_Pituitary$phenotype_id==overlap_one2oneF_Pituitary_pig[i]),]
  a$human_id<-overlap_one2oneF_Pituitary[i]
  SNP_pig_Pituitary<-rbind(SNP_pig_Pituitary,a)
}

SNP_hum_Pituitary<-NULL
for(i in 1:length(overlap_one2oneF_Pituitary)){
  a<-FM_Pituitary_hum[which(FM_Pituitary_hum$gene_id==overlap_one2oneF_Pituitary[i]),]
  SNP_hum_Pituitary<-rbind(SNP_hum_Pituitary,a)
}

variant_Pituitary<-strsplit(SNP_pig_Pituitary$variant_id,split = "_")
variant_Pituitary<-as.data.frame(variant_Pituitary)
variant_Pituitary<-as.data.frame(t(variant_Pituitary))
SNP_pig_Pituitary$chr <- variant_Pituitary[,1]
SNP_pig_Pituitary$SNP<-variant_Pituitary[,2]
SNP_pig_Pituitary$ref <- variant_Pituitary[,3]
SNP_pig_Pituitary$alt <- variant_Pituitary[,4]

SNP_Pituitary_pig<-NULL
SNP_Pituitary_hum<-NULL
for(i in 1:length(overlap_one2oneF_Pituitary)){
  SNP_hum<-NULL
  SNP_pig<-NULL
  pig_variant<-SNP_pig_Pituitary[which(SNP_pig_Pituitary$human_id==overlap_one2oneF_Pituitary[i]),]
  hum_variant<-SNP_hum_Pituitary[which(SNP_hum_Pituitary$gene_id==overlap_one2oneF_Pituitary[i]),]
  for(j in 1:length(nrow(pig_variant))){
    pair_pig<-paste0(pig_variant$ref[j],pig_variant$alt[j])
    pair_pig_trans<-paste0(pig_variant$alt[j],pig_variant$ref[j])
    for(k in 1:length(nrow(hum_variant))){
      pair_hum<-paste0(hum_variant$ref[k],hum_variant$alt[k])
      if(pair_pig==pair_hum|pair_hum==pair_pig_trans){
        SNP_hum<-rbind(SNP_hum,hum_variant[k,])
        SNP_pig<-rbind(SNP_pig,pig_variant[j,])
      }
    }
    SNP_Pituitary_pig<-rbind(SNP_Pituitary_pig,SNP_pig)
    SNP_Pituitary_hum<-rbind(SNP_Pituitary_hum,SNP_hum)
  }
}

#Spleen#
FM_Spleen<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/fine_mapping_results/common tissues/Spleen.cis_independent_qtl.txt"))
FM_Spleen_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Spleen.v8.egenes.txt"))
FM_Spleen_hum<-rbind(FM_Spleen_hum1)
FM_Spleen_hum$gene_id<-substr(FM_Spleen_hum$gene_id,1,15)
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
FMgene_pig_Spleen<-annotation$`Gene stable ID`[match(unique(FM_Spleen$phenotype_id),annotation$`Pig gene stable ID`)]
FMgene_one2one_Spleen<-intersect(FMgene_pig_Spleen,colnames(expression))
overlap_one2oneF_Spleen<-intersect(FMgene_one2one_Spleen,unique(FM_Spleen_hum$gene_id))
overlap_one2oneF_Spleen_pig<-annotation$`Pig gene stable ID`[match(overlap_one2oneF_Spleen,annotation$`Gene stable ID`)]

SNP_pig_Spleen<-NULL
for(i in 1:length(overlap_one2oneF_Spleen)){
  a<-FM_Spleen[which(FM_Spleen$phenotype_id==overlap_one2oneF_Spleen_pig[i]),]
  a$human_id<-overlap_one2oneF_Spleen[i]
  SNP_pig_Spleen<-rbind(SNP_pig_Spleen,a)
}

SNP_hum_Spleen<-NULL
for(i in 1:length(overlap_one2oneF_Spleen)){
  a<-FM_Spleen_hum[which(FM_Spleen_hum$gene_id==overlap_one2oneF_Spleen[i]),]
  SNP_hum_Spleen<-rbind(SNP_hum_Spleen,a)
}

variant_Spleen<-strsplit(SNP_pig_Spleen$variant_id,split = "_")
variant_Spleen<-as.data.frame(variant_Spleen)
variant_Spleen<-as.data.frame(t(variant_Spleen))
SNP_pig_Spleen$chr <- variant_Spleen[,1]
SNP_pig_Spleen$SNP<-variant_Spleen[,2]
SNP_pig_Spleen$ref <- variant_Spleen[,3]
SNP_pig_Spleen$alt <- variant_Spleen[,4]

SNP_Spleen_pig<-NULL
SNP_Spleen_hum<-NULL
for(i in 1:length(overlap_one2oneF_Spleen)){
  SNP_hum<-NULL
  SNP_pig<-NULL
  pig_variant<-SNP_pig_Spleen[which(SNP_pig_Spleen$human_id==overlap_one2oneF_Spleen[i]),]
  hum_variant<-SNP_hum_Spleen[which(SNP_hum_Spleen$gene_id==overlap_one2oneF_Spleen[i]),]
  for(j in 1:length(nrow(pig_variant))){
    pair_pig<-paste0(pig_variant$ref[j],pig_variant$alt[j])
    pair_pig_trans<-paste0(pig_variant$alt[j],pig_variant$ref[j])
    for(k in 1:length(nrow(hum_variant))){
      pair_hum<-paste0(hum_variant$ref[k],hum_variant$alt[k])
      if(pair_pig==pair_hum|pair_hum==pair_pig_trans){
        SNP_hum<-rbind(SNP_hum,hum_variant[k,])
        SNP_pig<-rbind(SNP_pig,pig_variant[j,])
      }
    }
    SNP_Spleen_pig<-rbind(SNP_Spleen_pig,SNP_pig)
    SNP_Spleen_hum<-rbind(SNP_Spleen_hum,SNP_hum)
  }
}

#Testis#
FM_Testis<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/fine_mapping_results/common tissues/Testis.cis_independent_qtl.txt"))
FM_Testis_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Testis.v8.egenes.txt"))
FM_Testis_hum<-rbind(FM_Testis_hum1)
FM_Testis_hum$gene_id<-substr(FM_Testis_hum$gene_id,1,15)
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
FMgene_pig_Testis<-annotation$`Gene stable ID`[match(unique(FM_Testis$phenotype_id),annotation$`Pig gene stable ID`)]
FMgene_one2one_Testis<-intersect(FMgene_pig_Testis,colnames(expression))
overlap_one2oneF_Testis<-intersect(FMgene_one2one_Testis,unique(FM_Testis_hum$gene_id))
overlap_one2oneF_Testis_pig<-annotation$`Pig gene stable ID`[match(overlap_one2oneF_Testis,annotation$`Gene stable ID`)]

SNP_pig_Testis<-NULL
for(i in 1:length(overlap_one2oneF_Testis)){
  a<-FM_Testis[which(FM_Testis$phenotype_id==overlap_one2oneF_Testis_pig[i]),]
  a$human_id<-overlap_one2oneF_Testis[i]
  SNP_pig_Testis<-rbind(SNP_pig_Testis,a)
}

SNP_hum_Testis<-NULL
for(i in 1:length(overlap_one2oneF_Testis)){
  a<-FM_Testis_hum[which(FM_Testis_hum$gene_id==overlap_one2oneF_Testis[i]),]
  SNP_hum_Testis<-rbind(SNP_hum_Testis,a)
}

variant_Testis<-strsplit(SNP_pig_Testis$variant_id,split = "_")
variant_Testis<-as.data.frame(variant_Testis)
variant_Testis<-as.data.frame(t(variant_Testis))
SNP_pig_Testis$chr <- variant_Testis[,1]
SNP_pig_Testis$SNP<-variant_Testis[,2]
SNP_pig_Testis$ref <- variant_Testis[,3]
SNP_pig_Testis$alt <- variant_Testis[,4]

SNP_Testis_pig<-NULL
SNP_Testis_hum<-NULL
for(i in 1:length(overlap_one2oneF_Testis)){
  SNP_hum<-NULL
  SNP_pig<-NULL
  pig_variant<-SNP_pig_Testis[which(SNP_pig_Testis$human_id==overlap_one2oneF_Testis[i]),]
  hum_variant<-SNP_hum_Testis[which(SNP_hum_Testis$gene_id==overlap_one2oneF_Testis[i]),]
  for(j in 1:length(nrow(pig_variant))){
    pair_pig<-paste0(pig_variant$ref[j],pig_variant$alt[j])
    pair_pig_trans<-paste0(pig_variant$alt[j],pig_variant$ref[j])
    for(k in 1:length(nrow(hum_variant))){
      pair_hum<-paste0(hum_variant$ref[k],hum_variant$alt[k])
      if(pair_pig==pair_hum|pair_hum==pair_pig_trans){
        SNP_hum<-rbind(SNP_hum,hum_variant[k,])
        SNP_pig<-rbind(SNP_pig,pig_variant[j,])
      }
    }
    SNP_Testis_pig<-rbind(SNP_Testis_pig,SNP_pig)
    SNP_Testis_hum<-rbind(SNP_Testis_hum,SNP_hum)
  }
}

#Uterus#
FM_Uterus<-as.data.frame(fread("/Users/baizhonghao/Downloads/Protein_coding_gene/fine_mapping_results/common tissues/Uterus.cis_independent_qtl.txt"))
FM_Uterus_hum1<-as.data.frame(fread("/Users/baizhonghao/Downloads/GTEx_Analysis_v8_eQTL/Uterus.v8.egenes.txt"))
FM_Uterus_hum<-rbind(FM_Uterus_hum1)
FM_Uterus_hum$gene_id<-substr(FM_Uterus_hum$gene_id,1,15)
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
FMgene_pig_Uterus<-annotation$`Gene stable ID`[match(unique(FM_Uterus$phenotype_id),annotation$`Pig gene stable ID`)]
FMgene_one2one_Uterus<-intersect(FMgene_pig_Uterus,colnames(expression))
overlap_one2oneF_Uterus<-intersect(FMgene_one2one_Uterus,unique(FM_Uterus_hum$gene_id))
overlap_one2oneF_Uterus_pig<-annotation$`Pig gene stable ID`[match(overlap_one2oneF_Uterus,annotation$`Gene stable ID`)]

SNP_pig_Uterus<-NULL
for(i in 1:length(overlap_one2oneF_Uterus)){
  a<-FM_Uterus[which(FM_Uterus$phenotype_id==overlap_one2oneF_Uterus_pig[i]),]
  a$human_id<-overlap_one2oneF_Uterus[i]
  SNP_pig_Uterus<-rbind(SNP_pig_Uterus,a)
}

SNP_hum_Uterus<-NULL
for(i in 1:length(overlap_one2oneF_Uterus)){
  a<-FM_Uterus_hum[which(FM_Uterus_hum$gene_id==overlap_one2oneF_Uterus[i]),]
  SNP_hum_Uterus<-rbind(SNP_hum_Uterus,a)
}

variant_Uterus<-strsplit(SNP_pig_Uterus$variant_id,split = "_")
variant_Uterus<-as.data.frame(variant_Uterus)
variant_Uterus<-as.data.frame(t(variant_Uterus))
SNP_pig_Uterus$chr <- variant_Uterus[,1]
SNP_pig_Uterus$SNP<-variant_Uterus[,2]
SNP_pig_Uterus$ref <- variant_Uterus[,3]
SNP_pig_Uterus$alt <- variant_Uterus[,4]

SNP_Uterus_pig<-NULL
SNP_Uterus_hum<-NULL
for(i in 1:length(overlap_one2oneF_Uterus)){
  SNP_hum<-NULL
  SNP_pig<-NULL
  pig_variant<-SNP_pig_Uterus[which(SNP_pig_Uterus$human_id==overlap_one2oneF_Uterus[i]),]
  hum_variant<-SNP_hum_Uterus[which(SNP_hum_Uterus$gene_id==overlap_one2oneF_Uterus[i]),]
  for(j in 1:length(nrow(pig_variant))){
    pair_pig<-paste0(pig_variant$ref[j],pig_variant$alt[j])
    pair_pig_trans<-paste0(pig_variant$alt[j],pig_variant$ref[j])
    for(k in 1:length(nrow(hum_variant))){
      pair_hum<-paste0(hum_variant$ref[k],hum_variant$alt[k])
      if(pair_pig==pair_hum|pair_hum==pair_pig_trans){
        SNP_hum<-rbind(SNP_hum,hum_variant[k,])
        SNP_pig<-rbind(SNP_pig,pig_variant[j,])
      }
    }
    SNP_Uterus_pig<-rbind(SNP_Uterus_pig,SNP_pig)
    SNP_Uterus_hum<-rbind(SNP_Uterus_hum,SNP_hum)
  }
}

SNP_sum<-array(NA,dim=c(length(tissues),2))
colnames(SNP_sum)<-c("tissues","SNP_number")
SNP_sum<-as.data.frame(SNP_sum)
SNP_sum$tissues<-tissues
for(i in 1:length(tissues)){
  SNP_sum$SNP_number[i]<-nrow(get(paste0("SNP_",tissues[i],"_hum")))
}
SNP_sum$Category<-"overlapSNP"

FM_sum_hum<-array(NA,dim=c(length(tissues),2))
colnames(FM_sum_hum)<-c("tissues","SNP_number")
FM_sum_hum<-as.data.frame(FM_sum_hum)
FM_sum_hum$tissues<-tissues
for(i in 1:length(tissues)){
  FM_sum_hum$SNP_number[i]<-nrow(get(paste0("FM_",tissues[i],"_hum")))
}
FM_sum_hum$Category<-"FM_hum"

FM_sum_pig<-array(NA,dim=c(length(tissues),2))
colnames(FM_sum_pig)<-c("tissues","SNP_number")
FM_sum_pig<-as.data.frame(FM_sum_pig)
FM_sum_pig$tissues<-tissues
for(i in 1:length(tissues)){
  FM_sum_pig$SNP_number[i]<-nrow(get(paste0("FM_",tissues[i])))
}
FM_sum_pig$Category<-"FM_pig"

correspond_hum<-array(NA,dim=c(length(tissues),2))
colnames(correspond_hum)<-c("tissues","SNP_number")
correspond_hum<-as.data.frame(correspond_hum)
correspond_hum$tissues<-tissues
for(i in 1:length(tissues)){
  correspond_hum$SNP_number[i]<-nrow(get(paste0("SNP_hum_",tissues[i])))
}
correspond_hum$Category<-"SNP_hum"

correspond_pig<-array(NA,dim=c(length(tissues),2))
colnames(correspond_pig)<-c("tissues","SNP_number")
correspond_pig<-as.data.frame(correspond_pig)
correspond_pig$tissues<-tissues
for(i in 1:length(tissues)){
  correspond_pig$SNP_number[i]<-nrow(get(paste0("SNP_pig_",tissues[i])))
}
correspond_pig$Category<-"SNP_pig"

SNP_sum<-rbind(FM_sum_hum,FM_sum_pig,correspond_hum,correspond_pig,SNP_sum)

SNP_sum$Category<-factor(SNP_sum$Category,level=c("FM_hum","FM_pig","SNP_hum","SNP_pig","overlapSNP"))

#same loci#
for(i in 1:length(tissues)){
  write.table(get(paste0('SNP_pig_',tissues[i]))[,c(19,20)],file=paste0('/Users/baizhonghao/Downloads/SNP_pig_',tissues[i],'.txt'),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
}
for(i in 1:length(tissues)){
  a<-strsplit(get(paste0("FM_",tissues[i]))$variant_id,split="_")
  a<-as.data.frame(a)
  a<-as.data.frame(t(a))
  write.table(a[,c(1,2)],file=paste0("/Users/baizhonghao/Downloads/SNP_",tissues[i],".txt"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
}


tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/eqtl/human_SNP_barplot.tiff",##reqiured to change
     res = 300, width = 4000, height = 2500,compression = "lzw")
par(mar=c(8,8,6,1))
ggplot(SNP_sum,mapping = aes(x=tissues,y=SNP_number,fill=Category))+
  geom_bar(stat='identity',position="dodge")+theme_bw()+  expand_limits(x = 0, y = 0)+
  labs(x = 'Tissues',y = 'Number of SNP') +
  theme(axis.title =element_text(size = 12),axis.text =element_text(size = 6, color = 'black'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=5.2))+
  theme(axis.title.y.left = element_text(vjust = 2))+
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(axis.text=element_text(color="black"))+
  theme(axis.title=element_text(face = "bold", color="black"))

dev.off()

for(i in 1:length(tissues)){
  i=1
  a<-array(NA,dim=c(nrow(get(paste0("SNP_",tissues[i],"_hum"))),3))
  colnames(a)<-c("ID","ref","alt")
  a<-as.data.frame(a)
  a$ID<-get(paste0("SNP_",tissues[i],"_hum"))$gene_id
  a$ref<-get(paste0("SNP_",tissues[i],"_hum"))$ref
  a$alt<-get(paste0("SNP_",tissues[i],"_hum"))$alt    
  a$pair[which]
}

#loci from human to pig#
#Adipose#
Adipose_chr1<-read.csv("/Users/baizhonghao/Downloads/liftover/Adipose/Adipose_chr1.csv")
Adipose_chr2<-read.csv("/Users/baizhonghao/Downloads/liftover/Adipose/Adipose_chr2.csv")
Adipose_chr3<-read.csv("/Users/baizhonghao/Downloads/liftover/Adipose/Adipose_chr3.csv")
Adipose_chr4<-read.csv("/Users/baizhonghao/Downloads/liftover/Adipose/Adipose_chr4.csv")
Adipose_chr5<-read.csv("/Users/baizhonghao/Downloads/liftover/Adipose/Adipose_chr5.csv")
Adipose_chr6<-read.csv("/Users/baizhonghao/Downloads/liftover/Adipose/Adipose_chr6.csv")
Adipose_chr7<-read.csv("/Users/baizhonghao/Downloads/liftover/Adipose/Adipose_chr7.csv")
Adipose_chr8<-read.csv("/Users/baizhonghao/Downloads/liftover/Adipose/Adipose_chr8.csv")
Adipose_chr9<-read.csv("/Users/baizhonghao/Downloads/liftover/Adipose/Adipose_chr9.csv")
Adipose_chr10<-read.csv("/Users/baizhonghao/Downloads/liftover/Adipose/Adipose_chr10.csv")
Adipose_chr11<-read.csv("/Users/baizhonghao/Downloads/liftover/Adipose/Adipose_chr11.csv")
Adipose_chr12<-read.csv("/Users/baizhonghao/Downloads/liftover/Adipose/Adipose_chr12.csv")
Adipose_chr13<-read.csv("/Users/baizhonghao/Downloads/liftover/Adipose/Adipose_chr13.csv")
Adipose_chr14<-read.csv("/Users/baizhonghao/Downloads/liftover/Adipose/Adipose_chr14.csv")
Adipose_chr15<-read.csv("/Users/baizhonghao/Downloads/liftover/Adipose/Adipose_chr15.csv")
Adipose_chr16<-read.csv("/Users/baizhonghao/Downloads/liftover/Adipose/Adipose_chr16.csv")
Adipose_chr17<-read.csv("/Users/baizhonghao/Downloads/liftover/Adipose/Adipose_chr17.csv")
Adipose_chr18<-read.csv("/Users/baizhonghao/Downloads/liftover/Adipose/Adipose_chr18.csv")
Adipose_loci<-rbind(Adipose_chr1,Adipose_chr2,Adipose_chr3,Adipose_chr4,Adipose_chr5,Adipose_chr6,Adipose_chr7,
                    Adipose_chr8,Adipose_chr9,Adipose_chr10,Adipose_chr11,Adipose_chr12,Adipose_chr13,Adipose_chr14,
                    Adipose_chr15,Adipose_chr16,Adipose_chr17,Adipose_chr18)


#Artery#
Artery_chr1<-read.csv("/Users/baizhonghao/Downloads/liftover/Artery/Artery_chr1.csv")
Artery_chr2<-read.csv("/Users/baizhonghao/Downloads/liftover/Artery/Artery_chr2.csv")
Artery_chr3<-read.csv("/Users/baizhonghao/Downloads/liftover/Artery/Artery_chr3.csv")
Artery_chr4<-read.csv("/Users/baizhonghao/Downloads/liftover/Artery/Artery_chr4.csv")
Artery_chr5<-read.csv("/Users/baizhonghao/Downloads/liftover/Artery/Artery_chr5.csv")
Artery_chr6<-read.csv("/Users/baizhonghao/Downloads/liftover/Artery/Artery_chr6.csv")
Artery_chr7<-read.csv("/Users/baizhonghao/Downloads/liftover/Artery/Artery_chr7.csv")
Artery_chr8<-read.csv("/Users/baizhonghao/Downloads/liftover/Artery/Artery_chr8.csv")
Artery_chr9<-read.csv("/Users/baizhonghao/Downloads/liftover/Artery/Artery_chr9.csv")
Artery_chr10<-read.csv("/Users/baizhonghao/Downloads/liftover/Artery/Artery_chr10.csv")
Artery_chr11<-read.csv("/Users/baizhonghao/Downloads/liftover/Artery/Artery_chr11.csv")
Artery_chr12<-read.csv("/Users/baizhonghao/Downloads/liftover/Artery/Artery_chr12.csv")
Artery_chr13<-read.csv("/Users/baizhonghao/Downloads/liftover/Artery/Artery_chr13.csv")
Artery_chr14<-read.csv("/Users/baizhonghao/Downloads/liftover/Artery/Artery_chr14.csv")
Artery_chr15<-read.csv("/Users/baizhonghao/Downloads/liftover/Artery/Artery_chr15.csv")
Artery_chr16<-read.csv("/Users/baizhonghao/Downloads/liftover/Artery/Artery_chr16.csv")
Artery_chr17<-read.csv("/Users/baizhonghao/Downloads/liftover/Artery/Artery_chr17.csv")
Artery_chr18<-read.csv("/Users/baizhonghao/Downloads/liftover/Artery/Artery_chr18.csv")
Artery_loci<-rbind(Artery_chr1,Artery_chr2,Artery_chr3,Artery_chr4,Artery_chr5,Artery_chr6,Artery_chr7,
                   Artery_chr8,Artery_chr9,Artery_chr10,Artery_chr11,Artery_chr12,Artery_chr13,Artery_chr14,
                   Artery_chr15,Artery_chr16,Artery_chr17,Artery_chr18)

#Blood#
Blood_chr1<-read.csv("/Users/baizhonghao/Downloads/liftover/Blood/Blood_chr1.csv")
Blood_chr2<-read.csv("/Users/baizhonghao/Downloads/liftover/Blood/Blood_chr2.csv")
Blood_chr3<-read.csv("/Users/baizhonghao/Downloads/liftover/Blood/Blood_chr3.csv")
Blood_chr4<-read.csv("/Users/baizhonghao/Downloads/liftover/Blood/Blood_chr4.csv")
Blood_chr5<-read.csv("/Users/baizhonghao/Downloads/liftover/Blood/Blood_chr5.csv")
Blood_chr6<-read.csv("/Users/baizhonghao/Downloads/liftover/Blood/Blood_chr6.csv")
Blood_chr7<-read.csv("/Users/baizhonghao/Downloads/liftover/Blood/Blood_chr7.csv")
Blood_chr8<-read.csv("/Users/baizhonghao/Downloads/liftover/Blood/Blood_chr8.csv")
Blood_chr9<-read.csv("/Users/baizhonghao/Downloads/liftover/Blood/Blood_chr9.csv")
Blood_chr10<-read.csv("/Users/baizhonghao/Downloads/liftover/Blood/Blood_chr10.csv")
Blood_chr11<-read.csv("/Users/baizhonghao/Downloads/liftover/Blood/Blood_chr11.csv")
Blood_chr12<-read.csv("/Users/baizhonghao/Downloads/liftover/Blood/Blood_chr12.csv")
Blood_chr13<-read.csv("/Users/baizhonghao/Downloads/liftover/Blood/Blood_chr13.csv")
Blood_chr14<-read.csv("/Users/baizhonghao/Downloads/liftover/Blood/Blood_chr14.csv")
Blood_chr15<-read.csv("/Users/baizhonghao/Downloads/liftover/Blood/Blood_chr15.csv")
Blood_chr16<-read.csv("/Users/baizhonghao/Downloads/liftover/Blood/Blood_chr16.csv")
Blood_chr17<-read.csv("/Users/baizhonghao/Downloads/liftover/Blood/Blood_chr17.csv")
Blood_chr18<-read.csv("/Users/baizhonghao/Downloads/liftover/Blood/Blood_chr18.csv")
Blood_loci<-rbind(Blood_chr1,Blood_chr2,Blood_chr3,Blood_chr4,Blood_chr5,Blood_chr6,Blood_chr7,
                  Blood_chr8,Blood_chr9,Blood_chr10,Blood_chr11,Blood_chr12,Blood_chr13,Blood_chr14,
                  Blood_chr15,Blood_chr16,Blood_chr17,Blood_chr18)

#Colon#
Colon_chr1<-read.csv("/Users/baizhonghao/Downloads/liftover/Colon/Colon_chr1.csv")
Colon_chr2<-read.csv("/Users/baizhonghao/Downloads/liftover/Colon/Colon_chr2.csv")
Colon_chr3<-read.csv("/Users/baizhonghao/Downloads/liftover/Colon/Colon_chr3.csv")
Colon_chr4<-read.csv("/Users/baizhonghao/Downloads/liftover/Colon/Colon_chr4.csv")
Colon_chr5<-read.csv("/Users/baizhonghao/Downloads/liftover/Colon/Colon_chr5.csv")
Colon_chr6<-read.csv("/Users/baizhonghao/Downloads/liftover/Colon/Colon_chr6.csv")
Colon_chr7<-read.csv("/Users/baizhonghao/Downloads/liftover/Colon/Colon_chr7.csv")
Colon_chr8<-read.csv("/Users/baizhonghao/Downloads/liftover/Colon/Colon_chr8.csv")
Colon_chr9<-read.csv("/Users/baizhonghao/Downloads/liftover/Colon/Colon_chr9.csv")
Colon_chr10<-read.csv("/Users/baizhonghao/Downloads/liftover/Colon/Colon_chr10.csv")
Colon_chr11<-read.csv("/Users/baizhonghao/Downloads/liftover/Colon/Colon_chr11.csv")
Colon_chr12<-read.csv("/Users/baizhonghao/Downloads/liftover/Colon/Colon_chr12.csv")
Colon_chr13<-read.csv("/Users/baizhonghao/Downloads/liftover/Colon/Colon_chr13.csv")
Colon_chr14<-read.csv("/Users/baizhonghao/Downloads/liftover/Colon/Colon_chr14.csv")
Colon_chr15<-read.csv("/Users/baizhonghao/Downloads/liftover/Colon/Colon_chr15.csv")
Colon_chr16<-read.csv("/Users/baizhonghao/Downloads/liftover/Colon/Colon_chr16.csv")
Colon_chr17<-read.csv("/Users/baizhonghao/Downloads/liftover/Colon/Colon_chr17.csv")
Colon_chr18<-read.csv("/Users/baizhonghao/Downloads/liftover/Colon/Colon_chr18.csv")
Colon_loci<-rbind(Colon_chr1,Colon_chr2,Colon_chr3,Colon_chr4,Colon_chr5,Colon_chr6,Colon_chr7,
                  Colon_chr8,Colon_chr9,Colon_chr10,Colon_chr11,Colon_chr12,Colon_chr13,Colon_chr14,
                  Colon_chr15,Colon_chr16,Colon_chr17,Colon_chr18)

#Frontal_cortex#
Frontal_cortex_chr1<-read.csv("/Users/baizhonghao/Downloads/liftover/Frontal_cortex/Frontal_cortex_chr1.csv")
Frontal_cortex_chr2<-read.csv("/Users/baizhonghao/Downloads/liftover/Frontal_cortex/Frontal_cortex_chr2.csv")
Frontal_cortex_chr3<-read.csv("/Users/baizhonghao/Downloads/liftover/Frontal_cortex/Frontal_cortex_chr3.csv")
Frontal_cortex_chr4<-read.csv("/Users/baizhonghao/Downloads/liftover/Frontal_cortex/Frontal_cortex_chr4.csv")
Frontal_cortex_chr5<-read.csv("/Users/baizhonghao/Downloads/liftover/Frontal_cortex/Frontal_cortex_chr5.csv")
Frontal_cortex_chr6<-read.csv("/Users/baizhonghao/Downloads/liftover/Frontal_cortex/Frontal_cortex_chr6.csv")
Frontal_cortex_chr7<-read.csv("/Users/baizhonghao/Downloads/liftover/Frontal_cortex/Frontal_cortex_chr7.csv")
Frontal_cortex_chr8<-read.csv("/Users/baizhonghao/Downloads/liftover/Frontal_cortex/Frontal_cortex_chr8.csv")
Frontal_cortex_chr9<-read.csv("/Users/baizhonghao/Downloads/liftover/Frontal_cortex/Frontal_cortex_chr9.csv")
Frontal_cortex_chr10<-read.csv("/Users/baizhonghao/Downloads/liftover/Frontal_cortex/Frontal_cortex_chr10.csv")
Frontal_cortex_chr11<-read.csv("/Users/baizhonghao/Downloads/liftover/Frontal_cortex/Frontal_cortex_chr11.csv")
Frontal_cortex_chr12<-read.csv("/Users/baizhonghao/Downloads/liftover/Frontal_cortex/Frontal_cortex_chr12.csv")
Frontal_cortex_chr13<-read.csv("/Users/baizhonghao/Downloads/liftover/Frontal_cortex/Frontal_cortex_chr13.csv")
Frontal_cortex_chr14<-read.csv("/Users/baizhonghao/Downloads/liftover/Frontal_cortex/Frontal_cortex_chr14.csv")
Frontal_cortex_chr15<-read.csv("/Users/baizhonghao/Downloads/liftover/Frontal_cortex/Frontal_cortex_chr15.csv")
Frontal_cortex_chr16<-read.csv("/Users/baizhonghao/Downloads/liftover/Frontal_cortex/Frontal_cortex_chr16.csv")
Frontal_cortex_chr17<-read.csv("/Users/baizhonghao/Downloads/liftover/Frontal_cortex/Frontal_cortex_chr17.csv")
Frontal_cortex_chr18<-read.csv("/Users/baizhonghao/Downloads/liftover/Frontal_cortex/Frontal_cortex_chr18.csv")
Frontal_cortex_loci<-rbind(Frontal_cortex_chr1,Frontal_cortex_chr2,Frontal_cortex_chr3,Frontal_cortex_chr4,Frontal_cortex_chr5,Frontal_cortex_chr6,Frontal_cortex_chr7,
                           Frontal_cortex_chr8,Frontal_cortex_chr9,Frontal_cortex_chr10,Frontal_cortex_chr11,Frontal_cortex_chr12,Frontal_cortex_chr13,Frontal_cortex_chr14,
                           Frontal_cortex_chr15,Frontal_cortex_chr16,Frontal_cortex_chr17,Frontal_cortex_chr18)

#Heart#
Heart_chr1<-read.csv("/Users/baizhonghao/Downloads/liftover/Heart/Heart_chr1.csv")
Heart_chr2<-read.csv("/Users/baizhonghao/Downloads/liftover/Heart/Heart_chr2.csv")
Heart_chr3<-read.csv("/Users/baizhonghao/Downloads/liftover/Heart/Heart_chr3.csv")
Heart_chr4<-read.csv("/Users/baizhonghao/Downloads/liftover/Heart/Heart_chr4.csv")
Heart_chr5<-read.csv("/Users/baizhonghao/Downloads/liftover/Heart/Heart_chr5.csv")
Heart_chr6<-read.csv("/Users/baizhonghao/Downloads/liftover/Heart/Heart_chr6.csv")
Heart_chr7<-read.csv("/Users/baizhonghao/Downloads/liftover/Heart/Heart_chr7.csv")
Heart_chr8<-read.csv("/Users/baizhonghao/Downloads/liftover/Heart/Heart_chr8.csv")
Heart_chr9<-read.csv("/Users/baizhonghao/Downloads/liftover/Heart/Heart_chr9.csv")
Heart_chr10<-read.csv("/Users/baizhonghao/Downloads/liftover/Heart/Heart_chr10.csv")
Heart_chr11<-read.csv("/Users/baizhonghao/Downloads/liftover/Heart/Heart_chr11.csv")
Heart_chr12<-read.csv("/Users/baizhonghao/Downloads/liftover/Heart/Heart_chr12.csv")
Heart_chr13<-read.csv("/Users/baizhonghao/Downloads/liftover/Heart/Heart_chr13.csv")
Heart_chr14<-read.csv("/Users/baizhonghao/Downloads/liftover/Heart/Heart_chr14.csv")
Heart_chr15<-read.csv("/Users/baizhonghao/Downloads/liftover/Heart/Heart_chr15.csv")
Heart_chr16<-read.csv("/Users/baizhonghao/Downloads/liftover/Heart/Heart_chr16.csv")
Heart_chr17<-read.csv("/Users/baizhonghao/Downloads/liftover/Heart/Heart_chr17.csv")
Heart_chr18<-read.csv("/Users/baizhonghao/Downloads/liftover/Heart/Heart_chr18.csv")
Heart_loci<-rbind(Heart_chr1,Heart_chr2,Heart_chr3,Heart_chr4,Heart_chr5,Heart_chr6,Heart_chr7,
                  Heart_chr8,Heart_chr9,Heart_chr10,Heart_chr11,Heart_chr12,Heart_chr13,Heart_chr14,
                  Heart_chr15,Heart_chr16,Heart_chr17,Heart_chr18)

#Hypothalamus#
Hypothalamus_chr1<-read.csv("/Users/baizhonghao/Downloads/liftover/Hypothalamus/Hypothalamus_chr1.csv")
Hypothalamus_chr2<-read.csv("/Users/baizhonghao/Downloads/liftover/Hypothalamus/Hypothalamus_chr2.csv")
Hypothalamus_chr3<-read.csv("/Users/baizhonghao/Downloads/liftover/Hypothalamus/Hypothalamus_chr3.csv")
Hypothalamus_chr4<-read.csv("/Users/baizhonghao/Downloads/liftover/Hypothalamus/Hypothalamus_chr4.csv")
Hypothalamus_chr5<-read.csv("/Users/baizhonghao/Downloads/liftover/Hypothalamus/Hypothalamus_chr5.csv")
Hypothalamus_chr6<-read.csv("/Users/baizhonghao/Downloads/liftover/Hypothalamus/Hypothalamus_chr6.csv")
Hypothalamus_chr7<-read.csv("/Users/baizhonghao/Downloads/liftover/Hypothalamus/Hypothalamus_chr7.csv")
Hypothalamus_chr8<-read.csv("/Users/baizhonghao/Downloads/liftover/Hypothalamus/Hypothalamus_chr8.csv")
Hypothalamus_chr9<-read.csv("/Users/baizhonghao/Downloads/liftover/Hypothalamus/Hypothalamus_chr9.csv")
Hypothalamus_chr10<-read.csv("/Users/baizhonghao/Downloads/liftover/Hypothalamus/Hypothalamus_chr10.csv")
Hypothalamus_chr11<-read.csv("/Users/baizhonghao/Downloads/liftover/Hypothalamus/Hypothalamus_chr11.csv")
Hypothalamus_chr12<-read.csv("/Users/baizhonghao/Downloads/liftover/Hypothalamus/Hypothalamus_chr12.csv")
Hypothalamus_chr13<-read.csv("/Users/baizhonghao/Downloads/liftover/Hypothalamus/Hypothalamus_chr13.csv")
Hypothalamus_chr14<-read.csv("/Users/baizhonghao/Downloads/liftover/Hypothalamus/Hypothalamus_chr14.csv")
Hypothalamus_chr15<-read.csv("/Users/baizhonghao/Downloads/liftover/Hypothalamus/Hypothalamus_chr15.csv")
Hypothalamus_chr16<-read.csv("/Users/baizhonghao/Downloads/liftover/Hypothalamus/Hypothalamus_chr16.csv")
Hypothalamus_chr17<-read.csv("/Users/baizhonghao/Downloads/liftover/Hypothalamus/Hypothalamus_chr17.csv")
Hypothalamus_chr18<-read.csv("/Users/baizhonghao/Downloads/liftover/Hypothalamus/Hypothalamus_chr18.csv")
Hypothalamus_loci<-rbind(Hypothalamus_chr1,Hypothalamus_chr2,Hypothalamus_chr3,Hypothalamus_chr4,Hypothalamus_chr5,Hypothalamus_chr6,Hypothalamus_chr7,
                         Hypothalamus_chr8,Hypothalamus_chr9,Hypothalamus_chr10,Hypothalamus_chr11,Hypothalamus_chr12,Hypothalamus_chr13,Hypothalamus_chr14,
                         Hypothalamus_chr15,Hypothalamus_chr16,Hypothalamus_chr17,Hypothalamus_chr18)

#Ileum#
Ileum_chr1<-read.csv("/Users/baizhonghao/Downloads/liftover/Ileum/Ileum_chr1.csv")
Ileum_chr2<-read.csv("/Users/baizhonghao/Downloads/liftover/Ileum/Ileum_chr2.csv")
Ileum_chr3<-read.csv("/Users/baizhonghao/Downloads/liftover/Ileum/Ileum_chr3.csv")
Ileum_chr4<-read.csv("/Users/baizhonghao/Downloads/liftover/Ileum/Ileum_chr4.csv")
Ileum_chr5<-read.csv("/Users/baizhonghao/Downloads/liftover/Ileum/Ileum_chr5.csv")
Ileum_chr6<-read.csv("/Users/baizhonghao/Downloads/liftover/Ileum/Ileum_chr6.csv")
Ileum_chr7<-read.csv("/Users/baizhonghao/Downloads/liftover/Ileum/Ileum_chr7.csv")
Ileum_chr8<-read.csv("/Users/baizhonghao/Downloads/liftover/Ileum/Ileum_chr8.csv")
Ileum_chr9<-read.csv("/Users/baizhonghao/Downloads/liftover/Ileum/Ileum_chr9.csv")
Ileum_chr10<-read.csv("/Users/baizhonghao/Downloads/liftover/Ileum/Ileum_chr10.csv")
Ileum_chr11<-read.csv("/Users/baizhonghao/Downloads/liftover/Ileum/Ileum_chr11.csv")
Ileum_chr12<-read.csv("/Users/baizhonghao/Downloads/liftover/Ileum/Ileum_chr12.csv")
Ileum_chr13<-read.csv("/Users/baizhonghao/Downloads/liftover/Ileum/Ileum_chr13.csv")
Ileum_chr14<-read.csv("/Users/baizhonghao/Downloads/liftover/Ileum/Ileum_chr14.csv")
Ileum_chr15<-read.csv("/Users/baizhonghao/Downloads/liftover/Ileum/Ileum_chr15.csv")
Ileum_chr16<-read.csv("/Users/baizhonghao/Downloads/liftover/Ileum/Ileum_chr16.csv")
Ileum_chr17<-read.csv("/Users/baizhonghao/Downloads/liftover/Ileum/Ileum_chr17.csv")
Ileum_chr18<-read.csv("/Users/baizhonghao/Downloads/liftover/Ileum/Ileum_chr18.csv")
Ileum_loci<-rbind(Ileum_chr1,Ileum_chr2,Ileum_chr3,Ileum_chr4,Ileum_chr5,Ileum_chr6,Ileum_chr7,
                  Ileum_chr8,Ileum_chr9,Ileum_chr10,Ileum_chr11,Ileum_chr12,Ileum_chr13,Ileum_chr14,
                  Ileum_chr15,Ileum_chr16,Ileum_chr17,Ileum_chr18)

#Kidney#
Kidney_chr1<-read.csv("/Users/baizhonghao/Downloads/liftover/Kidney/Kidney_chr1.csv")
Kidney_chr2<-read.csv("/Users/baizhonghao/Downloads/liftover/Kidney/Kidney_chr2.csv")
Kidney_chr3<-read.csv("/Users/baizhonghao/Downloads/liftover/Kidney/Kidney_chr3.csv")
Kidney_chr4<-read.csv("/Users/baizhonghao/Downloads/liftover/Kidney/Kidney_chr4.csv")
Kidney_chr5<-read.csv("/Users/baizhonghao/Downloads/liftover/Kidney/Kidney_chr5.csv")
Kidney_chr6<-read.csv("/Users/baizhonghao/Downloads/liftover/Kidney/Kidney_chr6.csv")
Kidney_chr7<-read.csv("/Users/baizhonghao/Downloads/liftover/Kidney/Kidney_chr7.csv")
Kidney_chr8<-read.csv("/Users/baizhonghao/Downloads/liftover/Kidney/Kidney_chr8.csv")
Kidney_chr9<-read.csv("/Users/baizhonghao/Downloads/liftover/Kidney/Kidney_chr9.csv")
Kidney_chr10<-read.csv("/Users/baizhonghao/Downloads/liftover/Kidney/Kidney_chr10.csv")
Kidney_chr11<-read.csv("/Users/baizhonghao/Downloads/liftover/Kidney/Kidney_chr11.csv")
Kidney_chr12<-read.csv("/Users/baizhonghao/Downloads/liftover/Kidney/Kidney_chr12.csv")
Kidney_chr13<-read.csv("/Users/baizhonghao/Downloads/liftover/Kidney/Kidney_chr13.csv")
Kidney_chr14<-read.csv("/Users/baizhonghao/Downloads/liftover/Kidney/Kidney_chr14.csv")
Kidney_chr15<-read.csv("/Users/baizhonghao/Downloads/liftover/Kidney/Kidney_chr15.csv")
Kidney_chr16<-read.csv("/Users/baizhonghao/Downloads/liftover/Kidney/Kidney_chr16.csv")
Kidney_chr17<-read.csv("/Users/baizhonghao/Downloads/liftover/Kidney/Kidney_chr17.csv")
Kidney_chr18<-read.csv("/Users/baizhonghao/Downloads/liftover/Kidney/Kidney_chr18.csv")
Kidney_loci<-rbind(Kidney_chr1,Kidney_chr2,Kidney_chr3,Kidney_chr4,Kidney_chr5,Kidney_chr6,Kidney_chr7,
                   Kidney_chr8,Kidney_chr9,Kidney_chr10,Kidney_chr11,Kidney_chr12,Kidney_chr13,Kidney_chr14,
                   Kidney_chr15,Kidney_chr16,Kidney_chr17,Kidney_chr18)

#Liver#
Liver_chr1<-read.csv("/Users/baizhonghao/Downloads/liftover/Liver/Liver_chr1.csv")
Liver_chr2<-read.csv("/Users/baizhonghao/Downloads/liftover/Liver/Liver_chr2.csv")
Liver_chr3<-read.csv("/Users/baizhonghao/Downloads/liftover/Liver/Liver_chr3.csv")
Liver_chr4<-read.csv("/Users/baizhonghao/Downloads/liftover/Liver/Liver_chr4.csv")
Liver_chr5<-read.csv("/Users/baizhonghao/Downloads/liftover/Liver/Liver_chr5.csv")
Liver_chr6<-read.csv("/Users/baizhonghao/Downloads/liftover/Liver/Liver_chr6.csv")
Liver_chr7<-read.csv("/Users/baizhonghao/Downloads/liftover/Liver/Liver_chr7.csv")
Liver_chr8<-read.csv("/Users/baizhonghao/Downloads/liftover/Liver/Liver_chr8.csv")
Liver_chr9<-read.csv("/Users/baizhonghao/Downloads/liftover/Liver/Liver_chr9.csv")
Liver_chr10<-read.csv("/Users/baizhonghao/Downloads/liftover/Liver/Liver_chr10.csv")
Liver_chr11<-read.csv("/Users/baizhonghao/Downloads/liftover/Liver/Liver_chr11.csv")
Liver_chr12<-read.csv("/Users/baizhonghao/Downloads/liftover/Liver/Liver_chr12.csv")
Liver_chr13<-read.csv("/Users/baizhonghao/Downloads/liftover/Liver/Liver_chr13.csv")
Liver_chr14<-read.csv("/Users/baizhonghao/Downloads/liftover/Liver/Liver_chr14.csv")
Liver_chr15<-read.csv("/Users/baizhonghao/Downloads/liftover/Liver/Liver_chr15.csv")
Liver_chr16<-read.csv("/Users/baizhonghao/Downloads/liftover/Liver/Liver_chr16.csv")
Liver_chr17<-read.csv("/Users/baizhonghao/Downloads/liftover/Liver/Liver_chr17.csv")
Liver_chr18<-read.csv("/Users/baizhonghao/Downloads/liftover/Liver/Liver_chr18.csv")
Liver_loci<-rbind(Liver_chr1,Liver_chr2,Liver_chr3,Liver_chr4,Liver_chr5,Liver_chr6,Liver_chr7,
                  Liver_chr8,Liver_chr9,Liver_chr10,Liver_chr11,Liver_chr12,Liver_chr13,Liver_chr14,
                  Liver_chr15,Liver_chr16,Liver_chr17,Liver_chr18)

#Lung#
Lung_chr1<-read.csv("/Users/baizhonghao/Downloads/liftover/Lung/Lung_chr1.csv")
Lung_chr2<-read.csv("/Users/baizhonghao/Downloads/liftover/Lung/Lung_chr2.csv")
Lung_chr3<-read.csv("/Users/baizhonghao/Downloads/liftover/Lung/Lung_chr3.csv")
Lung_chr4<-read.csv("/Users/baizhonghao/Downloads/liftover/Lung/Lung_chr4.csv")
Lung_chr5<-read.csv("/Users/baizhonghao/Downloads/liftover/Lung/Lung_chr5.csv")
Lung_chr6<-read.csv("/Users/baizhonghao/Downloads/liftover/Lung/Lung_chr6.csv")
Lung_chr7<-read.csv("/Users/baizhonghao/Downloads/liftover/Lung/Lung_chr7.csv")
Lung_chr8<-read.csv("/Users/baizhonghao/Downloads/liftover/Lung/Lung_chr8.csv")
Lung_chr9<-read.csv("/Users/baizhonghao/Downloads/liftover/Lung/Lung_chr9.csv")
Lung_chr10<-read.csv("/Users/baizhonghao/Downloads/liftover/Lung/Lung_chr10.csv")
Lung_chr11<-read.csv("/Users/baizhonghao/Downloads/liftover/Lung/Lung_chr11.csv")
Lung_chr12<-read.csv("/Users/baizhonghao/Downloads/liftover/Lung/Lung_chr12.csv")
Lung_chr13<-read.csv("/Users/baizhonghao/Downloads/liftover/Lung/Lung_chr13.csv")
Lung_chr14<-read.csv("/Users/baizhonghao/Downloads/liftover/Lung/Lung_chr14.csv")
Lung_chr15<-read.csv("/Users/baizhonghao/Downloads/liftover/Lung/Lung_chr15.csv")
Lung_chr16<-read.csv("/Users/baizhonghao/Downloads/liftover/Lung/Lung_chr16.csv")
Lung_chr17<-read.csv("/Users/baizhonghao/Downloads/liftover/Lung/Lung_chr17.csv")
Lung_chr18<-read.csv("/Users/baizhonghao/Downloads/liftover/Lung/Lung_chr18.csv")
Lung_loci<-rbind(Lung_chr1,Lung_chr2,Lung_chr3,Lung_chr4,Lung_chr5,Lung_chr6,Lung_chr7,
                 Lung_chr8,Lung_chr9,Lung_chr10,Lung_chr11,Lung_chr12,Lung_chr13,Lung_chr14,
                 Lung_chr15,Lung_chr16,Lung_chr17,Lung_chr18)

#Muscle#
Muscle_chr1<-read.csv("/Users/baizhonghao/Downloads/liftover/Muscle/Muscle_chr1.csv")
Muscle_chr2<-read.csv("/Users/baizhonghao/Downloads/liftover/Muscle/Muscle_chr2.csv")
Muscle_chr3<-read.csv("/Users/baizhonghao/Downloads/liftover/Muscle/Muscle_chr3.csv")
Muscle_chr4<-read.csv("/Users/baizhonghao/Downloads/liftover/Muscle/Muscle_chr4.csv")
Muscle_chr5<-read.csv("/Users/baizhonghao/Downloads/liftover/Muscle/Muscle_chr5.csv")
Muscle_chr6<-read.csv("/Users/baizhonghao/Downloads/liftover/Muscle/Muscle_chr6.csv")
Muscle_chr7<-read.csv("/Users/baizhonghao/Downloads/liftover/Muscle/Muscle_chr7.csv")
Muscle_chr8<-read.csv("/Users/baizhonghao/Downloads/liftover/Muscle/Muscle_chr8.csv")
Muscle_chr9<-read.csv("/Users/baizhonghao/Downloads/liftover/Muscle/Muscle_chr9.csv")
Muscle_chr10<-read.csv("/Users/baizhonghao/Downloads/liftover/Muscle/Muscle_chr10.csv")
Muscle_chr11<-read.csv("/Users/baizhonghao/Downloads/liftover/Muscle/Muscle_chr11.csv")
Muscle_chr12<-read.csv("/Users/baizhonghao/Downloads/liftover/Muscle/Muscle_chr12.csv")
Muscle_chr13<-read.csv("/Users/baizhonghao/Downloads/liftover/Muscle/Muscle_chr13.csv")
Muscle_chr14<-read.csv("/Users/baizhonghao/Downloads/liftover/Muscle/Muscle_chr14.csv")
Muscle_chr15<-read.csv("/Users/baizhonghao/Downloads/liftover/Muscle/Muscle_chr15.csv")
Muscle_chr16<-read.csv("/Users/baizhonghao/Downloads/liftover/Muscle/Muscle_chr16.csv")
Muscle_chr17<-read.csv("/Users/baizhonghao/Downloads/liftover/Muscle/Muscle_chr17.csv")
Muscle_chr18<-read.csv("/Users/baizhonghao/Downloads/liftover/Muscle/Muscle_chr18.csv")
Muscle_loci<-rbind(Muscle_chr1,Muscle_chr2,Muscle_chr3,Muscle_chr4,Muscle_chr5,Muscle_chr6,Muscle_chr7,
                   Muscle_chr8,Muscle_chr9,Muscle_chr10,Muscle_chr11,Muscle_chr12,Muscle_chr13,Muscle_chr14,
                   Muscle_chr15,Muscle_chr16,Muscle_chr17,Muscle_chr18)

eqtl_Muscle_chr1<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Muscle/Muscle.cis_qtl_pairs.1.txt"))
eqtl_Muscle_chr2<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Muscle/Muscle.cis_qtl_pairs.2.txt"))
eqtl_Muscle_chr3<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Muscle/Muscle.cis_qtl_pairs.3.txt"))
eqtl_Muscle_chr4<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Muscle/Muscle.cis_qtl_pairs.4.txt"))
eqtl_Muscle_chr5<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Muscle/Muscle.cis_qtl_pairs.5.txt"))
eqtl_Muscle_chr6<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Muscle/Muscle.cis_qtl_pairs.6.txt"))
eqtl_Muscle_chr7<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Muscle/Muscle.cis_qtl_pairs.7.txt"))
eqtl_Muscle_chr8<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Muscle/Muscle.cis_qtl_pairs.8.txt"))
eqtl_Muscle_chr9<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Muscle/Muscle.cis_qtl_pairs.9.txt"))
eqtl_Muscle_chr10<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Muscle/Muscle.cis_qtl_pairs.10.txt"))
eqtl_Muscle_chr11<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Muscle/Muscle.cis_qtl_pairs.11.txt"))
eqtl_Muscle_chr12<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Muscle/Muscle.cis_qtl_pairs.12.txt"))
eqtl_Muscle_chr13<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Muscle/Muscle.cis_qtl_pairs.13.txt"))
eqtl_Muscle_chr14<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Muscle/Muscle.cis_qtl_pairs.14.txt"))
eqtl_Muscle_chr15<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Muscle/Muscle.cis_qtl_pairs.15.txt"))
eqtl_Muscle_chr16<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Muscle/Muscle.cis_qtl_pairs.16.txt"))
eqtl_Muscle_chr17<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Muscle/Muscle.cis_qtl_pairs.17.txt"))
eqtl_Muscle_chr18<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Muscle/Muscle.cis_qtl_pairs.18.txt"))
eqtl_Muscle_loci<-rbind(eqtl_Muscle_chr1,eqtl_Muscle_chr2,eqtl_Muscle_chr3,eqtl_Muscle_chr4,eqtl_Muscle_chr5,eqtl_Muscle_chr6,eqtl_Muscle_chr7,
                        eqtl_Muscle_chr8,eqtl_Muscle_chr9,eqtl_Muscle_chr10,eqtl_Muscle_chr11,eqtl_Muscle_chr12,eqtl_Muscle_chr13,eqtl_Muscle_chr14,
                        eqtl_Muscle_chr15,eqtl_Muscle_chr16,eqtl_Muscle_chr17,eqtl_Muscle_chr18)


#Ovary#
Ovary_chr1<-read.csv("/Users/baizhonghao/Downloads/liftover/Ovary/Ovary_chr1.csv")
Ovary_chr2<-read.csv("/Users/baizhonghao/Downloads/liftover/Ovary/Ovary_chr2.csv")
Ovary_chr3<-read.csv("/Users/baizhonghao/Downloads/liftover/Ovary/Ovary_chr3.csv")
Ovary_chr4<-read.csv("/Users/baizhonghao/Downloads/liftover/Ovary/Ovary_chr4.csv")
Ovary_chr5<-read.csv("/Users/baizhonghao/Downloads/liftover/Ovary/Ovary_chr5.csv")
Ovary_chr6<-read.csv("/Users/baizhonghao/Downloads/liftover/Ovary/Ovary_chr6.csv")
Ovary_chr7<-read.csv("/Users/baizhonghao/Downloads/liftover/Ovary/Ovary_chr7.csv")
Ovary_chr8<-read.csv("/Users/baizhonghao/Downloads/liftover/Ovary/Ovary_chr8.csv")
Ovary_chr9<-read.csv("/Users/baizhonghao/Downloads/liftover/Ovary/Ovary_chr9.csv")
Ovary_chr10<-read.csv("/Users/baizhonghao/Downloads/liftover/Ovary/Ovary_chr10.csv")
Ovary_chr11<-read.csv("/Users/baizhonghao/Downloads/liftover/Ovary/Ovary_chr11.csv")
Ovary_chr12<-read.csv("/Users/baizhonghao/Downloads/liftover/Ovary/Ovary_chr12.csv")
Ovary_chr13<-read.csv("/Users/baizhonghao/Downloads/liftover/Ovary/Ovary_chr13.csv")
Ovary_chr14<-read.csv("/Users/baizhonghao/Downloads/liftover/Ovary/Ovary_chr14.csv")
Ovary_chr15<-read.csv("/Users/baizhonghao/Downloads/liftover/Ovary/Ovary_chr15.csv")
Ovary_chr16<-read.csv("/Users/baizhonghao/Downloads/liftover/Ovary/Ovary_chr16.csv")
Ovary_chr17<-read.csv("/Users/baizhonghao/Downloads/liftover/Ovary/Ovary_chr17.csv")
Ovary_chr18<-read.csv("/Users/baizhonghao/Downloads/liftover/Ovary/Ovary_chr18.csv")
Ovary_loci<-rbind(Ovary_chr1,Ovary_chr2,Ovary_chr3,Ovary_chr4,Ovary_chr5,Ovary_chr6,Ovary_chr7,
                  Ovary_chr8,Ovary_chr9,Ovary_chr10,Ovary_chr11,Ovary_chr12,Ovary_chr13,Ovary_chr14,
                  Ovary_chr15,Ovary_chr16,Ovary_chr17,Ovary_chr18)

#Pituitary#
Pituitary_chr1<-read.csv("/Users/baizhonghao/Downloads/liftover/Pituitary/Pituitary_chr1.csv")
Pituitary_chr2<-read.csv("/Users/baizhonghao/Downloads/liftover/Pituitary/Pituitary_chr2.csv")
Pituitary_chr3<-read.csv("/Users/baizhonghao/Downloads/liftover/Pituitary/Pituitary_chr3.csv")
Pituitary_chr4<-read.csv("/Users/baizhonghao/Downloads/liftover/Pituitary/Pituitary_chr4.csv")
Pituitary_chr5<-read.csv("/Users/baizhonghao/Downloads/liftover/Pituitary/Pituitary_chr5.csv")
Pituitary_chr6<-read.csv("/Users/baizhonghao/Downloads/liftover/Pituitary/Pituitary_chr6.csv")
Pituitary_chr7<-read.csv("/Users/baizhonghao/Downloads/liftover/Pituitary/Pituitary_chr7.csv")
Pituitary_chr8<-read.csv("/Users/baizhonghao/Downloads/liftover/Pituitary/Pituitary_chr8.csv")
Pituitary_chr9<-read.csv("/Users/baizhonghao/Downloads/liftover/Pituitary/Pituitary_chr9.csv")
Pituitary_chr10<-read.csv("/Users/baizhonghao/Downloads/liftover/Pituitary/Pituitary_chr10.csv")
Pituitary_chr11<-read.csv("/Users/baizhonghao/Downloads/liftover/Pituitary/Pituitary_chr11.csv")
Pituitary_chr12<-read.csv("/Users/baizhonghao/Downloads/liftover/Pituitary/Pituitary_chr12.csv")
Pituitary_chr13<-read.csv("/Users/baizhonghao/Downloads/liftover/Pituitary/Pituitary_chr13.csv")
Pituitary_chr14<-read.csv("/Users/baizhonghao/Downloads/liftover/Pituitary/Pituitary_chr14.csv")
Pituitary_chr15<-read.csv("/Users/baizhonghao/Downloads/liftover/Pituitary/Pituitary_chr15.csv")
Pituitary_chr16<-read.csv("/Users/baizhonghao/Downloads/liftover/Pituitary/Pituitary_chr16.csv")
Pituitary_chr17<-read.csv("/Users/baizhonghao/Downloads/liftover/Pituitary/Pituitary_chr17.csv")
Pituitary_chr18<-read.csv("/Users/baizhonghao/Downloads/liftover/Pituitary/Pituitary_chr18.csv")
Pituitary_loci<-rbind(Pituitary_chr1,Pituitary_chr2,Pituitary_chr3,Pituitary_chr4,Pituitary_chr5,Pituitary_chr6,Pituitary_chr7,
                      Pituitary_chr8,Pituitary_chr9,Pituitary_chr10,Pituitary_chr11,Pituitary_chr12,Pituitary_chr13,Pituitary_chr14,
                      Pituitary_chr15,Pituitary_chr16,Pituitary_chr17,Pituitary_chr18)

#Spleen#
Spleen_chr1<-read.csv("/Users/baizhonghao/Downloads/liftover/Spleen/Spleen_chr1.csv")
Spleen_chr2<-read.csv("/Users/baizhonghao/Downloads/liftover/Spleen/Spleen_chr2.csv")
Spleen_chr3<-read.csv("/Users/baizhonghao/Downloads/liftover/Spleen/Spleen_chr3.csv")
Spleen_chr4<-read.csv("/Users/baizhonghao/Downloads/liftover/Spleen/Spleen_chr4.csv")
Spleen_chr5<-read.csv("/Users/baizhonghao/Downloads/liftover/Spleen/Spleen_chr5.csv")
Spleen_chr6<-read.csv("/Users/baizhonghao/Downloads/liftover/Spleen/Spleen_chr6.csv")
Spleen_chr7<-read.csv("/Users/baizhonghao/Downloads/liftover/Spleen/Spleen_chr7.csv")
Spleen_chr8<-read.csv("/Users/baizhonghao/Downloads/liftover/Spleen/Spleen_chr8.csv")
Spleen_chr9<-read.csv("/Users/baizhonghao/Downloads/liftover/Spleen/Spleen_chr9.csv")
Spleen_chr10<-read.csv("/Users/baizhonghao/Downloads/liftover/Spleen/Spleen_chr10.csv")
Spleen_chr11<-read.csv("/Users/baizhonghao/Downloads/liftover/Spleen/Spleen_chr11.csv")
Spleen_chr12<-read.csv("/Users/baizhonghao/Downloads/liftover/Spleen/Spleen_chr12.csv")
Spleen_chr13<-read.csv("/Users/baizhonghao/Downloads/liftover/Spleen/Spleen_chr13.csv")
Spleen_chr14<-read.csv("/Users/baizhonghao/Downloads/liftover/Spleen/Spleen_chr14.csv")
Spleen_chr15<-read.csv("/Users/baizhonghao/Downloads/liftover/Spleen/Spleen_chr15.csv")
Spleen_chr16<-read.csv("/Users/baizhonghao/Downloads/liftover/Spleen/Spleen_chr16.csv")
Spleen_chr17<-read.csv("/Users/baizhonghao/Downloads/liftover/Spleen/Spleen_chr17.csv")
Spleen_chr18<-read.csv("/Users/baizhonghao/Downloads/liftover/Spleen/Spleen_chr18.csv")
Spleen_loci<-rbind(Spleen_chr1,Spleen_chr2,Spleen_chr3,Spleen_chr4,Spleen_chr5,Spleen_chr6,Spleen_chr7,
                   Spleen_chr8,Spleen_chr9,Spleen_chr10,Spleen_chr11,Spleen_chr12,Spleen_chr13,Spleen_chr14,
                   Spleen_chr15,Spleen_chr16,Spleen_chr17,Spleen_chr18)

#Testis#
Testis_chr1<-read.csv("/Users/baizhonghao/Downloads/liftover/Testis/Testis_chr1.csv")
Testis_chr2<-read.csv("/Users/baizhonghao/Downloads/liftover/Testis/Testis_chr2.csv")
Testis_chr3<-read.csv("/Users/baizhonghao/Downloads/liftover/Testis/Testis_chr3.csv")
Testis_chr4<-read.csv("/Users/baizhonghao/Downloads/liftover/Testis/Testis_chr4.csv")
Testis_chr5<-read.csv("/Users/baizhonghao/Downloads/liftover/Testis/Testis_chr5.csv")
Testis_chr6<-read.csv("/Users/baizhonghao/Downloads/liftover/Testis/Testis_chr6.csv")
Testis_chr7<-read.csv("/Users/baizhonghao/Downloads/liftover/Testis/Testis_chr7.csv")
Testis_chr8<-read.csv("/Users/baizhonghao/Downloads/liftover/Testis/Testis_chr8.csv")
Testis_chr9<-read.csv("/Users/baizhonghao/Downloads/liftover/Testis/Testis_chr9.csv")
Testis_chr10<-read.csv("/Users/baizhonghao/Downloads/liftover/Testis/Testis_chr10.csv")
Testis_chr11<-read.csv("/Users/baizhonghao/Downloads/liftover/Testis/Testis_chr11.csv")
Testis_chr12<-read.csv("/Users/baizhonghao/Downloads/liftover/Testis/Testis_chr12.csv")
Testis_chr13<-read.csv("/Users/baizhonghao/Downloads/liftover/Testis/Testis_chr13.csv")
Testis_chr14<-read.csv("/Users/baizhonghao/Downloads/liftover/Testis/Testis_chr14.csv")
Testis_chr15<-read.csv("/Users/baizhonghao/Downloads/liftover/Testis/Testis_chr15.csv")
Testis_chr16<-read.csv("/Users/baizhonghao/Downloads/liftover/Testis/Testis_chr16.csv")
Testis_chr17<-read.csv("/Users/baizhonghao/Downloads/liftover/Testis/Testis_chr17.csv")
Testis_chr18<-read.csv("/Users/baizhonghao/Downloads/liftover/Testis/Testis_chr18.csv")
Testis_loci<-rbind(Testis_chr1,Testis_chr2,Testis_chr3,Testis_chr4,Testis_chr5,Testis_chr6,Testis_chr7,
                   Testis_chr8,Testis_chr9,Testis_chr10,Testis_chr11,Testis_chr12,Testis_chr13,Testis_chr14,
                   Testis_chr15,Testis_chr16,Testis_chr17,Testis_chr18)

#Uterus#
Uterus_chr1<-read.csv("/Users/baizhonghao/Downloads/liftover/Uterus/Uterus_chr1.csv")
Uterus_chr2<-read.csv("/Users/baizhonghao/Downloads/liftover/Uterus/Uterus_chr2.csv")
Uterus_chr3<-read.csv("/Users/baizhonghao/Downloads/liftover/Uterus/Uterus_chr3.csv")
Uterus_chr4<-read.csv("/Users/baizhonghao/Downloads/liftover/Uterus/Uterus_chr4.csv")
Uterus_chr5<-read.csv("/Users/baizhonghao/Downloads/liftover/Uterus/Uterus_chr5.csv")
Uterus_chr6<-read.csv("/Users/baizhonghao/Downloads/liftover/Uterus/Uterus_chr6.csv")
Uterus_chr7<-read.csv("/Users/baizhonghao/Downloads/liftover/Uterus/Uterus_chr7.csv")
Uterus_chr8<-read.csv("/Users/baizhonghao/Downloads/liftover/Uterus/Uterus_chr8.csv")
Uterus_chr9<-read.csv("/Users/baizhonghao/Downloads/liftover/Uterus/Uterus_chr9.csv")
Uterus_chr10<-read.csv("/Users/baizhonghao/Downloads/liftover/Uterus/Uterus_chr10.csv")
Uterus_chr11<-read.csv("/Users/baizhonghao/Downloads/liftover/Uterus/Uterus_chr11.csv")
Uterus_chr12<-read.csv("/Users/baizhonghao/Downloads/liftover/Uterus/Uterus_chr12.csv")
Uterus_chr13<-read.csv("/Users/baizhonghao/Downloads/liftover/Uterus/Uterus_chr13.csv")
Uterus_chr14<-read.csv("/Users/baizhonghao/Downloads/liftover/Uterus/Uterus_chr14.csv")
Uterus_chr15<-read.csv("/Users/baizhonghao/Downloads/liftover/Uterus/Uterus_chr15.csv")
Uterus_chr16<-read.csv("/Users/baizhonghao/Downloads/liftover/Uterus/Uterus_chr16.csv")
Uterus_chr17<-read.csv("/Users/baizhonghao/Downloads/liftover/Uterus/Uterus_chr17.csv")
Uterus_chr18<-read.csv("/Users/baizhonghao/Downloads/liftover/Uterus/Uterus_chr18.csv")
Uterus_loci<-rbind(Uterus_chr1,Uterus_chr2,Uterus_chr3,Uterus_chr4,Uterus_chr5,Uterus_chr6,Uterus_chr7,
                   Uterus_chr8,Uterus_chr9,Uterus_chr10,Uterus_chr11,Uterus_chr12,Uterus_chr13,Uterus_chr14,
                   Uterus_chr15,Uterus_chr16,Uterus_chr17,Uterus_chr18)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose$variant_id,split="_"))))
eqtl_Adipose$chr<-tmp$V1
eqtl_Adipose$loci<-tmp$V2

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery$variant_id,split="_"))))
eqtl_Artery$chr<-tmp$V1
eqtl_Artery$loci<-tmp$V2

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood$variant_id,split="_"))))
eqtl_Blood$chr<-tmp$V1
eqtl_Blood$loci<-tmp$V2

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon$variant_id,split="_"))))
eqtl_Colon$chr<-tmp$V1
eqtl_Colon$loci<-tmp$V2

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex$variant_id,split="_"))))
eqtl_Frontal_cortex$chr<-tmp$V1
eqtl_Frontal_cortex$loci<-tmp$V2

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart$variant_id,split="_"))))
eqtl_Heart$chr<-tmp$V1
eqtl_Heart$loci<-tmp$V2

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus$variant_id,split="_"))))
eqtl_Hypothalamus$chr<-tmp$V1
eqtl_Hypothalamus$loci<-tmp$V2

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum$variant_id,split="_"))))
eqtl_Ileum$chr<-tmp$V1
eqtl_Ileum$loci<-tmp$V2

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney$variant_id,split="_"))))
eqtl_Kidney$chr<-tmp$V1
eqtl_Kidney$loci<-tmp$V2

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver$variant_id,split="_"))))
eqtl_Liver$chr<-tmp$V1
eqtl_Liver$loci<-tmp$V2

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung$variant_id,split="_"))))
eqtl_Lung$chr<-tmp$V1
eqtl_Lung$loci<-tmp$V2

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_loci$variant_id,split="_"))))
eqtl_Muscle_loci$chr<-tmp$V1
eqtl_Muscle_loci$loci<-tmp$V2

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary$variant_id,split="_"))))
eqtl_Ovary$chr<-tmp$V1
eqtl_Ovary$loci<-tmp$V2

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary$variant_id,split="_"))))
eqtl_Pituitary$chr<-tmp$V1
eqtl_Pituitary$loci<-tmp$V2

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen$variant_id,split="_"))))
eqtl_Spleen$chr<-tmp$V1
eqtl_Spleen$loci<-tmp$V2

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis$variant_id,split="_"))))
eqtl_Testis$chr<-tmp$V1
eqtl_Testis$loci<-tmp$V2

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus$variant_id,split="_"))))
eqtl_Uterus$chr<-tmp$V1
eqtl_Uterus$loci<-tmp$V2

eqtl_Adipose$index<-paste0(eqtl_Adipose$chr,"-",eqtl_Adipose$loci)
Adipose_loci$index<-paste0(Adipose_loci$chr,"-",Adipose_loci$pig_pos)
same_Adipose<-intersect(Adipose_loci$index,eqtl_Adipose$index)

eqtl_Artery$index<-paste0(eqtl_Artery$chr,"-",eqtl_Artery$loci)
Artery_loci$index<-paste0(Artery_loci$chr,"-",Artery_loci$pig_pos)
same_Artery<-intersect(Artery_loci$index,eqtl_Artery$index)

eqtl_Blood$index<-paste0(eqtl_Blood$chr,"-",eqtl_Blood$loci)
Blood_loci$index<-paste0(Blood_loci$chr,"-",Blood_loci$pig_pos)
same_Blood<-intersect(Blood_loci$index,eqtl_Blood$index)

eqtl_Colon$index<-paste0(eqtl_Colon$chr,"-",eqtl_Colon$loci)
Colon_loci$index<-paste0(Colon_loci$chr,"-",Colon_loci$pig_pos)
same_Colon<-intersect(Colon_loci$index,eqtl_Colon$index)

eqtl_Frontal_cortex$index<-paste0(eqtl_Frontal_cortex$chr,"-",eqtl_Frontal_cortex$loci)
Frontal_cortex_loci$index<-paste0(Frontal_cortex_loci$chr,"-",Frontal_cortex_loci$pig_pos)
same_Frontal_cortex<-intersect(Frontal_cortex_loci$index,eqtl_Frontal_cortex$index)

eqtl_Heart$index<-paste0(eqtl_Heart$chr,"-",eqtl_Heart$loci)
Heart_loci$index<-paste0(Heart_loci$chr,"-",Heart_loci$pig_pos)
same_Heart<-intersect(Heart_loci$index,eqtl_Heart$index)

eqtl_Hypothalamus$index<-paste0(eqtl_Hypothalamus$chr,"-",eqtl_Hypothalamus$loci)
Hypothalamus_loci$index<-paste0(Hypothalamus_loci$chr,"-",Hypothalamus_loci$pig_pos)
same_Hypothalamus<-intersect(Hypothalamus_loci$index,eqtl_Hypothalamus$index)

eqtl_Ileum$index<-paste0(eqtl_Ileum$chr,"-",eqtl_Ileum$loci)
Ileum_loci$index<-paste0(Ileum_loci$chr,"-",Ileum_loci$pig_pos)
same_Ileum<-intersect(Ileum_loci$index,eqtl_Ileum$index)

eqtl_Kidney$index<-paste0(eqtl_Kidney$chr,"-",eqtl_Kidney$loci)
Kidney_loci$index<-paste0(Kidney_loci$chr,"-",Kidney_loci$pig_pos)
same_Kidney<-intersect(Kidney_loci$index,eqtl_Kidney$index)

eqtl_Liver$index<-paste0(eqtl_Liver$chr,"-",eqtl_Liver$loci)
Liver_loci$index<-paste0(Liver_loci$chr,"-",Liver_loci$pig_pos)
same_Liver<-intersect(Liver_loci$index,eqtl_Liver$index)

eqtl_Lung$index<-paste0(eqtl_Lung$chr,"-",eqtl_Lung$loci)
Lung_loci$index<-paste0(Lung_loci$chr,"-",Lung_loci$pig_pos)
same_Lung<-intersect(Lung_loci$index,eqtl_Lung$index)

eqtl_Muscle$index<-paste0(eqtl_Muscle$chr,"-",eqtl_Muscle$loci)
Muscle_loci$index<-paste0(Muscle_loci$chr,"-",Muscle_loci$pig_pos)
same_Muscle<-intersect(Muscle_loci$index,eqtl_Muscle$index)

eqtl_Ovary$index<-paste0(eqtl_Ovary$chr,"-",eqtl_Ovary$loci)
Ovary_loci$index<-paste0(Ovary_loci$chr,"-",Ovary_loci$pig_pos)
same_Ovary<-intersect(Ovary_loci$index,eqtl_Ovary$index)

eqtl_Pituitary$index<-paste0(eqtl_Pituitary$chr,"-",eqtl_Pituitary$loci)
Pituitary_loci$index<-paste0(Pituitary_loci$chr,"-",Pituitary_loci$pig_pos)
same_Pituitary<-intersect(Pituitary_loci$index,eqtl_Pituitary$index)

eqtl_Spleen$index<-paste0(eqtl_Spleen$chr,"-",eqtl_Spleen$loci)
Spleen_loci$index<-paste0(Spleen_loci$chr,"-",Spleen_loci$pig_pos)
same_Spleen<-intersect(Spleen_loci$index,eqtl_Spleen$index)

eqtl_Testis$index<-paste0(eqtl_Testis$chr,"-",eqtl_Testis$loci)
Testis_loci$index<-paste0(Testis_loci$chr,"-",Testis_loci$pig_pos)
same_Testis<-intersect(Testis_loci$index,eqtl_Testis$index)

eqtl_Uterus$index<-paste0(eqtl_Uterus$chr,"-",eqtl_Uterus$loci)
Uterus_loci$index<-paste0(Uterus_loci$chr,"-",Uterus_loci$pig_pos)
same_Uterus<-intersect(Uterus_loci$index,eqtl_Uterus$index)

###############################################################################
#sequencing conservation -- PhastCons#
library(data.table)
eqtlfm<-as.data.frame(fread("/Users/baizhonghao/Downloads/Pig.chrAll.gene_coverage_phastCons.txt"))
orthologous_anno<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/human-pig.txt"))
eqtlfm$PhastCons[which(eqtlfm$PhastCons=='.')]<-0
eqtlfm$PhastCons<-as.numeric(eqtlfm$PhastCons)
length(eqtlfm$gene_id)
eqtlfm<-eqtlfm[order(eqtlfm$PhastCons,decreasing = T),]
top10<-eqtlfm[1:2900,]
top10_20<-eqtlfm[2901:5800,]
top20_30<-eqtlfm[5801:8700,]
top30_40<-eqtlfm[8701:11600,]
top40_50<-eqtlfm[11601:14500,]
top50_60<-eqtlfm[14501:17400,]
top60_70<-eqtlfm[17401:20300,]
top70_80<-eqtlfm[20301:23200,]
top80_90<-eqtlfm[23201:26100,]
last10<-eqtlfm[26101:29000,]

#Adipose#
top10_Adipose_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Adipose)
top10_Adipose_Pig<-intersect(intersect(top10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Adipose)
top10_Adipose_overlap<-intersect(annotation$`Gene stable ID`[match(top10_Adipose_Pig,annotation$`Pig gene stable ID`)],top10_Adipose_Human)

top10_20_Adipose_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Adipose)
top10_20_Adipose_Pig<-intersect(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),pig_egene_Adipose)
top10_20_Adipose_overlap<-intersect(annotation$`Gene stable ID`[match(top10_20_Adipose_Pig,annotation$`Pig gene stable ID`)],top10_20_Adipose_Human)

top20_30_Adipose_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Adipose)
top20_30_Adipose_Pig<-intersect(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),pig_egene_Adipose)
top20_30_Adipose_overlap<-intersect(annotation$`Gene stable ID`[match(top20_30_Adipose_Pig,annotation$`Pig gene stable ID`)],top20_30_Adipose_Human)

top30_40_Adipose_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Adipose)
top30_40_Adipose_Pig<-intersect(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),pig_egene_Adipose)
top30_40_Adipose_overlap<-intersect(annotation$`Gene stable ID`[match(top30_40_Adipose_Pig,annotation$`Pig gene stable ID`)],top30_40_Adipose_Human)

top40_50_Adipose_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Adipose)
top40_50_Adipose_Pig<-intersect(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),pig_egene_Adipose)
top40_50_Adipose_overlap<-intersect(annotation$`Gene stable ID`[match(top40_50_Adipose_Pig,annotation$`Pig gene stable ID`)],top40_50_Adipose_Human)

top50_60_Adipose_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Adipose)
top50_60_Adipose_Pig<-intersect(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),pig_egene_Adipose)
top50_60_Adipose_overlap<-intersect(annotation$`Gene stable ID`[match(top50_60_Adipose_Pig,annotation$`Pig gene stable ID`)],top50_60_Adipose_Human)

top60_70_Adipose_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Adipose)
top60_70_Adipose_Pig<-intersect(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),pig_egene_Adipose)
top60_70_Adipose_overlap<-intersect(annotation$`Gene stable ID`[match(top60_70_Adipose_Pig,annotation$`Pig gene stable ID`)],top60_70_Adipose_Human)

top70_80_Adipose_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Adipose)
top70_80_Adipose_Pig<-intersect(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),pig_egene_Adipose)
top70_80_Adipose_overlap<-intersect(annotation$`Gene stable ID`[match(top70_80_Adipose_Pig,annotation$`Pig gene stable ID`)],top70_80_Adipose_Human)

top80_90_Adipose_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Adipose)
top80_90_Adipose_Pig<-intersect(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),pig_egene_Adipose)
top80_90_Adipose_overlap<-intersect(annotation$`Gene stable ID`[match(top80_90_Adipose_Pig,annotation$`Pig gene stable ID`)],top80_90_Adipose_Human)

last10_Adipose_Human<-intersect(annotation$`Gene stable ID`[match(intersect(last10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Adipose)
last10_Adipose_Pig<-intersect(intersect(last10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Adipose)
last10_Adipose_overlap<-intersect(annotation$`Gene stable ID`[match(last10_Adipose_Pig,annotation$`Pig gene stable ID`)],last10_Adipose_Human)

#Artery#
top10_Artery_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Artery)
top10_Artery_Pig<-intersect(intersect(top10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Artery)
top10_Artery_overlap<-intersect(annotation$`Gene stable ID`[match(top10_Artery_Pig,annotation$`Pig gene stable ID`)],top10_Artery_Human)

top10_20_Artery_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Artery)
top10_20_Artery_Pig<-intersect(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),pig_egene_Artery)
top10_20_Artery_overlap<-intersect(annotation$`Gene stable ID`[match(top10_20_Artery_Pig,annotation$`Pig gene stable ID`)],top10_20_Artery_Human)

top20_30_Artery_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Artery)
top20_30_Artery_Pig<-intersect(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),pig_egene_Artery)
top20_30_Artery_overlap<-intersect(annotation$`Gene stable ID`[match(top20_30_Artery_Pig,annotation$`Pig gene stable ID`)],top20_30_Artery_Human)

top30_40_Artery_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Artery)
top30_40_Artery_Pig<-intersect(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),pig_egene_Artery)
top30_40_Artery_overlap<-intersect(annotation$`Gene stable ID`[match(top30_40_Artery_Pig,annotation$`Pig gene stable ID`)],top30_40_Artery_Human)

top40_50_Artery_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Artery)
top40_50_Artery_Pig<-intersect(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),pig_egene_Artery)
top40_50_Artery_overlap<-intersect(annotation$`Gene stable ID`[match(top40_50_Artery_Pig,annotation$`Pig gene stable ID`)],top40_50_Artery_Human)

top50_60_Artery_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Artery)
top50_60_Artery_Pig<-intersect(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),pig_egene_Artery)
top50_60_Artery_overlap<-intersect(annotation$`Gene stable ID`[match(top50_60_Artery_Pig,annotation$`Pig gene stable ID`)],top50_60_Artery_Human)

top60_70_Artery_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Artery)
top60_70_Artery_Pig<-intersect(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),pig_egene_Artery)
top60_70_Artery_overlap<-intersect(annotation$`Gene stable ID`[match(top60_70_Artery_Pig,annotation$`Pig gene stable ID`)],top60_70_Artery_Human)

top70_80_Artery_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Artery)
top70_80_Artery_Pig<-intersect(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),pig_egene_Artery)
top70_80_Artery_overlap<-intersect(annotation$`Gene stable ID`[match(top70_80_Artery_Pig,annotation$`Pig gene stable ID`)],top70_80_Artery_Human)

top80_90_Artery_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Artery)
top80_90_Artery_Pig<-intersect(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),pig_egene_Artery)
top80_90_Artery_overlap<-intersect(annotation$`Gene stable ID`[match(top80_90_Artery_Pig,annotation$`Pig gene stable ID`)],top80_90_Artery_Human)

last10_Artery_Human<-intersect(annotation$`Gene stable ID`[match(intersect(last10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Artery)
last10_Artery_Pig<-intersect(intersect(last10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Artery)
last10_Artery_overlap<-intersect(annotation$`Gene stable ID`[match(last10_Artery_Pig,annotation$`Pig gene stable ID`)],last10_Artery_Human)

#Blood#
top10_Blood_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Blood)
top10_Blood_Pig<-intersect(intersect(top10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Blood)
top10_Blood_overlap<-intersect(annotation$`Gene stable ID`[match(top10_Blood_Pig,annotation$`Pig gene stable ID`)],top10_Blood_Human)

top10_20_Blood_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Blood)
top10_20_Blood_Pig<-intersect(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),pig_egene_Blood)
top10_20_Blood_overlap<-intersect(annotation$`Gene stable ID`[match(top10_20_Blood_Pig,annotation$`Pig gene stable ID`)],top10_20_Blood_Human)

top20_30_Blood_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Blood)
top20_30_Blood_Pig<-intersect(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),pig_egene_Blood)
top20_30_Blood_overlap<-intersect(annotation$`Gene stable ID`[match(top20_30_Blood_Pig,annotation$`Pig gene stable ID`)],top20_30_Blood_Human)

top30_40_Blood_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Blood)
top30_40_Blood_Pig<-intersect(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),pig_egene_Blood)
top30_40_Blood_overlap<-intersect(annotation$`Gene stable ID`[match(top30_40_Blood_Pig,annotation$`Pig gene stable ID`)],top30_40_Blood_Human)

top40_50_Blood_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Blood)
top40_50_Blood_Pig<-intersect(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),pig_egene_Blood)
top40_50_Blood_overlap<-intersect(annotation$`Gene stable ID`[match(top40_50_Blood_Pig,annotation$`Pig gene stable ID`)],top40_50_Blood_Human)

top50_60_Blood_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Blood)
top50_60_Blood_Pig<-intersect(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),pig_egene_Blood)
top50_60_Blood_overlap<-intersect(annotation$`Gene stable ID`[match(top50_60_Blood_Pig,annotation$`Pig gene stable ID`)],top50_60_Blood_Human)

top60_70_Blood_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Blood)
top60_70_Blood_Pig<-intersect(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),pig_egene_Blood)
top60_70_Blood_overlap<-intersect(annotation$`Gene stable ID`[match(top60_70_Blood_Pig,annotation$`Pig gene stable ID`)],top60_70_Blood_Human)

top70_80_Blood_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Blood)
top70_80_Blood_Pig<-intersect(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),pig_egene_Blood)
top70_80_Blood_overlap<-intersect(annotation$`Gene stable ID`[match(top70_80_Blood_Pig,annotation$`Pig gene stable ID`)],top70_80_Blood_Human)

top80_90_Blood_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Blood)
top80_90_Blood_Pig<-intersect(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),pig_egene_Blood)
top80_90_Blood_overlap<-intersect(annotation$`Gene stable ID`[match(top80_90_Blood_Pig,annotation$`Pig gene stable ID`)],top80_90_Blood_Human)

last10_Blood_Human<-intersect(annotation$`Gene stable ID`[match(intersect(last10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Blood)
last10_Blood_Pig<-intersect(intersect(last10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Blood)
last10_Blood_overlap<-intersect(annotation$`Gene stable ID`[match(last10_Blood_Pig,annotation$`Pig gene stable ID`)],last10_Blood_Human)

#Colon#
top10_Colon_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Colon)
top10_Colon_Pig<-intersect(intersect(top10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Colon)
top10_Colon_overlap<-intersect(annotation$`Gene stable ID`[match(top10_Colon_Pig,annotation$`Pig gene stable ID`)],top10_Colon_Human)

top10_20_Colon_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Colon)
top10_20_Colon_Pig<-intersect(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),pig_egene_Colon)
top10_20_Colon_overlap<-intersect(annotation$`Gene stable ID`[match(top10_20_Colon_Pig,annotation$`Pig gene stable ID`)],top10_20_Colon_Human)

top20_30_Colon_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Colon)
top20_30_Colon_Pig<-intersect(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),pig_egene_Colon)
top20_30_Colon_overlap<-intersect(annotation$`Gene stable ID`[match(top20_30_Colon_Pig,annotation$`Pig gene stable ID`)],top20_30_Colon_Human)

top30_40_Colon_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Colon)
top30_40_Colon_Pig<-intersect(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),pig_egene_Colon)
top30_40_Colon_overlap<-intersect(annotation$`Gene stable ID`[match(top30_40_Colon_Pig,annotation$`Pig gene stable ID`)],top30_40_Colon_Human)

top40_50_Colon_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Colon)
top40_50_Colon_Pig<-intersect(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),pig_egene_Colon)
top40_50_Colon_overlap<-intersect(annotation$`Gene stable ID`[match(top40_50_Colon_Pig,annotation$`Pig gene stable ID`)],top40_50_Colon_Human)

top50_60_Colon_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Colon)
top50_60_Colon_Pig<-intersect(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),pig_egene_Colon)
top50_60_Colon_overlap<-intersect(annotation$`Gene stable ID`[match(top50_60_Colon_Pig,annotation$`Pig gene stable ID`)],top50_60_Colon_Human)

top60_70_Colon_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Colon)
top60_70_Colon_Pig<-intersect(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),pig_egene_Colon)
top60_70_Colon_overlap<-intersect(annotation$`Gene stable ID`[match(top60_70_Colon_Pig,annotation$`Pig gene stable ID`)],top60_70_Colon_Human)

top70_80_Colon_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Colon)
top70_80_Colon_Pig<-intersect(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),pig_egene_Colon)
top70_80_Colon_overlap<-intersect(annotation$`Gene stable ID`[match(top70_80_Colon_Pig,annotation$`Pig gene stable ID`)],top70_80_Colon_Human)

top80_90_Colon_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Colon)
top80_90_Colon_Pig<-intersect(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),pig_egene_Colon)
top80_90_Colon_overlap<-intersect(annotation$`Gene stable ID`[match(top80_90_Colon_Pig,annotation$`Pig gene stable ID`)],top80_90_Colon_Human)

last10_Colon_Human<-intersect(annotation$`Gene stable ID`[match(intersect(last10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Colon)
last10_Colon_Pig<-intersect(intersect(last10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Colon)
last10_Colon_overlap<-intersect(annotation$`Gene stable ID`[match(last10_Colon_Pig,annotation$`Pig gene stable ID`)],last10_Colon_Human)

#Frontal_cortex#
top10_Frontal_cortex_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Frontal_cortex)
top10_Frontal_cortex_Pig<-intersect(intersect(top10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Frontal_cortex)
top10_Frontal_cortex_overlap<-intersect(annotation$`Gene stable ID`[match(top10_Frontal_cortex_Pig,annotation$`Pig gene stable ID`)],top10_Frontal_cortex_Human)

top10_20_Frontal_cortex_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Frontal_cortex)
top10_20_Frontal_cortex_Pig<-intersect(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),pig_egene_Frontal_cortex)
top10_20_Frontal_cortex_overlap<-intersect(annotation$`Gene stable ID`[match(top10_20_Frontal_cortex_Pig,annotation$`Pig gene stable ID`)],top10_20_Frontal_cortex_Human)

top20_30_Frontal_cortex_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Frontal_cortex)
top20_30_Frontal_cortex_Pig<-intersect(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),pig_egene_Frontal_cortex)
top20_30_Frontal_cortex_overlap<-intersect(annotation$`Gene stable ID`[match(top20_30_Frontal_cortex_Pig,annotation$`Pig gene stable ID`)],top20_30_Frontal_cortex_Human)

top30_40_Frontal_cortex_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Frontal_cortex)
top30_40_Frontal_cortex_Pig<-intersect(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),pig_egene_Frontal_cortex)
top30_40_Frontal_cortex_overlap<-intersect(annotation$`Gene stable ID`[match(top30_40_Frontal_cortex_Pig,annotation$`Pig gene stable ID`)],top30_40_Frontal_cortex_Human)

top40_50_Frontal_cortex_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Frontal_cortex)
top40_50_Frontal_cortex_Pig<-intersect(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),pig_egene_Frontal_cortex)
top40_50_Frontal_cortex_overlap<-intersect(annotation$`Gene stable ID`[match(top40_50_Frontal_cortex_Pig,annotation$`Pig gene stable ID`)],top40_50_Frontal_cortex_Human)

top50_60_Frontal_cortex_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Frontal_cortex)
top50_60_Frontal_cortex_Pig<-intersect(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),pig_egene_Frontal_cortex)
top50_60_Frontal_cortex_overlap<-intersect(annotation$`Gene stable ID`[match(top50_60_Frontal_cortex_Pig,annotation$`Pig gene stable ID`)],top50_60_Frontal_cortex_Human)

top60_70_Frontal_cortex_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Frontal_cortex)
top60_70_Frontal_cortex_Pig<-intersect(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),pig_egene_Frontal_cortex)
top60_70_Frontal_cortex_overlap<-intersect(annotation$`Gene stable ID`[match(top60_70_Frontal_cortex_Pig,annotation$`Pig gene stable ID`)],top60_70_Frontal_cortex_Human)

top70_80_Frontal_cortex_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Frontal_cortex)
top70_80_Frontal_cortex_Pig<-intersect(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),pig_egene_Frontal_cortex)
top70_80_Frontal_cortex_overlap<-intersect(annotation$`Gene stable ID`[match(top70_80_Frontal_cortex_Pig,annotation$`Pig gene stable ID`)],top70_80_Frontal_cortex_Human)

top80_90_Frontal_cortex_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Frontal_cortex)
top80_90_Frontal_cortex_Pig<-intersect(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),pig_egene_Frontal_cortex)
top80_90_Frontal_cortex_overlap<-intersect(annotation$`Gene stable ID`[match(top80_90_Frontal_cortex_Pig,annotation$`Pig gene stable ID`)],top80_90_Frontal_cortex_Human)

last10_Frontal_cortex_Human<-intersect(annotation$`Gene stable ID`[match(intersect(last10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Frontal_cortex)
last10_Frontal_cortex_Pig<-intersect(intersect(last10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Frontal_cortex)
last10_Frontal_cortex_overlap<-intersect(annotation$`Gene stable ID`[match(last10_Frontal_cortex_Pig,annotation$`Pig gene stable ID`)],last10_Frontal_cortex_Human)

#Heart#
top10_Heart_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Heart)
top10_Heart_Pig<-intersect(intersect(top10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Heart)
top10_Heart_overlap<-intersect(annotation$`Gene stable ID`[match(top10_Heart_Pig,annotation$`Pig gene stable ID`)],top10_Heart_Human)

top10_20_Heart_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Heart)
top10_20_Heart_Pig<-intersect(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),pig_egene_Heart)
top10_20_Heart_overlap<-intersect(annotation$`Gene stable ID`[match(top10_20_Heart_Pig,annotation$`Pig gene stable ID`)],top10_20_Heart_Human)

top20_30_Heart_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Heart)
top20_30_Heart_Pig<-intersect(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),pig_egene_Heart)
top20_30_Heart_overlap<-intersect(annotation$`Gene stable ID`[match(top20_30_Heart_Pig,annotation$`Pig gene stable ID`)],top20_30_Heart_Human)

top30_40_Heart_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Heart)
top30_40_Heart_Pig<-intersect(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),pig_egene_Heart)
top30_40_Heart_overlap<-intersect(annotation$`Gene stable ID`[match(top30_40_Heart_Pig,annotation$`Pig gene stable ID`)],top30_40_Heart_Human)

top40_50_Heart_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Heart)
top40_50_Heart_Pig<-intersect(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),pig_egene_Heart)
top40_50_Heart_overlap<-intersect(annotation$`Gene stable ID`[match(top40_50_Heart_Pig,annotation$`Pig gene stable ID`)],top40_50_Heart_Human)

top50_60_Heart_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Heart)
top50_60_Heart_Pig<-intersect(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),pig_egene_Heart)
top50_60_Heart_overlap<-intersect(annotation$`Gene stable ID`[match(top50_60_Heart_Pig,annotation$`Pig gene stable ID`)],top50_60_Heart_Human)

top60_70_Heart_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Heart)
top60_70_Heart_Pig<-intersect(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),pig_egene_Heart)
top60_70_Heart_overlap<-intersect(annotation$`Gene stable ID`[match(top60_70_Heart_Pig,annotation$`Pig gene stable ID`)],top60_70_Heart_Human)

top70_80_Heart_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Heart)
top70_80_Heart_Pig<-intersect(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),pig_egene_Heart)
top70_80_Heart_overlap<-intersect(annotation$`Gene stable ID`[match(top70_80_Heart_Pig,annotation$`Pig gene stable ID`)],top70_80_Heart_Human)

top80_90_Heart_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Heart)
top80_90_Heart_Pig<-intersect(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),pig_egene_Heart)
top80_90_Heart_overlap<-intersect(annotation$`Gene stable ID`[match(top80_90_Heart_Pig,annotation$`Pig gene stable ID`)],top80_90_Heart_Human)

last10_Heart_Human<-intersect(annotation$`Gene stable ID`[match(intersect(last10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Heart)
last10_Heart_Pig<-intersect(intersect(last10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Heart)
last10_Heart_overlap<-intersect(annotation$`Gene stable ID`[match(last10_Heart_Pig,annotation$`Pig gene stable ID`)],last10_Heart_Human)

#Hypothalamus#
top10_Hypothalamus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Hypothalamus)
top10_Hypothalamus_Pig<-intersect(intersect(top10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Hypothalamus)
top10_Hypothalamus_overlap<-intersect(annotation$`Gene stable ID`[match(top10_Hypothalamus_Pig,annotation$`Pig gene stable ID`)],top10_Hypothalamus_Human)

top10_20_Hypothalamus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Hypothalamus)
top10_20_Hypothalamus_Pig<-intersect(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),pig_egene_Hypothalamus)
top10_20_Hypothalamus_overlap<-intersect(annotation$`Gene stable ID`[match(top10_20_Hypothalamus_Pig,annotation$`Pig gene stable ID`)],top10_20_Hypothalamus_Human)

top20_30_Hypothalamus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Hypothalamus)
top20_30_Hypothalamus_Pig<-intersect(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),pig_egene_Hypothalamus)
top20_30_Hypothalamus_overlap<-intersect(annotation$`Gene stable ID`[match(top20_30_Hypothalamus_Pig,annotation$`Pig gene stable ID`)],top20_30_Hypothalamus_Human)

top30_40_Hypothalamus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Hypothalamus)
top30_40_Hypothalamus_Pig<-intersect(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),pig_egene_Hypothalamus)
top30_40_Hypothalamus_overlap<-intersect(annotation$`Gene stable ID`[match(top30_40_Hypothalamus_Pig,annotation$`Pig gene stable ID`)],top30_40_Hypothalamus_Human)

top40_50_Hypothalamus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Hypothalamus)
top40_50_Hypothalamus_Pig<-intersect(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),pig_egene_Hypothalamus)
top40_50_Hypothalamus_overlap<-intersect(annotation$`Gene stable ID`[match(top40_50_Hypothalamus_Pig,annotation$`Pig gene stable ID`)],top40_50_Hypothalamus_Human)

top50_60_Hypothalamus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Hypothalamus)
top50_60_Hypothalamus_Pig<-intersect(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),pig_egene_Hypothalamus)
top50_60_Hypothalamus_overlap<-intersect(annotation$`Gene stable ID`[match(top50_60_Hypothalamus_Pig,annotation$`Pig gene stable ID`)],top50_60_Hypothalamus_Human)

top60_70_Hypothalamus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Hypothalamus)
top60_70_Hypothalamus_Pig<-intersect(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),pig_egene_Hypothalamus)
top60_70_Hypothalamus_overlap<-intersect(annotation$`Gene stable ID`[match(top60_70_Hypothalamus_Pig,annotation$`Pig gene stable ID`)],top60_70_Hypothalamus_Human)

top70_80_Hypothalamus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Hypothalamus)
top70_80_Hypothalamus_Pig<-intersect(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),pig_egene_Hypothalamus)
top70_80_Hypothalamus_overlap<-intersect(annotation$`Gene stable ID`[match(top70_80_Hypothalamus_Pig,annotation$`Pig gene stable ID`)],top70_80_Hypothalamus_Human)

top80_90_Hypothalamus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Hypothalamus)
top80_90_Hypothalamus_Pig<-intersect(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),pig_egene_Hypothalamus)
top80_90_Hypothalamus_overlap<-intersect(annotation$`Gene stable ID`[match(top80_90_Hypothalamus_Pig,annotation$`Pig gene stable ID`)],top80_90_Hypothalamus_Human)

last10_Hypothalamus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(last10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Hypothalamus)
last10_Hypothalamus_Pig<-intersect(intersect(last10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Hypothalamus)
last10_Hypothalamus_overlap<-intersect(annotation$`Gene stable ID`[match(last10_Hypothalamus_Pig,annotation$`Pig gene stable ID`)],last10_Hypothalamus_Human)

#Ileum#
top10_Ileum_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ileum)
top10_Ileum_Pig<-intersect(intersect(top10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ileum)
top10_Ileum_overlap<-intersect(annotation$`Gene stable ID`[match(top10_Ileum_Pig,annotation$`Pig gene stable ID`)],top10_Ileum_Human)

top10_20_Ileum_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ileum)
top10_20_Ileum_Pig<-intersect(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ileum)
top10_20_Ileum_overlap<-intersect(annotation$`Gene stable ID`[match(top10_20_Ileum_Pig,annotation$`Pig gene stable ID`)],top10_20_Ileum_Human)

top20_30_Ileum_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ileum)
top20_30_Ileum_Pig<-intersect(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ileum)
top20_30_Ileum_overlap<-intersect(annotation$`Gene stable ID`[match(top20_30_Ileum_Pig,annotation$`Pig gene stable ID`)],top20_30_Ileum_Human)

top30_40_Ileum_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ileum)
top30_40_Ileum_Pig<-intersect(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ileum)
top30_40_Ileum_overlap<-intersect(annotation$`Gene stable ID`[match(top30_40_Ileum_Pig,annotation$`Pig gene stable ID`)],top30_40_Ileum_Human)

top40_50_Ileum_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ileum)
top40_50_Ileum_Pig<-intersect(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ileum)
top40_50_Ileum_overlap<-intersect(annotation$`Gene stable ID`[match(top40_50_Ileum_Pig,annotation$`Pig gene stable ID`)],top40_50_Ileum_Human)

top50_60_Ileum_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ileum)
top50_60_Ileum_Pig<-intersect(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ileum)
top50_60_Ileum_overlap<-intersect(annotation$`Gene stable ID`[match(top50_60_Ileum_Pig,annotation$`Pig gene stable ID`)],top50_60_Ileum_Human)

top60_70_Ileum_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ileum)
top60_70_Ileum_Pig<-intersect(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ileum)
top60_70_Ileum_overlap<-intersect(annotation$`Gene stable ID`[match(top60_70_Ileum_Pig,annotation$`Pig gene stable ID`)],top60_70_Ileum_Human)

top70_80_Ileum_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ileum)
top70_80_Ileum_Pig<-intersect(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ileum)
top70_80_Ileum_overlap<-intersect(annotation$`Gene stable ID`[match(top70_80_Ileum_Pig,annotation$`Pig gene stable ID`)],top70_80_Ileum_Human)

top80_90_Ileum_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ileum)
top80_90_Ileum_Pig<-intersect(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ileum)
top80_90_Ileum_overlap<-intersect(annotation$`Gene stable ID`[match(top80_90_Ileum_Pig,annotation$`Pig gene stable ID`)],top80_90_Ileum_Human)

last10_Ileum_Human<-intersect(annotation$`Gene stable ID`[match(intersect(last10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ileum)
last10_Ileum_Pig<-intersect(intersect(last10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ileum)
last10_Ileum_overlap<-intersect(annotation$`Gene stable ID`[match(last10_Ileum_Pig,annotation$`Pig gene stable ID`)],last10_Ileum_Human)

#Kidney#
top10_Kidney_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Kidney)
top10_Kidney_Pig<-intersect(intersect(top10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Kidney)
top10_Kidney_overlap<-intersect(annotation$`Gene stable ID`[match(top10_Kidney_Pig,annotation$`Pig gene stable ID`)],top10_Kidney_Human)

top10_20_Kidney_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Kidney)
top10_20_Kidney_Pig<-intersect(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),pig_egene_Kidney)
top10_20_Kidney_overlap<-intersect(annotation$`Gene stable ID`[match(top10_20_Kidney_Pig,annotation$`Pig gene stable ID`)],top10_20_Kidney_Human)

top20_30_Kidney_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Kidney)
top20_30_Kidney_Pig<-intersect(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),pig_egene_Kidney)
top20_30_Kidney_overlap<-intersect(annotation$`Gene stable ID`[match(top20_30_Kidney_Pig,annotation$`Pig gene stable ID`)],top20_30_Kidney_Human)

top30_40_Kidney_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Kidney)
top30_40_Kidney_Pig<-intersect(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),pig_egene_Kidney)
top30_40_Kidney_overlap<-intersect(annotation$`Gene stable ID`[match(top30_40_Kidney_Pig,annotation$`Pig gene stable ID`)],top30_40_Kidney_Human)

top40_50_Kidney_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Kidney)
top40_50_Kidney_Pig<-intersect(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),pig_egene_Kidney)
top40_50_Kidney_overlap<-intersect(annotation$`Gene stable ID`[match(top40_50_Kidney_Pig,annotation$`Pig gene stable ID`)],top40_50_Kidney_Human)

top50_60_Kidney_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Kidney)
top50_60_Kidney_Pig<-intersect(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),pig_egene_Kidney)
top50_60_Kidney_overlap<-intersect(annotation$`Gene stable ID`[match(top50_60_Kidney_Pig,annotation$`Pig gene stable ID`)],top50_60_Kidney_Human)

top60_70_Kidney_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Kidney)
top60_70_Kidney_Pig<-intersect(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),pig_egene_Kidney)
top60_70_Kidney_overlap<-intersect(annotation$`Gene stable ID`[match(top60_70_Kidney_Pig,annotation$`Pig gene stable ID`)],top60_70_Kidney_Human)

top70_80_Kidney_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Kidney)
top70_80_Kidney_Pig<-intersect(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),pig_egene_Kidney)
top70_80_Kidney_overlap<-intersect(annotation$`Gene stable ID`[match(top70_80_Kidney_Pig,annotation$`Pig gene stable ID`)],top70_80_Kidney_Human)

top80_90_Kidney_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Kidney)
top80_90_Kidney_Pig<-intersect(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),pig_egene_Kidney)
top80_90_Kidney_overlap<-intersect(annotation$`Gene stable ID`[match(top80_90_Kidney_Pig,annotation$`Pig gene stable ID`)],top80_90_Kidney_Human)

last10_Kidney_Human<-intersect(annotation$`Gene stable ID`[match(intersect(last10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Kidney)
last10_Kidney_Pig<-intersect(intersect(last10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Kidney)
last10_Kidney_overlap<-intersect(annotation$`Gene stable ID`[match(last10_Kidney_Pig,annotation$`Pig gene stable ID`)],last10_Kidney_Human)

#Liver#
top10_Liver_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Liver)
top10_Liver_Pig<-intersect(intersect(top10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Liver)
top10_Liver_overlap<-intersect(annotation$`Gene stable ID`[match(top10_Liver_Pig,annotation$`Pig gene stable ID`)],top10_Liver_Human)

top10_20_Liver_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Liver)
top10_20_Liver_Pig<-intersect(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),pig_egene_Liver)
top10_20_Liver_overlap<-intersect(annotation$`Gene stable ID`[match(top10_20_Liver_Pig,annotation$`Pig gene stable ID`)],top10_20_Liver_Human)

top20_30_Liver_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Liver)
top20_30_Liver_Pig<-intersect(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),pig_egene_Liver)
top20_30_Liver_overlap<-intersect(annotation$`Gene stable ID`[match(top20_30_Liver_Pig,annotation$`Pig gene stable ID`)],top20_30_Liver_Human)

top30_40_Liver_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Liver)
top30_40_Liver_Pig<-intersect(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),pig_egene_Liver)
top30_40_Liver_overlap<-intersect(annotation$`Gene stable ID`[match(top30_40_Liver_Pig,annotation$`Pig gene stable ID`)],top30_40_Liver_Human)

top40_50_Liver_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Liver)
top40_50_Liver_Pig<-intersect(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),pig_egene_Liver)
top40_50_Liver_overlap<-intersect(annotation$`Gene stable ID`[match(top40_50_Liver_Pig,annotation$`Pig gene stable ID`)],top40_50_Liver_Human)

top50_60_Liver_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Liver)
top50_60_Liver_Pig<-intersect(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),pig_egene_Liver)
top50_60_Liver_overlap<-intersect(annotation$`Gene stable ID`[match(top50_60_Liver_Pig,annotation$`Pig gene stable ID`)],top50_60_Liver_Human)

top60_70_Liver_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Liver)
top60_70_Liver_Pig<-intersect(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),pig_egene_Liver)
top60_70_Liver_overlap<-intersect(annotation$`Gene stable ID`[match(top60_70_Liver_Pig,annotation$`Pig gene stable ID`)],top60_70_Liver_Human)

top70_80_Liver_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Liver)
top70_80_Liver_Pig<-intersect(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),pig_egene_Liver)
top70_80_Liver_overlap<-intersect(annotation$`Gene stable ID`[match(top70_80_Liver_Pig,annotation$`Pig gene stable ID`)],top70_80_Liver_Human)

top80_90_Liver_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Liver)
top80_90_Liver_Pig<-intersect(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),pig_egene_Liver)
top80_90_Liver_overlap<-intersect(annotation$`Gene stable ID`[match(top80_90_Liver_Pig,annotation$`Pig gene stable ID`)],top80_90_Liver_Human)

last10_Liver_Human<-intersect(annotation$`Gene stable ID`[match(intersect(last10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Liver)
last10_Liver_Pig<-intersect(intersect(last10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Liver)
last10_Liver_overlap<-intersect(annotation$`Gene stable ID`[match(last10_Liver_Pig,annotation$`Pig gene stable ID`)],last10_Liver_Human)

#Lung#
top10_Lung_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Lung)
top10_Lung_Pig<-intersect(intersect(top10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Lung)
top10_Lung_overlap<-intersect(annotation$`Gene stable ID`[match(top10_Lung_Pig,annotation$`Pig gene stable ID`)],top10_Lung_Human)

top10_20_Lung_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Lung)
top10_20_Lung_Pig<-intersect(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),pig_egene_Lung)
top10_20_Lung_overlap<-intersect(annotation$`Gene stable ID`[match(top10_20_Lung_Pig,annotation$`Pig gene stable ID`)],top10_20_Lung_Human)

top20_30_Lung_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Lung)
top20_30_Lung_Pig<-intersect(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),pig_egene_Lung)
top20_30_Lung_overlap<-intersect(annotation$`Gene stable ID`[match(top20_30_Lung_Pig,annotation$`Pig gene stable ID`)],top20_30_Lung_Human)

top30_40_Lung_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Lung)
top30_40_Lung_Pig<-intersect(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),pig_egene_Lung)
top30_40_Lung_overlap<-intersect(annotation$`Gene stable ID`[match(top30_40_Lung_Pig,annotation$`Pig gene stable ID`)],top30_40_Lung_Human)

top40_50_Lung_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Lung)
top40_50_Lung_Pig<-intersect(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),pig_egene_Lung)
top40_50_Lung_overlap<-intersect(annotation$`Gene stable ID`[match(top40_50_Lung_Pig,annotation$`Pig gene stable ID`)],top40_50_Lung_Human)

top50_60_Lung_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Lung)
top50_60_Lung_Pig<-intersect(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),pig_egene_Lung)
top50_60_Lung_overlap<-intersect(annotation$`Gene stable ID`[match(top50_60_Lung_Pig,annotation$`Pig gene stable ID`)],top50_60_Lung_Human)

top60_70_Lung_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Lung)
top60_70_Lung_Pig<-intersect(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),pig_egene_Lung)
top60_70_Lung_overlap<-intersect(annotation$`Gene stable ID`[match(top60_70_Lung_Pig,annotation$`Pig gene stable ID`)],top60_70_Lung_Human)

top70_80_Lung_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Lung)
top70_80_Lung_Pig<-intersect(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),pig_egene_Lung)
top70_80_Lung_overlap<-intersect(annotation$`Gene stable ID`[match(top70_80_Lung_Pig,annotation$`Pig gene stable ID`)],top70_80_Lung_Human)

top80_90_Lung_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Lung)
top80_90_Lung_Pig<-intersect(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),pig_egene_Lung)
top80_90_Lung_overlap<-intersect(annotation$`Gene stable ID`[match(top80_90_Lung_Pig,annotation$`Pig gene stable ID`)],top80_90_Lung_Human)

last10_Lung_Human<-intersect(annotation$`Gene stable ID`[match(intersect(last10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Lung)
last10_Lung_Pig<-intersect(intersect(last10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Lung)
last10_Lung_overlap<-intersect(annotation$`Gene stable ID`[match(last10_Lung_Pig,annotation$`Pig gene stable ID`)],last10_Lung_Human)

#Muscle#
top10_Muscle_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Muscle)
top10_Muscle_Pig<-intersect(intersect(top10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Muscle)
top10_Muscle_overlap<-intersect(annotation$`Gene stable ID`[match(top10_Muscle_Pig,annotation$`Pig gene stable ID`)],top10_Muscle_Human)

top10_20_Muscle_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Muscle)
top10_20_Muscle_Pig<-intersect(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),pig_egene_Muscle)
top10_20_Muscle_overlap<-intersect(annotation$`Gene stable ID`[match(top10_20_Muscle_Pig,annotation$`Pig gene stable ID`)],top10_20_Muscle_Human)

top20_30_Muscle_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Muscle)
top20_30_Muscle_Pig<-intersect(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),pig_egene_Muscle)
top20_30_Muscle_overlap<-intersect(annotation$`Gene stable ID`[match(top20_30_Muscle_Pig,annotation$`Pig gene stable ID`)],top20_30_Muscle_Human)

top30_40_Muscle_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Muscle)
top30_40_Muscle_Pig<-intersect(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),pig_egene_Muscle)
top30_40_Muscle_overlap<-intersect(annotation$`Gene stable ID`[match(top30_40_Muscle_Pig,annotation$`Pig gene stable ID`)],top30_40_Muscle_Human)

top40_50_Muscle_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Muscle)
top40_50_Muscle_Pig<-intersect(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),pig_egene_Muscle)
top40_50_Muscle_overlap<-intersect(annotation$`Gene stable ID`[match(top40_50_Muscle_Pig,annotation$`Pig gene stable ID`)],top40_50_Muscle_Human)

top50_60_Muscle_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Muscle)
top50_60_Muscle_Pig<-intersect(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),pig_egene_Muscle)
top50_60_Muscle_overlap<-intersect(annotation$`Gene stable ID`[match(top50_60_Muscle_Pig,annotation$`Pig gene stable ID`)],top50_60_Muscle_Human)

top60_70_Muscle_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Muscle)
top60_70_Muscle_Pig<-intersect(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),pig_egene_Muscle)
top60_70_Muscle_overlap<-intersect(annotation$`Gene stable ID`[match(top60_70_Muscle_Pig,annotation$`Pig gene stable ID`)],top60_70_Muscle_Human)

top70_80_Muscle_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Muscle)
top70_80_Muscle_Pig<-intersect(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),pig_egene_Muscle)
top70_80_Muscle_overlap<-intersect(annotation$`Gene stable ID`[match(top70_80_Muscle_Pig,annotation$`Pig gene stable ID`)],top70_80_Muscle_Human)

top80_90_Muscle_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Muscle)
top80_90_Muscle_Pig<-intersect(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),pig_egene_Muscle)
top80_90_Muscle_overlap<-intersect(annotation$`Gene stable ID`[match(top80_90_Muscle_Pig,annotation$`Pig gene stable ID`)],top80_90_Muscle_Human)

last10_Muscle_Human<-intersect(annotation$`Gene stable ID`[match(intersect(last10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Muscle)
last10_Muscle_Pig<-intersect(intersect(last10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Muscle)
last10_Muscle_overlap<-intersect(annotation$`Gene stable ID`[match(last10_Muscle_Pig,annotation$`Pig gene stable ID`)],last10_Muscle_Human)

#Ovary#
top10_Ovary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ovary)
top10_Ovary_Pig<-intersect(intersect(top10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ovary)
top10_Ovary_overlap<-intersect(annotation$`Gene stable ID`[match(top10_Ovary_Pig,annotation$`Pig gene stable ID`)],top10_Ovary_Human)

top10_20_Ovary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ovary)
top10_20_Ovary_Pig<-intersect(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ovary)
top10_20_Ovary_overlap<-intersect(annotation$`Gene stable ID`[match(top10_20_Ovary_Pig,annotation$`Pig gene stable ID`)],top10_20_Ovary_Human)

top20_30_Ovary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ovary)
top20_30_Ovary_Pig<-intersect(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ovary)
top20_30_Ovary_overlap<-intersect(annotation$`Gene stable ID`[match(top20_30_Ovary_Pig,annotation$`Pig gene stable ID`)],top20_30_Ovary_Human)

top30_40_Ovary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ovary)
top30_40_Ovary_Pig<-intersect(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ovary)
top30_40_Ovary_overlap<-intersect(annotation$`Gene stable ID`[match(top30_40_Ovary_Pig,annotation$`Pig gene stable ID`)],top30_40_Ovary_Human)

top40_50_Ovary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ovary)
top40_50_Ovary_Pig<-intersect(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ovary)
top40_50_Ovary_overlap<-intersect(annotation$`Gene stable ID`[match(top40_50_Ovary_Pig,annotation$`Pig gene stable ID`)],top40_50_Ovary_Human)

top50_60_Ovary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ovary)
top50_60_Ovary_Pig<-intersect(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ovary)
top50_60_Ovary_overlap<-intersect(annotation$`Gene stable ID`[match(top50_60_Ovary_Pig,annotation$`Pig gene stable ID`)],top50_60_Ovary_Human)

top60_70_Ovary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ovary)
top60_70_Ovary_Pig<-intersect(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ovary)
top60_70_Ovary_overlap<-intersect(annotation$`Gene stable ID`[match(top60_70_Ovary_Pig,annotation$`Pig gene stable ID`)],top60_70_Ovary_Human)

top70_80_Ovary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ovary)
top70_80_Ovary_Pig<-intersect(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ovary)
top70_80_Ovary_overlap<-intersect(annotation$`Gene stable ID`[match(top70_80_Ovary_Pig,annotation$`Pig gene stable ID`)],top70_80_Ovary_Human)

top80_90_Ovary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ovary)
top80_90_Ovary_Pig<-intersect(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ovary)
top80_90_Ovary_overlap<-intersect(annotation$`Gene stable ID`[match(top80_90_Ovary_Pig,annotation$`Pig gene stable ID`)],top80_90_Ovary_Human)

last10_Ovary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(last10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Ovary)
last10_Ovary_Pig<-intersect(intersect(last10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Ovary)
last10_Ovary_overlap<-intersect(annotation$`Gene stable ID`[match(last10_Ovary_Pig,annotation$`Pig gene stable ID`)],last10_Ovary_Human)

#Pituitary#
top10_Pituitary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Pituitary)
top10_Pituitary_Pig<-intersect(intersect(top10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Pituitary)
top10_Pituitary_overlap<-intersect(annotation$`Gene stable ID`[match(top10_Pituitary_Pig,annotation$`Pig gene stable ID`)],top10_Pituitary_Human)

top10_20_Pituitary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Pituitary)
top10_20_Pituitary_Pig<-intersect(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),pig_egene_Pituitary)
top10_20_Pituitary_overlap<-intersect(annotation$`Gene stable ID`[match(top10_20_Pituitary_Pig,annotation$`Pig gene stable ID`)],top10_20_Pituitary_Human)

top20_30_Pituitary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Pituitary)
top20_30_Pituitary_Pig<-intersect(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),pig_egene_Pituitary)
top20_30_Pituitary_overlap<-intersect(annotation$`Gene stable ID`[match(top20_30_Pituitary_Pig,annotation$`Pig gene stable ID`)],top20_30_Pituitary_Human)

top30_40_Pituitary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Pituitary)
top30_40_Pituitary_Pig<-intersect(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),pig_egene_Pituitary)
top30_40_Pituitary_overlap<-intersect(annotation$`Gene stable ID`[match(top30_40_Pituitary_Pig,annotation$`Pig gene stable ID`)],top30_40_Pituitary_Human)

top40_50_Pituitary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Pituitary)
top40_50_Pituitary_Pig<-intersect(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),pig_egene_Pituitary)
top40_50_Pituitary_overlap<-intersect(annotation$`Gene stable ID`[match(top40_50_Pituitary_Pig,annotation$`Pig gene stable ID`)],top40_50_Pituitary_Human)

top50_60_Pituitary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Pituitary)
top50_60_Pituitary_Pig<-intersect(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),pig_egene_Pituitary)
top50_60_Pituitary_overlap<-intersect(annotation$`Gene stable ID`[match(top50_60_Pituitary_Pig,annotation$`Pig gene stable ID`)],top50_60_Pituitary_Human)

top60_70_Pituitary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Pituitary)
top60_70_Pituitary_Pig<-intersect(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),pig_egene_Pituitary)
top60_70_Pituitary_overlap<-intersect(annotation$`Gene stable ID`[match(top60_70_Pituitary_Pig,annotation$`Pig gene stable ID`)],top60_70_Pituitary_Human)

top70_80_Pituitary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Pituitary)
top70_80_Pituitary_Pig<-intersect(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),pig_egene_Pituitary)
top70_80_Pituitary_overlap<-intersect(annotation$`Gene stable ID`[match(top70_80_Pituitary_Pig,annotation$`Pig gene stable ID`)],top70_80_Pituitary_Human)

top80_90_Pituitary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Pituitary)
top80_90_Pituitary_Pig<-intersect(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),pig_egene_Pituitary)
top80_90_Pituitary_overlap<-intersect(annotation$`Gene stable ID`[match(top80_90_Pituitary_Pig,annotation$`Pig gene stable ID`)],top80_90_Pituitary_Human)

last10_Pituitary_Human<-intersect(annotation$`Gene stable ID`[match(intersect(last10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Pituitary)
last10_Pituitary_Pig<-intersect(intersect(last10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Pituitary)
last10_Pituitary_overlap<-intersect(annotation$`Gene stable ID`[match(last10_Pituitary_Pig,annotation$`Pig gene stable ID`)],last10_Pituitary_Human)

#Spleen#
top10_Spleen_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Spleen)
top10_Spleen_Pig<-intersect(intersect(top10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Spleen)
top10_Spleen_overlap<-intersect(annotation$`Gene stable ID`[match(top10_Spleen_Pig,annotation$`Pig gene stable ID`)],top10_Spleen_Human)

top10_20_Spleen_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Spleen)
top10_20_Spleen_Pig<-intersect(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),pig_egene_Spleen)
top10_20_Spleen_overlap<-intersect(annotation$`Gene stable ID`[match(top10_20_Spleen_Pig,annotation$`Pig gene stable ID`)],top10_20_Spleen_Human)

top20_30_Spleen_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Spleen)
top20_30_Spleen_Pig<-intersect(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),pig_egene_Spleen)
top20_30_Spleen_overlap<-intersect(annotation$`Gene stable ID`[match(top20_30_Spleen_Pig,annotation$`Pig gene stable ID`)],top20_30_Spleen_Human)

top30_40_Spleen_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Spleen)
top30_40_Spleen_Pig<-intersect(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),pig_egene_Spleen)
top30_40_Spleen_overlap<-intersect(annotation$`Gene stable ID`[match(top30_40_Spleen_Pig,annotation$`Pig gene stable ID`)],top30_40_Spleen_Human)

top40_50_Spleen_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Spleen)
top40_50_Spleen_Pig<-intersect(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),pig_egene_Spleen)
top40_50_Spleen_overlap<-intersect(annotation$`Gene stable ID`[match(top40_50_Spleen_Pig,annotation$`Pig gene stable ID`)],top40_50_Spleen_Human)

top50_60_Spleen_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Spleen)
top50_60_Spleen_Pig<-intersect(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),pig_egene_Spleen)
top50_60_Spleen_overlap<-intersect(annotation$`Gene stable ID`[match(top50_60_Spleen_Pig,annotation$`Pig gene stable ID`)],top50_60_Spleen_Human)

top60_70_Spleen_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Spleen)
top60_70_Spleen_Pig<-intersect(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),pig_egene_Spleen)
top60_70_Spleen_overlap<-intersect(annotation$`Gene stable ID`[match(top60_70_Spleen_Pig,annotation$`Pig gene stable ID`)],top60_70_Spleen_Human)

top70_80_Spleen_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Spleen)
top70_80_Spleen_Pig<-intersect(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),pig_egene_Spleen)
top70_80_Spleen_overlap<-intersect(annotation$`Gene stable ID`[match(top70_80_Spleen_Pig,annotation$`Pig gene stable ID`)],top70_80_Spleen_Human)

top80_90_Spleen_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Spleen)
top80_90_Spleen_Pig<-intersect(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),pig_egene_Spleen)
top80_90_Spleen_overlap<-intersect(annotation$`Gene stable ID`[match(top80_90_Spleen_Pig,annotation$`Pig gene stable ID`)],top80_90_Spleen_Human)

last10_Spleen_Human<-intersect(annotation$`Gene stable ID`[match(intersect(last10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Spleen)
last10_Spleen_Pig<-intersect(intersect(last10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Spleen)
last10_Spleen_overlap<-intersect(annotation$`Gene stable ID`[match(last10_Spleen_Pig,annotation$`Pig gene stable ID`)],last10_Spleen_Human)

#Testis#
top10_Testis_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Testis)
top10_Testis_Pig<-intersect(intersect(top10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Testis)
top10_Testis_overlap<-intersect(annotation$`Gene stable ID`[match(top10_Testis_Pig,annotation$`Pig gene stable ID`)],top10_Testis_Human)

top10_20_Testis_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Testis)
top10_20_Testis_Pig<-intersect(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),pig_egene_Testis)
top10_20_Testis_overlap<-intersect(annotation$`Gene stable ID`[match(top10_20_Testis_Pig,annotation$`Pig gene stable ID`)],top10_20_Testis_Human)

top20_30_Testis_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Testis)
top20_30_Testis_Pig<-intersect(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),pig_egene_Testis)
top20_30_Testis_overlap<-intersect(annotation$`Gene stable ID`[match(top20_30_Testis_Pig,annotation$`Pig gene stable ID`)],top20_30_Testis_Human)

top30_40_Testis_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Testis)
top30_40_Testis_Pig<-intersect(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),pig_egene_Testis)
top30_40_Testis_overlap<-intersect(annotation$`Gene stable ID`[match(top30_40_Testis_Pig,annotation$`Pig gene stable ID`)],top30_40_Testis_Human)

top40_50_Testis_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Testis)
top40_50_Testis_Pig<-intersect(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),pig_egene_Testis)
top40_50_Testis_overlap<-intersect(annotation$`Gene stable ID`[match(top40_50_Testis_Pig,annotation$`Pig gene stable ID`)],top40_50_Testis_Human)

top50_60_Testis_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Testis)
top50_60_Testis_Pig<-intersect(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),pig_egene_Testis)
top50_60_Testis_overlap<-intersect(annotation$`Gene stable ID`[match(top50_60_Testis_Pig,annotation$`Pig gene stable ID`)],top50_60_Testis_Human)

top60_70_Testis_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Testis)
top60_70_Testis_Pig<-intersect(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),pig_egene_Testis)
top60_70_Testis_overlap<-intersect(annotation$`Gene stable ID`[match(top60_70_Testis_Pig,annotation$`Pig gene stable ID`)],top60_70_Testis_Human)

top70_80_Testis_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Testis)
top70_80_Testis_Pig<-intersect(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),pig_egene_Testis)
top70_80_Testis_overlap<-intersect(annotation$`Gene stable ID`[match(top70_80_Testis_Pig,annotation$`Pig gene stable ID`)],top70_80_Testis_Human)

top80_90_Testis_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Testis)
top80_90_Testis_Pig<-intersect(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),pig_egene_Testis)
top80_90_Testis_overlap<-intersect(annotation$`Gene stable ID`[match(top80_90_Testis_Pig,annotation$`Pig gene stable ID`)],top80_90_Testis_Human)

last10_Testis_Human<-intersect(annotation$`Gene stable ID`[match(intersect(last10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Testis)
last10_Testis_Pig<-intersect(intersect(last10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Testis)
last10_Testis_overlap<-intersect(annotation$`Gene stable ID`[match(last10_Testis_Pig,annotation$`Pig gene stable ID`)],last10_Testis_Human)

#Uterus#
top10_Uterus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Uterus)
top10_Uterus_Pig<-intersect(intersect(top10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Uterus)
top10_Uterus_overlap<-intersect(annotation$`Gene stable ID`[match(top10_Uterus_Pig,annotation$`Pig gene stable ID`)],top10_Uterus_Human)

top10_20_Uterus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Uterus)
top10_20_Uterus_Pig<-intersect(intersect(top10_20$gene_id,annotation$`Pig gene stable ID`),pig_egene_Uterus)
top10_20_Uterus_overlap<-intersect(annotation$`Gene stable ID`[match(top10_20_Uterus_Pig,annotation$`Pig gene stable ID`)],top10_20_Uterus_Human)

top20_30_Uterus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Uterus)
top20_30_Uterus_Pig<-intersect(intersect(top20_30$gene_id,annotation$`Pig gene stable ID`),pig_egene_Uterus)
top20_30_Uterus_overlap<-intersect(annotation$`Gene stable ID`[match(top20_30_Uterus_Pig,annotation$`Pig gene stable ID`)],top20_30_Uterus_Human)

top30_40_Uterus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Uterus)
top30_40_Uterus_Pig<-intersect(intersect(top30_40$gene_id,annotation$`Pig gene stable ID`),pig_egene_Uterus)
top30_40_Uterus_overlap<-intersect(annotation$`Gene stable ID`[match(top30_40_Uterus_Pig,annotation$`Pig gene stable ID`)],top30_40_Uterus_Human)

top40_50_Uterus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Uterus)
top40_50_Uterus_Pig<-intersect(intersect(top40_50$gene_id,annotation$`Pig gene stable ID`),pig_egene_Uterus)
top40_50_Uterus_overlap<-intersect(annotation$`Gene stable ID`[match(top40_50_Uterus_Pig,annotation$`Pig gene stable ID`)],top40_50_Uterus_Human)

top50_60_Uterus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Uterus)
top50_60_Uterus_Pig<-intersect(intersect(top50_60$gene_id,annotation$`Pig gene stable ID`),pig_egene_Uterus)
top50_60_Uterus_overlap<-intersect(annotation$`Gene stable ID`[match(top50_60_Uterus_Pig,annotation$`Pig gene stable ID`)],top50_60_Uterus_Human)

top60_70_Uterus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Uterus)
top60_70_Uterus_Pig<-intersect(intersect(top60_70$gene_id,annotation$`Pig gene stable ID`),pig_egene_Uterus)
top60_70_Uterus_overlap<-intersect(annotation$`Gene stable ID`[match(top60_70_Uterus_Pig,annotation$`Pig gene stable ID`)],top60_70_Uterus_Human)

top70_80_Uterus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Uterus)
top70_80_Uterus_Pig<-intersect(intersect(top70_80$gene_id,annotation$`Pig gene stable ID`),pig_egene_Uterus)
top70_80_Uterus_overlap<-intersect(annotation$`Gene stable ID`[match(top70_80_Uterus_Pig,annotation$`Pig gene stable ID`)],top70_80_Uterus_Human)

top80_90_Uterus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Uterus)
top80_90_Uterus_Pig<-intersect(intersect(top80_90$gene_id,annotation$`Pig gene stable ID`),pig_egene_Uterus)
top80_90_Uterus_overlap<-intersect(annotation$`Gene stable ID`[match(top80_90_Uterus_Pig,annotation$`Pig gene stable ID`)],top80_90_Uterus_Human)

last10_Uterus_Human<-intersect(annotation$`Gene stable ID`[match(intersect(last10$gene_id,annotation$`Pig gene stable ID`),annotation$`Pig gene stable ID`)],hum_egene_Uterus)
last10_Uterus_Pig<-intersect(intersect(last10$gene_id,annotation$`Pig gene stable ID`),pig_egene_Uterus)
last10_Uterus_overlap<-intersect(annotation$`Gene stable ID`[match(last10_Uterus_Pig,annotation$`Pig gene stable ID`)],last10_Uterus_Human)

Result<-array(NA, dim=c(17,10))
rownames(Result)<-tissues
colnames(Result)<-c("top10%","10%-20%","20%-30%","30%-40%","40%-50%","50%-60%","60%-70%","70%-80%","80%-90%","bottom10%")
Result<-as.data.frame(Result)
for(i in 1:length(tissues)){
  Result[i,1]<-length(get(paste0('top10_',tissues[i],'_overlap')))
  Result[i,2]<-length(get(paste0('top10_20_',tissues[i],'_overlap')))
  Result[i,3]<-length(get(paste0('top20_30_',tissues[i],'_overlap')))
  Result[i,4]<-length(get(paste0('top30_40_',tissues[i],'_overlap')))
  Result[i,5]<-length(get(paste0('top40_50_',tissues[i],'_overlap')))
  Result[i,6]<-length(get(paste0('top50_60_',tissues[i],'_overlap')))
  Result[i,7]<-length(get(paste0('top60_70_',tissues[i],'_overlap')))
  Result[i,8]<-length(get(paste0('top70_80_',tissues[i],'_overlap')))
  Result[i,9]<-length(get(paste0('top80_90_',tissues[i],'_overlap')))
  Result[i,10]<-length(get(paste0('last10_',tissues[i],'_overlap')))
}
sum_tmp<-apply(Result,1,sum)
Result_normalized<-Result
for(i in 1:17){
  Result_normalized[i,]<-Result_normalized[i,]/sum_tmp[i]
}
Result_normalized<-as.matrix(Result_normalized)
#Result_percent<-Result/2900
#Result_percent
#Result_percent<-as.matrix(Result_percent)
col<-c("#CC66FF","#AAAAFF","#FF0000","#8EABD2","#FFD700","#33CCCC","#FDFDBF", "#8EA9DB","#8B0F55",
       "#E2EFDA","#7570B3","#AAEEFF","#99BB88","#FFDD99","#FF6600","#A6CEE3","#A6761D")
names(col)<-c("Adipose","Artery","Blood","Colon","Frontal_cortex","Heart","Hypothalamus","Ileum",
              "Kidney","Liver","Lung","Muscle","Ovary", "Pituitary","Spleen","Testis","Uterus")
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/Phastcons_sorted_normalized.tiff",##reqiured to change
     res = 300, width = 4500, height = 2500,compression = "lzw")
par(mar=c(8,8,6,1))
my<-barplot(Result_normalized,beside = T,col=col,cex.axis = 1.8,cex.names = 1.8,ylim = c(0,0.5))
title(ylab="Percentage of shared eGenes", xlab = "Windows of genes sorted by Phastcons",line=4, cex.lab=2.5)
legend("top", legend = c("Adipose","Artery","Blood","Colon","Frontal_cortex","Heart","Hypothalamus","Ileum",
                         "Kidney","Liver","Lung","Muscle","Ovary", "Pituitary","Spleen","Testis","Uterus"), 
       col = col, ncol = 4,xpd=T,
       bty = "n", pch=15 , pt.cex = 1.5, cex = 1.5, horiz = FALSE, inset = c(0.05, 0.05))
dev.off()


for(i in 1:length(tissues)){
  a<-get(paste0("eqtl_",tissues[i],"_hum"))$chr
  a<-as.data.frame(a)
  a$pos<-get(paste0("eqtl_",tissues[i],"_hum"))$variant_pos
  a$pos_1<-a$pos - 1
  a$chr<-substr(a$a,4,5)
  tmplist<-paste0(a$chr,"_",a$pos_1,"_",a$pos)
  write.table(tmplist,file=paste0("/Users/baizhonghao/Downloads/liftover/",tissues[i],".txt"),append = F,quote = F,sep = "\t",row.names = F,col.names = F)
}

#violin Phastcons#

eqtlfm<-as.data.frame(fread("/Users/baizhonghao/Downloads/Pig.chrAll.gene_coverage_phastCons.txt"))
eqtlfm$PhastCons[which(eqtlfm$PhastCons==".")]<-0
eqtlfm$PhastCons<-as.numeric(eqtlfm$PhastCons)
one2one_pig<-annotation[which(annotation$`Pig homology type`=='ortholog_one2one'),]
eqtlfm<-eqtlfm[match(intersect(eqtlfm$gene_id,one2one_pig$`Pig gene stable ID`),eqtlfm$gene_id),]
eqtlfm$gene_id<-one2one_pig$`Gene stable ID`[match(eqtlfm$gene_id, one2one_pig$`Pig gene stable ID`)]
Phastcons_sum<-NULL
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
  a_human$Phastcons<-eqtlfm$PhastCons[match(a_human$V1,eqtlfm$gene_id)]
  a_pig$Phastcons<-eqtlfm$PhastCons[match(a_pig$V1,eqtlfm$gene_id)]
  a_both$Phastcons<-eqtlfm$PhastCons[match(a_both$V1,eqtlfm$gene_id)]
  a_neither$Phastcons<-eqtlfm$PhastCons[match(a_neither$V1,eqtlfm$gene_id)]
  a_human$category<-'human-specific'
  a_pig$category<-'pig-specific'
  a_both$category<-'shared'
  a_neither$category<-'neither'
  a_sum<-rbind(a_human,a_pig,a_both,a_neither)
  a_sum$Tissues<-tissues[i]
  Phastcons_sum<-rbind(Phastcons_sum,a_sum)
}
Phastcons_sum$Phastcons<-as.numeric(round(Phastcons_sum$Phastcons,4))

Phastcons_mean<-NULL
for ( i in 1:length(tissues)){
  a_human<-read.table(paste0("/Users/baizhonghao/Downloads/human-pig GTEx/egene list/human_one2one_egenes_",tissues[i],".txt"))
  a_pig<-read.table(paste0("/Users/baizhonghao/Downloads/human-pig GTEx/egene list/pig_one2one_egenes_",tissues[i],".txt"))
  a_both<-read.table(paste0("/Users/baizhonghao/Downloads/human-pig GTEx/egene list/overlap_one2one_egenes_",tissues[i],".txt"))
  a_neither<-read.table(paste0("/Users/baizhonghao/Downloads/human-pig GTEx/egene list/neither_one2one_",tissues[i],".txt"))
  a_mean<-array(NA, dim=c(4,3))
  colnames(a_mean)<-c("Phastcons","category","Tissues")
  a_mean<-as.data.frame(a_mean)
  humanonly<-setdiff(a_human$V1,intersect(a_human$V1,a_both$V1))
  pigonly<-setdiff(a_pig$V1,intersect(a_pig$V1,a_both$V1))
  a_human<-as.data.frame(a_human[match(humanonly,a_human$V1),])
  a_pig<-as.data.frame(a_pig[match(pigonly,a_pig$V1),])
  colnames(a_human)<-'V1'
  colnames(a_pig)<-'V1'
  a_mean$Phastcons[1]<-mean(eqtlfm$PhastCons[match(intersect(a_human$V1,eqtlfm$gene_id),eqtlfm$gene_id)])
  a_mean$Phastcons[2]<-mean(eqtlfm$PhastCons[match(intersect(a_pig$V1,eqtlfm$gene_id),eqtlfm$gene_id)])
  a_mean$Phastcons[3]<-mean(eqtlfm$PhastCons[match(intersect(a_both$V1,eqtlfm$gene_id),eqtlfm$gene_id)])
  a_mean$Phastcons[4]<-mean(eqtlfm$PhastCons[match(intersect(a_neither$V1,eqtlfm$gene_id),eqtlfm$gene_id)])
  a_mean$category[1]<-'human-specific'
  a_mean$category[2]<-'pig-specific'
  a_mean$category[3]<-'shared'
  a_mean$category[4]<-'neither'
  a_mean$Tissues<-tissues[i]
  Phastcons_mean<-rbind(Phastcons_mean,a_mean)
}
phastcons_mean
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/eqtl/Phastconsviolinplot-category.tiff",##reqiured to change
     res = 300, width =6000, height = 2000,compression = "lzw")

ggplot(Phastcons_sum,mapping = aes(x=Tissues,y=Phastcons,fill=category))+
  geom_violin(aes(fill = category), trim = FALSE) +
  geom_boxplot(width = 0.25, aes(fill=category),position=position_dodge(0.9))+
  theme_bw()+
  labs(x = 'Tissues',y = 'Phastcons') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title =element_text(size = 12),axis.text =element_text(size = 6, color = 'black'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=5.2))+
  theme(axis.title.y.left = element_text(vjust = 2))+
  theme(axis.text.x = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(axis.text=element_text(color="black"))+
  theme(axis.title=element_text(face = "bold", color="black"))
dev.off()

col<-c("#CC66FF","#AAAAFF","#FF0000","#8EABD2","#FFD700","#33CCCC","#FDFDBF", "#8EA9DB","#8B0F55",
       "#E2EFDA","#7570B3","#AAEEFF","#99BB88","#FFDD99","#FF6600","#A6CEE3","#A6761D")
names(col)<-c("Adipose","Artery","Blood","Colon","Frontal_cortex","Heart","Hypothalamus","Ileum",
              "Kidney","Liver","Lung","Muscle","Ovary", "Pituitary","Spleen","Testis","Uterus")
Phastcons_sum$category<-factor(Phastcons_sum$category,levels=c("neither","human_only","pig_only","both"))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/eqtl/Phastconsviolinplot-tissues-shorted.tiff",##reqiured to change
     res = 300, width =4000, height = 1000,compression = "lzw")

ggplot(Phastcons_sum,mapping = aes(x=category,y=Phastcons,fill=Tissues))+
  geom_violin(aes(fill = Tissues), trim = FALSE) +
  geom_boxplot(width = 0.25, aes(fill=Tissues),position=position_dodge(0.9))+
  theme_bw()+
  labs(x = 'Categories',y = 'Phastcons') +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title =element_text(size = 12,vjust=2),axis.text =element_text(size = 6, color = 'black'))+
  theme(axis.text.x = element_text(size=5.2))+
  theme(axis.title.y.left = element_text(vjust = 2))+
  theme(axis.text.x = element_text(size=16))+
  theme(axis.text.y = element_text(size=16))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text=element_text(color="black"))+
  theme(axis.title=element_text(face = "bold", color="black"))+
  theme(legend.position = "none")+
  scale_fill_manual(values = col)
dev.off()

BiocManager::install("ggsignif")
library(ggsignif)
Phastcons_mean$category<-factor(Phastcons_mean$category,levels = c("neither","human-specific","pig-specific","shared"))

Phastconsvio<-ggplot(Phastcons_mean,mapping = aes(x=category,y=Phastcons,fill='red'))+
  geom_violin(aes(fill = 'red'), trim = FALSE) +
  geom_boxplot(width = 0.15)+
  geom_signif(comparisons = list(c("neither","human-specific"),
                                 c("neither","pig-specific"),
                                 c("neither","shared"),
                                 #c("human_only","pig_only"),
                                 c("human-specific","shared"),
                                 c("pig-specific","shared")),
              y_position = c(0.24,0.25,0.26,0.27,0.28),
              map_signif_level = TRUE)+
  theme_classic()+
  labs(x = NULL,y = 'Phastcons') +
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
ggsave(Phastconsvio,file='/Users/baizhonghao/Downloads/Phastconsviolin.pdf',dpi=300,width=6,height=6)