#Muscle_SNP_overlaploci#
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

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_chr1$variant_id,split="_"))))
eqtl_Muscle_chr1$chr<-tmp$V1
eqtl_Muscle_chr1$loci<-tmp$V2

eqtl_Muscle_chr1$index<-paste0(eqtl_Muscle_chr1$chr,"-",eqtl_Muscle_chr1$loci)
Muscle_chr1$index<-paste0(Muscle_chr1$chr,"-",Muscle_chr1$pig_pos)
same_Muscle_chr1<-intersect(Muscle_chr1$pig_pos,eqtl_Muscle_chr1$loci)


tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_chr2$variant_id,split="_"))))
eqtl_Muscle_chr2$chr<-tmp$V1
eqtl_Muscle_chr2$loci<-tmp$V2

eqtl_Muscle_chr2$index<-paste0(eqtl_Muscle_chr2$chr,"-",eqtl_Muscle_chr2$loci)
Muscle_chr2$index<-paste0(Muscle_chr2$chr,"-",Muscle_chr2$pig_pos)
same_Muscle_chr2<-intersect(Muscle_chr2$pig_pos,eqtl_Muscle_chr2$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_chr3$variant_id,split="_"))))
eqtl_Muscle_chr3$chr<-tmp$V1
eqtl_Muscle_chr3$loci<-tmp$V2

eqtl_Muscle_chr3$index<-paste0(eqtl_Muscle_chr3$chr,"-",eqtl_Muscle_chr3$loci)
Muscle_chr3$index<-paste0(Muscle_chr3$chr,"-",Muscle_chr3$pig_pos)
same_Muscle_chr3<-intersect(Muscle_chr3$pig_pos,eqtl_Muscle_chr3$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_chr4$variant_id,split="_"))))
eqtl_Muscle_chr4$chr<-tmp$V1
eqtl_Muscle_chr4$loci<-tmp$V2

eqtl_Muscle_chr4$index<-paste0(eqtl_Muscle_chr4$chr,"-",eqtl_Muscle_chr4$loci)
Muscle_chr4$index<-paste0(Muscle_chr4$chr,"-",Muscle_chr4$pig_pos)
same_Muscle_chr4<-intersect(Muscle_chr4$pig_pos,eqtl_Muscle_chr4$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_chr5$variant_id,split="_"))))
eqtl_Muscle_chr5$loci<-tmp$V2
same_Muscle_chr5<-intersect(Muscle_chr5$pig_pos,eqtl_Muscle_chr5$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_chr6$variant_id,split="_"))))
eqtl_Muscle_chr6$loci<-tmp$V2
same_Muscle_chr6<-intersect(Muscle_chr6$pig_pos,eqtl_Muscle_chr6$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_chr7$variant_id,split="_"))))
eqtl_Muscle_chr7$loci<-tmp$V2
same_Muscle_chr7<-intersect(Muscle_chr7$pig_pos,eqtl_Muscle_chr7$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_chr8$variant_id,split="_"))))
eqtl_Muscle_chr8$loci<-tmp$V2
same_Muscle_chr8<-intersect(Muscle_chr8$pig_pos,eqtl_Muscle_chr8$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_chr9$variant_id,split="_"))))
eqtl_Muscle_chr9$loci<-tmp$V2
same_Muscle_chr9<-intersect(Muscle_chr9$pig_pos,eqtl_Muscle_chr9$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_chr10$variant_id,split="_"))))
eqtl_Muscle_chr10$loci<-tmp$V2
same_Muscle_chr10<-intersect(Muscle_chr10$pig_pos,eqtl_Muscle_chr10$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_chr11$variant_id,split="_"))))
eqtl_Muscle_chr11$loci<-tmp$V2
same_Muscle_chr11<-intersect(Muscle_chr11$pig_pos,eqtl_Muscle_chr11$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_chr12$variant_id,split="_"))))
eqtl_Muscle_chr12$loci<-tmp$V2
same_Muscle_chr12<-intersect(Muscle_chr12$pig_pos,eqtl_Muscle_chr12$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_chr13$variant_id,split="_"))))
eqtl_Muscle_chr13$loci<-tmp$V2
same_Muscle_chr13<-intersect(Muscle_chr13$pig_pos,eqtl_Muscle_chr13$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_chr14$variant_id,split="_"))))
eqtl_Muscle_chr14$loci<-tmp$V2
same_Muscle_chr14<-intersect(Muscle_chr14$pig_pos,eqtl_Muscle_chr14$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_chr15$variant_id,split="_"))))
eqtl_Muscle_chr15$loci<-tmp$V2
same_Muscle_chr15<-intersect(Muscle_chr15$pig_pos,eqtl_Muscle_chr15$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_chr16$variant_id,split="_"))))
eqtl_Muscle_chr16$loci<-tmp$V2
same_Muscle_chr16<-intersect(Muscle_chr16$pig_pos,eqtl_Muscle_chr16$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_chr17$variant_id,split="_"))))
eqtl_Muscle_chr17$loci<-tmp$V2
same_Muscle_chr17<-intersect(Muscle_chr17$pig_pos,eqtl_Muscle_chr17$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Muscle_chr18$variant_id,split="_"))))
eqtl_Muscle_chr18$loci<-tmp$V2
same_Muscle_chr18<-intersect(Muscle_chr18$pig_pos,eqtl_Muscle_chr18$loci)

chr1_Muscle_genes<-NULL
chr1_Muscle_overloci<-NULL
if(length(same_Muscle_chr1!=0)){
  for(i in 1:length(same_Muscle_chr1)){
    a<-eqtl_Muscle_chr1$phenotype_id[grep(same_Muscle_chr1[i],eqtl_Muscle_chr1$loci)]
    b<-eqtl_Muscle_chr1[grep(same_Muscle_chr1[i],eqtl_Muscle_chr1$loci),]
    chr1_Muscle_genes<-c(chr1_Muscle_genes,a)
    chr1_Muscle_overloci<-rbind(chr1_Muscle_overloci,b)
  }
}

if(length(same_Muscle_chr2!=0)){
  chr2_Muscle_genes<-NULL
  chr2_Muscle_overloci<-NULL
  for(i in 1:length(same_Muscle_chr2)){
    a<-eqtl_Muscle_chr2$phenotype_id[grep(same_Muscle_chr2[i],eqtl_Muscle_chr2$loci)]
    b<-eqtl_Muscle_chr2[grep(same_Muscle_chr2[i],eqtl_Muscle_chr2$loci),]
    chr2_Muscle_genes<-c(chr2_Muscle_genes,a)
    chr2_Muscle_overloci<-rbind(chr2_Muscle_overloci,b)
  }
}

if(length(same_Muscle_chr3!=0)){
  chr3_Muscle_genes<-NULL
  chr3_Muscle_overloci<-NULL
  for(i in 1:length(same_Muscle_chr3)){
    a<-eqtl_Muscle_chr3$phenotype_id[grep(same_Muscle_chr3[i],eqtl_Muscle_chr3$loci)]
    b<-eqtl_Muscle_chr3[grep(same_Muscle_chr3[i],eqtl_Muscle_chr3$loci),]
    chr3_Muscle_genes<-c(chr3_Muscle_genes,a)
    chr3_Muscle_overloci<-rbind(chr3_Muscle_overloci,b)
  }
}

if(length(same_Muscle_chr4!=0)){
  chr4_Muscle_genes<-NULL
  chr4_Muscle_overloci<-NULL
  for(i in 1:length(same_Muscle_chr4)){
    a<-eqtl_Muscle_chr4$phenotype_id[grep(same_Muscle_chr4[i],eqtl_Muscle_chr4$loci)]
    b<-eqtl_Muscle_chr4[grep(same_Muscle_chr4[i],eqtl_Muscle_chr4$loci),]
    chr4_Muscle_genes<-c(chr4_Muscle_genes,a)
    chr4_Muscle_overloci<-rbind(chr4_Muscle_overloci,b)
  }
}

if(length(same_Muscle_chr5!=0)){
  chr5_Muscle_genes<-NULL
  chr5_Muscle_overloci<-NULL
  for(i in 1:length(same_Muscle_chr5)){
    a<-eqtl_Muscle_chr5$phenotype_id[grep(same_Muscle_chr5[i],eqtl_Muscle_chr5$loci)]
    b<-eqtl_Muscle_chr5[grep(same_Muscle_chr5[i],eqtl_Muscle_chr5$loci),]
    chr5_Muscle_genes<-c(chr5_Muscle_genes,a)
    chr5_Muscle_overloci<-rbind(chr5_Muscle_overloci,b)
  }
}

if(length(same_Muscle_chr6!=0)){
  chr6_Muscle_genes<-NULL
  chr6_Muscle_overloci<-NULL
  for(i in 1:length(same_Muscle_chr6)){
    a<-eqtl_Muscle_chr6$phenotype_id[grep(same_Muscle_chr6[i],eqtl_Muscle_chr6$loci)]
    b<-eqtl_Muscle_chr6[grep(same_Muscle_chr6[i],eqtl_Muscle_chr6$loci),]
    chr6_Muscle_genes<-c(chr6_Muscle_genes,a)
    chr6_Muscle_overloci<-rbind(chr6_Muscle_overloci,b)
  }
}

if(length(same_Muscle_chr7!=0)){
  chr7_Muscle_genes<-NULL
  chr7_Muscle_overloci<-NULL
  for(i in 1:length(same_Muscle_chr7)){
    a<-eqtl_Muscle_chr7$phenotype_id[grep(same_Muscle_chr7[i],eqtl_Muscle_chr7$loci)]
    b<-eqtl_Muscle_chr7[grep(same_Muscle_chr7[i],eqtl_Muscle_chr7$loci),]
    chr7_Muscle_genes<-c(chr7_Muscle_genes,a)
    chr7_Muscle_overloci<-rbind(chr7_Muscle_overloci,b)
  }
}

if(length(same_Muscle_chr8!=0)){
  chr8_Muscle_genes<-NULL
  chr8_Muscle_overloci<-NULL
  for(i in 1:length(same_Muscle_chr8)){
    a<-eqtl_Muscle_chr8$phenotype_id[grep(same_Muscle_chr8[i],eqtl_Muscle_chr8$loci)]
    b<-eqtl_Muscle_chr8[grep(same_Muscle_chr8[i],eqtl_Muscle_chr8$loci),]
    chr8_Muscle_genes<-c(chr8_Muscle_genes,a)
    chr8_Muscle_overloci<-rbind(chr8_Muscle_overloci,b)
  }
}

if(length(same_Muscle_chr9!=0)){
  chr9_Muscle_genes<-NULL
  chr9_Muscle_overloci<-NULL
  for(i in 1:length(same_Muscle_chr9)){
    a<-eqtl_Muscle_chr9$phenotype_id[grep(same_Muscle_chr9[i],eqtl_Muscle_chr9$loci)]
    b<-eqtl_Muscle_chr9[grep(same_Muscle_chr9[i],eqtl_Muscle_chr9$loci),]
    chr9_Muscle_genes<-c(chr9_Muscle_genes,a)
    chr9_Muscle_overloci<-rbind(chr9_Muscle_overloci,b)
  }
}

if(length(same_Muscle_chr10!=0)){
  chr10_Muscle_genes<-NULL
  chr10_Muscle_overloci<-NULL
  for(i in 1:length(same_Muscle_chr10)){
    a<-eqtl_Muscle_chr10$phenotype_id[grep(same_Muscle_chr10[i],eqtl_Muscle_chr10$loci)]
    b<-eqtl_Muscle_chr10[grep(same_Muscle_chr10[i],eqtl_Muscle_chr10$loci),]
    chr10_Muscle_genes<-c(chr10_Muscle_genes,a)
    chr10_Muscle_overloci<-rbind(chr10_Muscle_overloci,b)
  }
}

if(length(same_Muscle_chr11!=0)){
  chr11_Muscle_genes<-NULL
  chr11_Muscle_overloci<-NULL
  for(i in 1:length(same_Muscle_chr11)){
    a<-eqtl_Muscle_chr11$phenotype_id[grep(same_Muscle_chr11[i],eqtl_Muscle_chr11$loci)]
    b<-eqtl_Muscle_chr11[grep(same_Muscle_chr11[i],eqtl_Muscle_chr11$loci),]
    chr11_Muscle_genes<-c(chr11_Muscle_genes,a)
    chr11_Muscle_overloci<-rbind(chr11_Muscle_overloci,b)
  }
}

chr12_Muscle_genes<-NULL
chr12_Muscle_overloci<-NULL
if(length(same_Muscle_chr12!=0)){
  for(i in 1:length(same_Muscle_chr12)){
    a<-eqtl_Muscle_chr12$phenotype_id[grep(same_Muscle_chr12[i],eqtl_Muscle_chr12$loci)]
    b<-eqtl_Muscle_chr12[grep(same_Muscle_chr12[i],eqtl_Muscle_chr12$loci),]
    chr12_Muscle_genes<-c(chr12_Muscle_genes,a)
    chr12_Muscle_overloci<-rbind(chr12_Muscle_overloci,b)
  }
}

chr13_Muscle_genes<-NULL
chr13_Muscle_overloci<-NULL
if(length(same_Muscle_chr13!=0)){
  for(i in 1:length(same_Muscle_chr13)){
    a<-eqtl_Muscle_chr13$phenotype_id[grep(same_Muscle_chr13[i],eqtl_Muscle_chr13$loci)]
    b<-eqtl_Muscle_chr13[grep(same_Muscle_chr13[i],eqtl_Muscle_chr13$loci),]
    chr13_Muscle_genes<-c(chr13_Muscle_genes,a)
    chr13_Muscle_overloci<-rbind(chr13_Muscle_overloci,b)
  }
}

chr14_Muscle_genes<-NULL
chr14_Muscle_overloci<-NULL
if(length(same_Muscle_chr14!=0)){
  for(i in 1:length(same_Muscle_chr14)){
    a<-eqtl_Muscle_chr14$phenotype_id[grep(same_Muscle_chr14[i],eqtl_Muscle_chr14$loci)]
    b<-eqtl_Muscle_chr14[grep(same_Muscle_chr14[i],eqtl_Muscle_chr14$loci),]
    chr14_Muscle_genes<-c(chr14_Muscle_genes,a)
    chr14_Muscle_overloci<-rbind(chr14_Muscle_overloci,b)
  }
}

chr15_Muscle_genes<-NULL
chr15_Muscle_overloci<-NULL
if(length(same_Muscle_chr15!=0)){
  for(i in 1:length(same_Muscle_chr15)){
    a<-eqtl_Muscle_chr15$phenotype_id[grep(same_Muscle_chr15[i],eqtl_Muscle_chr15$loci)]
    b<-eqtl_Muscle_chr15[grep(same_Muscle_chr15[i],eqtl_Muscle_chr15$loci),]
    chr15_Muscle_genes<-c(chr15_Muscle_genes,a)
    chr15_Muscle_overloci<-rbind(chr15_Muscle_overloci,b)
  }
}

chr16_Muscle_genes<-NULL
chr16_Muscle_overloci<-NULL
if(length(same_Muscle_chr16!=0)){
  for(i in 1:length(same_Muscle_chr16)){
    a<-eqtl_Muscle_chr16$phenotype_id[grep(same_Muscle_chr16[i],eqtl_Muscle_chr16$loci)]
    b<-eqtl_Muscle_chr16[grep(same_Muscle_chr16[i],eqtl_Muscle_chr16$loci),]
    chr16_Muscle_genes<-c(chr16_Muscle_genes,a)
    chr16_Muscle_overloci<-rbind(chr16_Muscle_overloci,b)
  }
}

chr17_Muscle_genes<-NULL
chr17_Muscle_overloci<-NULL
if(length(same_Muscle_chr17!=0)){
  for(i in 1:length(same_Muscle_chr17)){
    a<-eqtl_Muscle_chr17$phenotype_id[grep(same_Muscle_chr17[i],eqtl_Muscle_chr17$loci)]
    b<-eqtl_Muscle_chr17[grep(same_Muscle_chr17[i],eqtl_Muscle_chr17$loci),]
    chr17_Muscle_genes<-c(chr17_Muscle_genes,a)
    chr17_Muscle_overloci<-rbind(chr17_Muscle_overloci,b)
  }
}

chr18_Muscle_genes<-NULL
chr18_Muscle_overloci<-NULL
if(length(same_Muscle_chr18!=0)){
  for(i in 1:length(same_Muscle_chr18)){
    a<-eqtl_Muscle_chr18$phenotype_id[grep(same_Muscle_chr18[i],eqtl_Muscle_chr18$loci)]
    b<-eqtl_Muscle_chr18[grep(same_Muscle_chr18[i],eqtl_Muscle_chr18$loci),]
    chr18_Muscle_genes<-c(chr18_Muscle_genes,a)
    chr18_Muscle_overloci<-rbind(chr18_Muscle_overloci,b)
  }
}

one2one_pig<-annotation$`Pig gene stable ID`[which(annotation$`Pig homology type`=="ortholog_one2one")]

chr1_Muscle_overloci_one2one<-chr1_Muscle_overloci[match(intersect(chr1_Muscle_genes,one2one_pig),chr1_Muscle_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Muscle_chr1$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Muscle_chr1$hum_chr<-tmp$V1
Muscle_chr1$hum_loci<-tmp$V3
FM_Muscle_hum$index<-paste0(FM_Muscle_hum$chr,"_",FM_Muscle_hum$variant_pos)
Muscle_chr1$index<-paste0("chr",Muscle_chr1$hum_chr,"_",Muscle_chr1$hum_loci)
chr1_Muscle_hum<-NULL
for(i in 1:length(same_Muscle_chr1)){
  a<-Muscle_chr1[match(same_Muscle_chr1[i],Muscle_chr1$pig_pos),]
  b<-FM_Muscle_hum[match(intersect(a$index,FM_Muscle_hum$index),FM_Muscle_hum$index),]
  chr1_Muscle_hum<-rbind(chr1_Muscle_hum,b)
}
chr1_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr1_Muscle_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr1_Muscle_hum_one2one<-chr1_Muscle_hum[match(intersect(chr1_pig2hum_one2one,chr1_Muscle_hum$gene_id),chr1_Muscle_hum$gene_id),]

chr2_Muscle_overloci_one2one<-chr2_Muscle_overloci[match(intersect(chr2_Muscle_genes,one2one_pig),chr2_Muscle_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Muscle_chr2$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Muscle_chr2$hum_chr<-tmp$V1
Muscle_chr2$hum_loci<-tmp$V3
FM_Muscle_hum$index<-paste0(FM_Muscle_hum$chr,"_",FM_Muscle_hum$variant_pos)
Muscle_chr2$index<-paste0("chr",Muscle_chr2$hum_chr,"_",Muscle_chr2$hum_loci)
chr2_Muscle_hum<-NULL
for(i in 1:length(same_Muscle_chr2)){
  a<-Muscle_chr2[match(same_Muscle_chr2[i],Muscle_chr2$pig_pos),]
  b<-FM_Muscle_hum[match(intersect(a$index,FM_Muscle_hum$index),FM_Muscle_hum$index),]
  chr2_Muscle_hum<-rbind(chr2_Muscle_hum,b)
}
chr2_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr2_Muscle_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr2_Muscle_hum_one2one<-chr2_Muscle_hum[match(intersect(chr2_pig2hum_one2one,chr2_Muscle_hum$gene_id),chr2_Muscle_hum$gene_id),]

chr3_Muscle_overloci_one2one<-chr3_Muscle_overloci[match(intersect(chr3_Muscle_genes,one2one_pig),chr3_Muscle_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Muscle_chr3$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Muscle_chr3$hum_chr<-tmp$V1
Muscle_chr3$hum_loci<-tmp$V3
FM_Muscle_hum$index<-paste0(FM_Muscle_hum$chr,"_",FM_Muscle_hum$variant_pos)
Muscle_chr3$index<-paste0("chr",Muscle_chr3$hum_chr,"_",Muscle_chr3$hum_loci)
chr3_Muscle_hum<-NULL
for(i in 1:length(same_Muscle_chr3)){
  a<-Muscle_chr3[match(same_Muscle_chr3[i],Muscle_chr3$pig_pos),]
  b<-FM_Muscle_hum[match(intersect(a$index,FM_Muscle_hum$index),FM_Muscle_hum$index),]
  chr3_Muscle_hum<-rbind(chr3_Muscle_hum,b)
}
chr3_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr3_Muscle_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr3_Muscle_hum_one2one<-chr3_Muscle_hum[match(intersect(chr3_pig2hum_one2one,chr3_Muscle_hum$gene_id),chr3_Muscle_hum$gene_id),]

chr4_Muscle_overloci_one2one<-chr4_Muscle_overloci[match(intersect(chr4_Muscle_genes,one2one_pig),chr4_Muscle_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Muscle_chr4$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Muscle_chr4$hum_chr<-tmp$V1
Muscle_chr4$hum_loci<-tmp$V3
FM_Muscle_hum$index<-paste0(FM_Muscle_hum$chr,"_",FM_Muscle_hum$variant_pos)
Muscle_chr4$index<-paste0("chr",Muscle_chr4$hum_chr,"_",Muscle_chr4$hum_loci)
chr4_Muscle_hum<-NULL
for(i in 1:length(same_Muscle_chr4)){
  a<-Muscle_chr4[match(same_Muscle_chr4[i],Muscle_chr4$pig_pos),]
  b<-FM_Muscle_hum[match(intersect(a$index,FM_Muscle_hum$index),FM_Muscle_hum$index),]
  chr4_Muscle_hum<-rbind(chr4_Muscle_hum,b)
}
chr4_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr4_Muscle_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr4_Muscle_hum_one2one<-chr4_Muscle_hum[match(intersect(chr4_pig2hum_one2one,chr4_Muscle_hum$gene_id),chr4_Muscle_hum$gene_id),]

chr5_Muscle_overloci_one2one<-chr5_Muscle_overloci[match(intersect(chr5_Muscle_genes,one2one_pig),chr5_Muscle_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Muscle_chr5$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Muscle_chr5$hum_chr<-tmp$V1
Muscle_chr5$hum_loci<-tmp$V3
FM_Muscle_hum$index<-paste0(FM_Muscle_hum$chr,"_",FM_Muscle_hum$variant_pos)
Muscle_chr5$index<-paste0("chr",Muscle_chr5$hum_chr,"_",Muscle_chr5$hum_loci)
chr5_Muscle_hum<-NULL
for(i in 1:length(same_Muscle_chr5)){
  a<-Muscle_chr5[match(same_Muscle_chr5[i],Muscle_chr5$pig_pos),]
  b<-FM_Muscle_hum[match(intersect(a$index,FM_Muscle_hum$index),FM_Muscle_hum$index),]
  chr5_Muscle_hum<-rbind(chr5_Muscle_hum,b)
}
chr5_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr5_Muscle_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr5_Muscle_hum_one2one<-chr5_Muscle_hum[match(intersect(chr5_pig2hum_one2one,chr5_Muscle_hum$gene_id),chr5_Muscle_hum$gene_id),]

chr6_Muscle_overloci_one2one<-chr6_Muscle_overloci[match(intersect(chr6_Muscle_genes,one2one_pig),chr6_Muscle_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Muscle_chr6$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Muscle_chr6$hum_chr<-tmp$V1
Muscle_chr6$hum_loci<-tmp$V3
FM_Muscle_hum$index<-paste0(FM_Muscle_hum$chr,"_",FM_Muscle_hum$variant_pos)
Muscle_chr6$index<-paste0("chr",Muscle_chr6$hum_chr,"_",Muscle_chr6$hum_loci)
chr6_Muscle_hum<-NULL
for(i in 1:length(same_Muscle_chr6)){
  a<-Muscle_chr6[match(same_Muscle_chr6[i],Muscle_chr6$pig_pos),]
  b<-FM_Muscle_hum[match(intersect(a$index,FM_Muscle_hum$index),FM_Muscle_hum$index),]
  chr6_Muscle_hum<-rbind(chr6_Muscle_hum,b)
}
chr6_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr6_Muscle_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr6_Muscle_hum_one2one<-chr6_Muscle_hum[match(intersect(chr6_pig2hum_one2one,chr6_Muscle_hum$gene_id),chr6_Muscle_hum$gene_id),]

chr7_Muscle_overloci_one2one<-chr7_Muscle_overloci[match(intersect(chr7_Muscle_genes,one2one_pig),chr7_Muscle_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Muscle_chr7$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Muscle_chr7$hum_chr<-tmp$V1
Muscle_chr7$hum_loci<-tmp$V3
FM_Muscle_hum$index<-paste0(FM_Muscle_hum$chr,"_",FM_Muscle_hum$variant_pos)
Muscle_chr7$index<-paste0("chr",Muscle_chr7$hum_chr,"_",Muscle_chr7$hum_loci)
chr7_Muscle_hum<-NULL
for(i in 1:length(same_Muscle_chr7)){
  a<-Muscle_chr7[match(same_Muscle_chr7[i],Muscle_chr7$pig_pos),]
  b<-FM_Muscle_hum[match(intersect(a$index,FM_Muscle_hum$index),FM_Muscle_hum$index),]
  chr7_Muscle_hum<-rbind(chr7_Muscle_hum,b)
}
chr7_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr7_Muscle_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr7_Muscle_hum_one2one<-chr7_Muscle_hum[match(intersect(chr7_pig2hum_one2one,chr7_Muscle_hum$gene_id),chr7_Muscle_hum$gene_id),]

chr8_Muscle_overloci_one2one<-chr8_Muscle_overloci[match(intersect(chr8_Muscle_genes,one2one_pig),chr8_Muscle_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Muscle_chr8$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Muscle_chr8$hum_chr<-tmp$V1
Muscle_chr8$hum_loci<-tmp$V3
FM_Muscle_hum$index<-paste0(FM_Muscle_hum$chr,"_",FM_Muscle_hum$variant_pos)
Muscle_chr8$index<-paste0("chr",Muscle_chr8$hum_chr,"_",Muscle_chr8$hum_loci)
chr8_Muscle_hum<-NULL
for(i in 1:length(same_Muscle_chr8)){
  a<-Muscle_chr8[match(same_Muscle_chr8[i],Muscle_chr8$pig_pos),]
  b<-FM_Muscle_hum[match(intersect(a$index,FM_Muscle_hum$index),FM_Muscle_hum$index),]
  chr8_Muscle_hum<-rbind(chr8_Muscle_hum,b)
}
chr8_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr8_Muscle_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr8_Muscle_hum_one2one<-chr8_Muscle_hum[match(intersect(chr8_pig2hum_one2one,chr8_Muscle_hum$gene_id),chr8_Muscle_hum$gene_id),]

chr9_Muscle_overloci_one2one<-chr9_Muscle_overloci[match(intersect(chr9_Muscle_genes,one2one_pig),chr9_Muscle_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Muscle_chr9$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Muscle_chr9$hum_chr<-tmp$V1
Muscle_chr9$hum_loci<-tmp$V3
FM_Muscle_hum$index<-paste0(FM_Muscle_hum$chr,"_",FM_Muscle_hum$variant_pos)
Muscle_chr9$index<-paste0("chr",Muscle_chr9$hum_chr,"_",Muscle_chr9$hum_loci)
chr9_Muscle_hum<-NULL
for(i in 1:length(same_Muscle_chr9)){
  a<-Muscle_chr9[match(same_Muscle_chr9[i],Muscle_chr9$pig_pos),]
  b<-FM_Muscle_hum[match(intersect(a$index,FM_Muscle_hum$index),FM_Muscle_hum$index),]
  chr9_Muscle_hum<-rbind(chr9_Muscle_hum,b)
}
chr9_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr9_Muscle_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr9_Muscle_hum_one2one<-chr9_Muscle_hum[match(intersect(chr9_pig2hum_one2one,chr9_Muscle_hum$gene_id),chr9_Muscle_hum$gene_id),]

chr10_Muscle_overloci_one2one<-chr10_Muscle_overloci[match(intersect(chr10_Muscle_genes,one2one_pig),chr10_Muscle_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Muscle_chr10$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Muscle_chr10$hum_chr<-tmp$V1
Muscle_chr10$hum_loci<-tmp$V3
FM_Muscle_hum$index<-paste0(FM_Muscle_hum$chr,"_",FM_Muscle_hum$variant_pos)
Muscle_chr10$index<-paste0("chr",Muscle_chr10$hum_chr,"_",Muscle_chr10$hum_loci)
chr10_Muscle_hum<-NULL
for(i in 1:length(same_Muscle_chr10)){
  a<-Muscle_chr10[match(same_Muscle_chr10[i],Muscle_chr10$pig_pos),]
  b<-FM_Muscle_hum[match(intersect(a$index,FM_Muscle_hum$index),FM_Muscle_hum$index),]
  chr10_Muscle_hum<-rbind(chr10_Muscle_hum,b)
}
chr10_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr10_Muscle_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr10_Muscle_hum_one2one<-chr10_Muscle_hum[match(intersect(chr10_pig2hum_one2one,chr10_Muscle_hum$gene_id),chr10_Muscle_hum$gene_id),]

chr11_Muscle_overloci_one2one<-chr11_Muscle_overloci[match(intersect(chr11_Muscle_genes,one2one_pig),chr11_Muscle_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Muscle_chr11$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Muscle_chr11$hum_chr<-tmp$V1
Muscle_chr11$hum_loci<-tmp$V3
FM_Muscle_hum$index<-paste0(FM_Muscle_hum$chr,"_",FM_Muscle_hum$variant_pos)
Muscle_chr11$index<-paste0("chr",Muscle_chr11$hum_chr,"_",Muscle_chr11$hum_loci)
chr11_Muscle_hum<-NULL
for(i in 1:length(same_Muscle_chr11)){
  a<-Muscle_chr11[match(same_Muscle_chr11[i],Muscle_chr11$pig_pos),]
  b<-FM_Muscle_hum[match(intersect(a$index,FM_Muscle_hum$index),FM_Muscle_hum$index),]
  chr11_Muscle_hum<-rbind(chr11_Muscle_hum,b)
}
chr11_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr11_Muscle_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr11_Muscle_hum_one2one<-chr11_Muscle_hum[match(intersect(chr11_pig2hum_one2one,chr11_Muscle_hum$gene_id),chr11_Muscle_hum$gene_id),]

chr12_Muscle_overloci_one2one<-chr12_Muscle_overloci[match(intersect(chr12_Muscle_genes,one2one_pig),chr12_Muscle_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Muscle_chr12$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Muscle_chr12$hum_chr<-tmp$V1
Muscle_chr12$hum_loci<-tmp$V3
FM_Muscle_hum$index<-paste0(FM_Muscle_hum$chr,"_",FM_Muscle_hum$variant_pos)
Muscle_chr12$index<-paste0("chr",Muscle_chr12$hum_chr,"_",Muscle_chr12$hum_loci)
chr12_Muscle_hum<-NULL
for(i in 1:length(same_Muscle_chr12)){
  a<-Muscle_chr12[match(same_Muscle_chr12[i],Muscle_chr12$pig_pos),]
  b<-FM_Muscle_hum[match(intersect(a$index,FM_Muscle_hum$index),FM_Muscle_hum$index),]
  chr12_Muscle_hum<-rbind(chr12_Muscle_hum,b)
}
chr12_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr12_Muscle_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr12_Muscle_hum_one2one<-chr12_Muscle_hum[match(intersect(chr12_pig2hum_one2one,chr12_Muscle_hum$gene_id),chr12_Muscle_hum$gene_id),]

chr13_Muscle_overloci_one2one<-chr13_Muscle_overloci[match(intersect(chr13_Muscle_genes,one2one_pig),chr13_Muscle_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Muscle_chr13$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Muscle_chr13$hum_chr<-tmp$V1
Muscle_chr13$hum_loci<-tmp$V3
FM_Muscle_hum$index<-paste0(FM_Muscle_hum$chr,"_",FM_Muscle_hum$variant_pos)
Muscle_chr13$index<-paste0("chr",Muscle_chr13$hum_chr,"_",Muscle_chr13$hum_loci)
chr13_Muscle_hum<-NULL
for(i in 1:length(same_Muscle_chr13)){
  a<-Muscle_chr13[match(same_Muscle_chr13[i],Muscle_chr13$pig_pos),]
  b<-FM_Muscle_hum[match(intersect(a$index,FM_Muscle_hum$index),FM_Muscle_hum$index),]
  chr13_Muscle_hum<-rbind(chr13_Muscle_hum,b)
}
chr13_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr13_Muscle_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr13_Muscle_hum_one2one<-chr13_Muscle_hum[match(intersect(chr13_pig2hum_one2one,chr13_Muscle_hum$gene_id),chr13_Muscle_hum$gene_id),]

chr14_Muscle_overloci_one2one<-chr14_Muscle_overloci[match(intersect(chr14_Muscle_genes,one2one_pig),chr14_Muscle_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Muscle_chr14$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Muscle_chr14$hum_chr<-tmp$V1
Muscle_chr14$hum_loci<-tmp$V3
FM_Muscle_hum$index<-paste0(FM_Muscle_hum$chr,"_",FM_Muscle_hum$variant_pos)
Muscle_chr14$index<-paste0("chr",Muscle_chr14$hum_chr,"_",Muscle_chr14$hum_loci)
chr14_Muscle_hum<-NULL
for(i in 1:length(same_Muscle_chr14)){
  a<-Muscle_chr14[match(same_Muscle_chr14[i],Muscle_chr14$pig_pos),]
  b<-FM_Muscle_hum[match(intersect(a$index,FM_Muscle_hum$index),FM_Muscle_hum$index),]
  chr14_Muscle_hum<-rbind(chr14_Muscle_hum,b)
}
chr14_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr14_Muscle_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr14_Muscle_hum_one2one<-chr14_Muscle_hum[match(intersect(chr14_pig2hum_one2one,chr14_Muscle_hum$gene_id),chr14_Muscle_hum$gene_id),]

chr15_Muscle_overloci_one2one<-chr15_Muscle_overloci[match(intersect(chr15_Muscle_genes,one2one_pig),chr15_Muscle_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Muscle_chr15$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Muscle_chr15$hum_chr<-tmp$V1
Muscle_chr15$hum_loci<-tmp$V3
FM_Muscle_hum$index<-paste0(FM_Muscle_hum$chr,"_",FM_Muscle_hum$variant_pos)
Muscle_chr15$index<-paste0("chr",Muscle_chr15$hum_chr,"_",Muscle_chr15$hum_loci)
chr15_Muscle_hum<-NULL
for(i in 1:length(same_Muscle_chr15)){
  a<-Muscle_chr15[match(same_Muscle_chr15[i],Muscle_chr15$pig_pos),]
  b<-FM_Muscle_hum[match(intersect(a$index,FM_Muscle_hum$index),FM_Muscle_hum$index),]
  chr15_Muscle_hum<-rbind(chr15_Muscle_hum,b)
}
chr15_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr15_Muscle_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr15_Muscle_hum_one2one<-chr15_Muscle_hum[match(intersect(chr15_pig2hum_one2one,chr15_Muscle_hum$gene_id),chr15_Muscle_hum$gene_id),]

chr16_Muscle_overloci_one2one<-chr16_Muscle_overloci[match(intersect(chr16_Muscle_genes,one2one_pig),chr16_Muscle_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Muscle_chr16$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Muscle_chr16$hum_chr<-tmp$V1
Muscle_chr16$hum_loci<-tmp$V3
FM_Muscle_hum$index<-paste0(FM_Muscle_hum$chr,"_",FM_Muscle_hum$variant_pos)
Muscle_chr16$index<-paste0("chr",Muscle_chr16$hum_chr,"_",Muscle_chr16$hum_loci)
chr16_Muscle_hum<-NULL
for(i in 1:length(same_Muscle_chr16)){
  a<-Muscle_chr16[match(same_Muscle_chr16[i],Muscle_chr16$pig_pos),]
  b<-FM_Muscle_hum[match(intersect(a$index,FM_Muscle_hum$index),FM_Muscle_hum$index),]
  chr16_Muscle_hum<-rbind(chr16_Muscle_hum,b)
}
chr16_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr16_Muscle_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr16_Muscle_hum_one2one<-chr16_Muscle_hum[match(intersect(chr16_pig2hum_one2one,chr16_Muscle_hum$gene_id),chr16_Muscle_hum$gene_id),]

chr17_Muscle_overloci_one2one<-chr17_Muscle_overloci[match(intersect(chr17_Muscle_genes,one2one_pig),chr17_Muscle_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Muscle_chr17$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Muscle_chr17$hum_chr<-tmp$V1
Muscle_chr17$hum_loci<-tmp$V3
FM_Muscle_hum$index<-paste0(FM_Muscle_hum$chr,"_",FM_Muscle_hum$variant_pos)
Muscle_chr17$index<-paste0("chr",Muscle_chr17$hum_chr,"_",Muscle_chr17$hum_loci)
chr17_Muscle_hum<-NULL
for(i in 1:length(same_Muscle_chr17)){
  a<-Muscle_chr17[match(same_Muscle_chr17[i],Muscle_chr17$pig_pos),]
  b<-FM_Muscle_hum[match(intersect(a$index,FM_Muscle_hum$index),FM_Muscle_hum$index),]
  chr17_Muscle_hum<-rbind(chr17_Muscle_hum,b)
}
chr17_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr17_Muscle_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr17_Muscle_hum_one2one<-chr17_Muscle_hum[match(intersect(chr17_pig2hum_one2one,chr17_Muscle_hum$gene_id),chr17_Muscle_hum$gene_id),]

chr18_Muscle_overloci_one2one<-chr18_Muscle_overloci[match(intersect(chr18_Muscle_genes,one2one_pig),chr18_Muscle_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Muscle_chr18$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Muscle_chr18$hum_chr<-tmp$V1
Muscle_chr18$hum_loci<-tmp$V3
FM_Muscle_hum$index<-paste0(FM_Muscle_hum$chr,"_",FM_Muscle_hum$variant_pos)
Muscle_chr18$index<-paste0("chr",Muscle_chr18$hum_chr,"_",Muscle_chr18$hum_loci)
chr18_Muscle_hum<-NULL
for(i in 1:length(same_Muscle_chr18)){
  a<-Muscle_chr18[match(same_Muscle_chr18[i],Muscle_chr18$pig_pos),]
  b<-FM_Muscle_hum[match(intersect(a$index,FM_Muscle_hum$index),FM_Muscle_hum$index),]
  chr18_Muscle_hum<-rbind(chr18_Muscle_hum,b)
}
chr18_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr18_Muscle_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr18_Muscle_hum_one2one<-chr18_Muscle_hum[match(intersect(chr18_pig2hum_one2one,chr18_Muscle_hum$gene_id),chr18_Muscle_hum$gene_id),]

Muscle_one2one_SNP_hum<-rbind(chr1_Muscle_hum_one2one,chr2_Muscle_hum_one2one,chr3_Muscle_hum_one2one,chr4_Muscle_hum_one2one,chr5_Muscle_hum_one2one,
                          chr6_Muscle_hum_one2one,chr7_Muscle_hum_one2one,chr8_Muscle_hum_one2one,chr9_Muscle_hum_one2one,chr10_Muscle_hum_one2one,
                          chr11_Muscle_hum_one2one,chr12_Muscle_hum_one2one,chr13_Muscle_hum_one2one,chr14_Muscle_hum_one2one,chr15_Muscle_hum_one2one,
                          chr16_Muscle_hum_one2one,chr17_Muscle_hum_one2one,chr18_Muscle_hum_one2one)

chr1_Muscle_pig_one2one<-chr1_Muscle_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr1_Muscle_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr1_Muscle_overloci_one2one$phenotype_id),1:9]
chr2_Muscle_pig_one2one<-chr2_Muscle_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr2_Muscle_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr2_Muscle_overloci_one2one$phenotype_id),1:9]
chr3_Muscle_pig_one2one<-chr3_Muscle_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr3_Muscle_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr3_Muscle_overloci_one2one$phenotype_id),1:9]
chr4_Muscle_pig_one2one<-chr4_Muscle_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr4_Muscle_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr4_Muscle_overloci_one2one$phenotype_id),1:9]
chr5_Muscle_pig_one2one<-chr5_Muscle_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr5_Muscle_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr5_Muscle_overloci_one2one$phenotype_id),1:9]
chr6_Muscle_pig_one2one<-chr6_Muscle_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr6_Muscle_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr6_Muscle_overloci_one2one$phenotype_id),1:9]
chr7_Muscle_pig_one2one<-chr7_Muscle_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr7_Muscle_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr7_Muscle_overloci_one2one$phenotype_id),1:9]
chr8_Muscle_pig_one2one<-chr8_Muscle_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr8_Muscle_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr8_Muscle_overloci_one2one$phenotype_id),1:9]
chr9_Muscle_pig_one2one<-chr9_Muscle_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr9_Muscle_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr9_Muscle_overloci_one2one$phenotype_id),1:9]
chr10_Muscle_pig_one2one<-chr10_Muscle_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr10_Muscle_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr10_Muscle_overloci_one2one$phenotype_id),1:9]
chr11_Muscle_pig_one2one<-chr11_Muscle_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr11_Muscle_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr11_Muscle_overloci_one2one$phenotype_id),1:9]
chr12_Muscle_pig_one2one<-chr12_Muscle_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr12_Muscle_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr12_Muscle_overloci_one2one$phenotype_id),1:9]
chr13_Muscle_pig_one2one<-chr13_Muscle_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr13_Muscle_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr13_Muscle_overloci_one2one$phenotype_id),1:9]
chr14_Muscle_pig_one2one<-chr14_Muscle_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr14_Muscle_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr14_Muscle_overloci_one2one$phenotype_id),1:9]
chr15_Muscle_pig_one2one<-chr15_Muscle_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr15_Muscle_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr15_Muscle_overloci_one2one$phenotype_id),1:9]
chr16_Muscle_pig_one2one<-chr16_Muscle_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr16_Muscle_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr16_Muscle_overloci_one2one$phenotype_id),1:9]
chr17_Muscle_pig_one2one<-chr17_Muscle_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr17_Muscle_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr17_Muscle_overloci_one2one$phenotype_id),1:9]
chr18_Muscle_pig_one2one<-chr18_Muscle_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr18_Muscle_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr18_Muscle_overloci_one2one$phenotype_id),1:9]

Muscle_one2one_SNP_pig<-rbind(chr1_Muscle_pig_one2one,chr2_Muscle_pig_one2one,chr3_Muscle_pig_one2one,chr4_Muscle_pig_one2one,chr5_Muscle_pig_one2one,
                          chr6_Muscle_pig_one2one,chr7_Muscle_pig_one2one,chr8_Muscle_pig_one2one,chr9_Muscle_pig_one2one,chr10_Muscle_pig_one2one,
                          chr11_Muscle_pig_one2one,chr12_Muscle_pig_one2one,chr13_Muscle_pig_one2one,chr14_Muscle_pig_one2one,chr15_Muscle_pig_one2one,
                          chr16_Muscle_pig_one2one,chr17_Muscle_pig_one2one,chr18_Muscle_pig_one2one)

Muscle_SNP_sum<-array(NA,dim=c(nrow(Muscle_one2one_SNP_hum),2))
colnames(Muscle_SNP_sum)<-c("Human","Pig")
Muscle_SNP_sum<-as.data.frame(Muscle_SNP_sum)
Muscle_SNP_sum$Human<-Muscle_one2one_SNP_hum$slope / Muscle_one2one_SNP_hum$slope_se
Muscle_SNP_sum$Pig<-Muscle_one2one_SNP_pig$slope / Muscle_one2one_SNP_pig$slope_se
cor<-cor(abs(Muscle_SNP_sum$Human),abs(Muscle_SNP_sum$Pig))
p_val<-t.test(abs(Muscle_SNP_sum$Human),abs(Muscle_SNP_sum$Pig))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/Muscle_SNP_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(Muscle_SNP_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The number of eGenes expressed in common tissues",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

save(Muscle_one2one_SNP_hum,Muscle_one2one_SNP_pig,file="/Users/baizhonghao/Downloads/human-pig GTEx/Muscle_SNP.Rdata")

#Adipose_SNP_overlaploci#
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

eqtl_Adipose_chr1<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Adipose/Adipose.cis_qtl_pairs.1.txt"))
eqtl_Adipose_chr2<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Adipose/Adipose.cis_qtl_pairs.2.txt"))
eqtl_Adipose_chr3<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Adipose/Adipose.cis_qtl_pairs.3.txt"))
eqtl_Adipose_chr4<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Adipose/Adipose.cis_qtl_pairs.4.txt"))
eqtl_Adipose_chr5<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Adipose/Adipose.cis_qtl_pairs.5.txt"))
eqtl_Adipose_chr6<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Adipose/Adipose.cis_qtl_pairs.6.txt"))
eqtl_Adipose_chr7<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Adipose/Adipose.cis_qtl_pairs.7.txt"))
eqtl_Adipose_chr8<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Adipose/Adipose.cis_qtl_pairs.8.txt"))
eqtl_Adipose_chr9<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Adipose/Adipose.cis_qtl_pairs.9.txt"))
eqtl_Adipose_chr10<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Adipose/Adipose.cis_qtl_pairs.10.txt"))
eqtl_Adipose_chr11<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Adipose/Adipose.cis_qtl_pairs.11.txt"))
eqtl_Adipose_chr12<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Adipose/Adipose.cis_qtl_pairs.12.txt"))
eqtl_Adipose_chr13<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Adipose/Adipose.cis_qtl_pairs.13.txt"))
eqtl_Adipose_chr14<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Adipose/Adipose.cis_qtl_pairs.14.txt"))
eqtl_Adipose_chr15<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Adipose/Adipose.cis_qtl_pairs.15.txt"))
eqtl_Adipose_chr16<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Adipose/Adipose.cis_qtl_pairs.16.txt"))
eqtl_Adipose_chr17<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Adipose/Adipose.cis_qtl_pairs.17.txt"))
eqtl_Adipose_chr18<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Adipose/Adipose.cis_qtl_pairs.18.txt"))

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose_chr1$variant_id,split="_"))))
eqtl_Adipose_chr1$chr<-tmp$V1
eqtl_Adipose_chr1$loci<-tmp$V2

eqtl_Adipose_chr1$index<-paste0(eqtl_Adipose_chr1$chr,"-",eqtl_Adipose_chr1$loci)
Adipose_chr1$index<-paste0(Adipose_chr1$chr,"-",Adipose_chr1$pig_pos)
same_Adipose_chr1<-intersect(Adipose_chr1$pig_pos,eqtl_Adipose_chr1$loci)


tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose_chr2$variant_id,split="_"))))
eqtl_Adipose_chr2$chr<-tmp$V1
eqtl_Adipose_chr2$loci<-tmp$V2

eqtl_Adipose_chr2$index<-paste0(eqtl_Adipose_chr2$chr,"-",eqtl_Adipose_chr2$loci)
Adipose_chr2$index<-paste0(Adipose_chr2$chr,"-",Adipose_chr2$pig_pos)
same_Adipose_chr2<-intersect(Adipose_chr2$pig_pos,eqtl_Adipose_chr2$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose_chr3$variant_id,split="_"))))
eqtl_Adipose_chr3$chr<-tmp$V1
eqtl_Adipose_chr3$loci<-tmp$V2

eqtl_Adipose_chr3$index<-paste0(eqtl_Adipose_chr3$chr,"-",eqtl_Adipose_chr3$loci)
Adipose_chr3$index<-paste0(Adipose_chr3$chr,"-",Adipose_chr3$pig_pos)
same_Adipose_chr3<-intersect(Adipose_chr3$pig_pos,eqtl_Adipose_chr3$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose_chr4$variant_id,split="_"))))
eqtl_Adipose_chr4$chr<-tmp$V1
eqtl_Adipose_chr4$loci<-tmp$V2

eqtl_Adipose_chr4$index<-paste0(eqtl_Adipose_chr4$chr,"-",eqtl_Adipose_chr4$loci)
Adipose_chr4$index<-paste0(Adipose_chr4$chr,"-",Adipose_chr4$pig_pos)
same_Adipose_chr4<-intersect(Adipose_chr4$pig_pos,eqtl_Adipose_chr4$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose_chr5$variant_id,split="_"))))
eqtl_Adipose_chr5$loci<-tmp$V2
same_Adipose_chr5<-intersect(Adipose_chr5$pig_pos,eqtl_Adipose_chr5$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose_chr6$variant_id,split="_"))))
eqtl_Adipose_chr6$loci<-tmp$V2
same_Adipose_chr6<-intersect(Adipose_chr6$pig_pos,eqtl_Adipose_chr6$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose_chr7$variant_id,split="_"))))
eqtl_Adipose_chr7$loci<-tmp$V2
same_Adipose_chr7<-intersect(Adipose_chr7$pig_pos,eqtl_Adipose_chr7$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose_chr8$variant_id,split="_"))))
eqtl_Adipose_chr8$loci<-tmp$V2
same_Adipose_chr8<-intersect(Adipose_chr8$pig_pos,eqtl_Adipose_chr8$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose_chr9$variant_id,split="_"))))
eqtl_Adipose_chr9$loci<-tmp$V2
same_Adipose_chr9<-intersect(Adipose_chr9$pig_pos,eqtl_Adipose_chr9$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose_chr10$variant_id,split="_"))))
eqtl_Adipose_chr10$loci<-tmp$V2
same_Adipose_chr10<-intersect(Adipose_chr10$pig_pos,eqtl_Adipose_chr10$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose_chr11$variant_id,split="_"))))
eqtl_Adipose_chr11$loci<-tmp$V2
same_Adipose_chr11<-intersect(Adipose_chr11$pig_pos,eqtl_Adipose_chr11$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose_chr12$variant_id,split="_"))))
eqtl_Adipose_chr12$loci<-tmp$V2
same_Adipose_chr12<-intersect(Adipose_chr12$pig_pos,eqtl_Adipose_chr12$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose_chr13$variant_id,split="_"))))
eqtl_Adipose_chr13$loci<-tmp$V2
same_Adipose_chr13<-intersect(Adipose_chr13$pig_pos,eqtl_Adipose_chr13$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose_chr14$variant_id,split="_"))))
eqtl_Adipose_chr14$loci<-tmp$V2
same_Adipose_chr14<-intersect(Adipose_chr14$pig_pos,eqtl_Adipose_chr14$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose_chr15$variant_id,split="_"))))
eqtl_Adipose_chr15$loci<-tmp$V2
same_Adipose_chr15<-intersect(Adipose_chr15$pig_pos,eqtl_Adipose_chr15$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose_chr16$variant_id,split="_"))))
eqtl_Adipose_chr16$loci<-tmp$V2
same_Adipose_chr16<-intersect(Adipose_chr16$pig_pos,eqtl_Adipose_chr16$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose_chr17$variant_id,split="_"))))
eqtl_Adipose_chr17$loci<-tmp$V2
same_Adipose_chr17<-intersect(Adipose_chr17$pig_pos,eqtl_Adipose_chr17$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Adipose_chr18$variant_id,split="_"))))
eqtl_Adipose_chr18$loci<-tmp$V2
same_Adipose_chr18<-intersect(Adipose_chr18$pig_pos,eqtl_Adipose_chr18$loci)

chr1_Adipose_genes<-NULL
chr1_Adipose_overloci<-NULL
if(length(same_Adipose_chr1!=0)){
  for(i in 1:length(same_Adipose_chr1)){
    a<-eqtl_Adipose_chr1$phenotype_id[grep(same_Adipose_chr1[i],eqtl_Adipose_chr1$loci)]
    b<-eqtl_Adipose_chr1[grep(same_Adipose_chr1[i],eqtl_Adipose_chr1$loci),]
    chr1_Adipose_genes<-c(chr1_Adipose_genes,a)
    chr1_Adipose_overloci<-rbind(chr1_Adipose_overloci,b)
  }
}

if(length(same_Adipose_chr2!=0)){
  chr2_Adipose_genes<-NULL
  chr2_Adipose_overloci<-NULL
  for(i in 1:length(same_Adipose_chr2)){
    a<-eqtl_Adipose_chr2$phenotype_id[grep(same_Adipose_chr2[i],eqtl_Adipose_chr2$loci)]
    b<-eqtl_Adipose_chr2[grep(same_Adipose_chr2[i],eqtl_Adipose_chr2$loci),]
    chr2_Adipose_genes<-c(chr2_Adipose_genes,a)
    chr2_Adipose_overloci<-rbind(chr2_Adipose_overloci,b)
  }
}

if(length(same_Adipose_chr3!=0)){
  chr3_Adipose_genes<-NULL
  chr3_Adipose_overloci<-NULL
  for(i in 1:length(same_Adipose_chr3)){
    a<-eqtl_Adipose_chr3$phenotype_id[grep(same_Adipose_chr3[i],eqtl_Adipose_chr3$loci)]
    b<-eqtl_Adipose_chr3[grep(same_Adipose_chr3[i],eqtl_Adipose_chr3$loci),]
    chr3_Adipose_genes<-c(chr3_Adipose_genes,a)
    chr3_Adipose_overloci<-rbind(chr3_Adipose_overloci,b)
  }
}

if(length(same_Adipose_chr4!=0)){
  chr4_Adipose_genes<-NULL
  chr4_Adipose_overloci<-NULL
  for(i in 1:length(same_Adipose_chr4)){
    a<-eqtl_Adipose_chr4$phenotype_id[grep(same_Adipose_chr4[i],eqtl_Adipose_chr4$loci)]
    b<-eqtl_Adipose_chr4[grep(same_Adipose_chr4[i],eqtl_Adipose_chr4$loci),]
    chr4_Adipose_genes<-c(chr4_Adipose_genes,a)
    chr4_Adipose_overloci<-rbind(chr4_Adipose_overloci,b)
  }
}

if(length(same_Adipose_chr5!=0)){
  chr5_Adipose_genes<-NULL
  chr5_Adipose_overloci<-NULL
  for(i in 1:length(same_Adipose_chr5)){
    a<-eqtl_Adipose_chr5$phenotype_id[grep(same_Adipose_chr5[i],eqtl_Adipose_chr5$loci)]
    b<-eqtl_Adipose_chr5[grep(same_Adipose_chr5[i],eqtl_Adipose_chr5$loci),]
    chr5_Adipose_genes<-c(chr5_Adipose_genes,a)
    chr5_Adipose_overloci<-rbind(chr5_Adipose_overloci,b)
  }
}

if(length(same_Adipose_chr6!=0)){
  chr6_Adipose_genes<-NULL
  chr6_Adipose_overloci<-NULL
  for(i in 1:length(same_Adipose_chr6)){
    a<-eqtl_Adipose_chr6$phenotype_id[grep(same_Adipose_chr6[i],eqtl_Adipose_chr6$loci)]
    b<-eqtl_Adipose_chr6[grep(same_Adipose_chr6[i],eqtl_Adipose_chr6$loci),]
    chr6_Adipose_genes<-c(chr6_Adipose_genes,a)
    chr6_Adipose_overloci<-rbind(chr6_Adipose_overloci,b)
  }
}

if(length(same_Adipose_chr7!=0)){
  chr7_Adipose_genes<-NULL
  chr7_Adipose_overloci<-NULL
  for(i in 1:length(same_Adipose_chr7)){
    a<-eqtl_Adipose_chr7$phenotype_id[grep(same_Adipose_chr7[i],eqtl_Adipose_chr7$loci)]
    b<-eqtl_Adipose_chr7[grep(same_Adipose_chr7[i],eqtl_Adipose_chr7$loci),]
    chr7_Adipose_genes<-c(chr7_Adipose_genes,a)
    chr7_Adipose_overloci<-rbind(chr7_Adipose_overloci,b)
  }
}

if(length(same_Adipose_chr8!=0)){
  chr8_Adipose_genes<-NULL
  chr8_Adipose_overloci<-NULL
  for(i in 1:length(same_Adipose_chr8)){
    a<-eqtl_Adipose_chr8$phenotype_id[grep(same_Adipose_chr8[i],eqtl_Adipose_chr8$loci)]
    b<-eqtl_Adipose_chr8[grep(same_Adipose_chr8[i],eqtl_Adipose_chr8$loci),]
    chr8_Adipose_genes<-c(chr8_Adipose_genes,a)
    chr8_Adipose_overloci<-rbind(chr8_Adipose_overloci,b)
  }
}

if(length(same_Adipose_chr9!=0)){
  chr9_Adipose_genes<-NULL
  chr9_Adipose_overloci<-NULL
  for(i in 1:length(same_Adipose_chr9)){
    a<-eqtl_Adipose_chr9$phenotype_id[grep(same_Adipose_chr9[i],eqtl_Adipose_chr9$loci)]
    b<-eqtl_Adipose_chr9[grep(same_Adipose_chr9[i],eqtl_Adipose_chr9$loci),]
    chr9_Adipose_genes<-c(chr9_Adipose_genes,a)
    chr9_Adipose_overloci<-rbind(chr9_Adipose_overloci,b)
  }
}

if(length(same_Adipose_chr10!=0)){
  chr10_Adipose_genes<-NULL
  chr10_Adipose_overloci<-NULL
  for(i in 1:length(same_Adipose_chr10)){
    a<-eqtl_Adipose_chr10$phenotype_id[grep(same_Adipose_chr10[i],eqtl_Adipose_chr10$loci)]
    b<-eqtl_Adipose_chr10[grep(same_Adipose_chr10[i],eqtl_Adipose_chr10$loci),]
    chr10_Adipose_genes<-c(chr10_Adipose_genes,a)
    chr10_Adipose_overloci<-rbind(chr10_Adipose_overloci,b)
  }
}

if(length(same_Adipose_chr11!=0)){
  chr11_Adipose_genes<-NULL
  chr11_Adipose_overloci<-NULL
  for(i in 1:length(same_Adipose_chr11)){
    a<-eqtl_Adipose_chr11$phenotype_id[grep(same_Adipose_chr11[i],eqtl_Adipose_chr11$loci)]
    b<-eqtl_Adipose_chr11[grep(same_Adipose_chr11[i],eqtl_Adipose_chr11$loci),]
    chr11_Adipose_genes<-c(chr11_Adipose_genes,a)
    chr11_Adipose_overloci<-rbind(chr11_Adipose_overloci,b)
  }
}

chr12_Adipose_genes<-NULL
chr12_Adipose_overloci<-NULL
if(length(same_Adipose_chr12!=0)){
  for(i in 1:length(same_Adipose_chr12)){
    a<-eqtl_Adipose_chr12$phenotype_id[grep(same_Adipose_chr12[i],eqtl_Adipose_chr12$loci)]
    b<-eqtl_Adipose_chr12[grep(same_Adipose_chr12[i],eqtl_Adipose_chr12$loci),]
    chr12_Adipose_genes<-c(chr12_Adipose_genes,a)
    chr12_Adipose_overloci<-rbind(chr12_Adipose_overloci,b)
  }
}

chr13_Adipose_genes<-NULL
chr13_Adipose_overloci<-NULL
if(length(same_Adipose_chr13!=0)){
  for(i in 1:length(same_Adipose_chr13)){
    a<-eqtl_Adipose_chr13$phenotype_id[grep(same_Adipose_chr13[i],eqtl_Adipose_chr13$loci)]
    b<-eqtl_Adipose_chr13[grep(same_Adipose_chr13[i],eqtl_Adipose_chr13$loci),]
    chr13_Adipose_genes<-c(chr13_Adipose_genes,a)
    chr13_Adipose_overloci<-rbind(chr13_Adipose_overloci,b)
  }
}

chr14_Adipose_genes<-NULL
chr14_Adipose_overloci<-NULL
if(length(same_Adipose_chr14!=0)){
  for(i in 1:length(same_Adipose_chr14)){
    a<-eqtl_Adipose_chr14$phenotype_id[grep(same_Adipose_chr14[i],eqtl_Adipose_chr14$loci)]
    b<-eqtl_Adipose_chr14[grep(same_Adipose_chr14[i],eqtl_Adipose_chr14$loci),]
    chr14_Adipose_genes<-c(chr14_Adipose_genes,a)
    chr14_Adipose_overloci<-rbind(chr14_Adipose_overloci,b)
  }
}

chr15_Adipose_genes<-NULL
chr15_Adipose_overloci<-NULL
if(length(same_Adipose_chr15!=0)){
  for(i in 1:length(same_Adipose_chr15)){
    a<-eqtl_Adipose_chr15$phenotype_id[grep(same_Adipose_chr15[i],eqtl_Adipose_chr15$loci)]
    b<-eqtl_Adipose_chr15[grep(same_Adipose_chr15[i],eqtl_Adipose_chr15$loci),]
    chr15_Adipose_genes<-c(chr15_Adipose_genes,a)
    chr15_Adipose_overloci<-rbind(chr15_Adipose_overloci,b)
  }
}

chr16_Adipose_genes<-NULL
chr16_Adipose_overloci<-NULL
if(length(same_Adipose_chr16!=0)){
  for(i in 1:length(same_Adipose_chr16)){
    a<-eqtl_Adipose_chr16$phenotype_id[grep(same_Adipose_chr16[i],eqtl_Adipose_chr16$loci)]
    b<-eqtl_Adipose_chr16[grep(same_Adipose_chr16[i],eqtl_Adipose_chr16$loci),]
    chr16_Adipose_genes<-c(chr16_Adipose_genes,a)
    chr16_Adipose_overloci<-rbind(chr16_Adipose_overloci,b)
  }
}

chr17_Adipose_genes<-NULL
chr17_Adipose_overloci<-NULL
if(length(same_Adipose_chr17!=0)){
  for(i in 1:length(same_Adipose_chr17)){
    a<-eqtl_Adipose_chr17$phenotype_id[grep(same_Adipose_chr17[i],eqtl_Adipose_chr17$loci)]
    b<-eqtl_Adipose_chr17[grep(same_Adipose_chr17[i],eqtl_Adipose_chr17$loci),]
    chr17_Adipose_genes<-c(chr17_Adipose_genes,a)
    chr17_Adipose_overloci<-rbind(chr17_Adipose_overloci,b)
  }
}

chr18_Adipose_genes<-NULL
chr18_Adipose_overloci<-NULL
if(length(same_Adipose_chr18!=0)){
  for(i in 1:length(same_Adipose_chr18)){
    a<-eqtl_Adipose_chr18$phenotype_id[grep(same_Adipose_chr18[i],eqtl_Adipose_chr18$loci)]
    b<-eqtl_Adipose_chr18[grep(same_Adipose_chr18[i],eqtl_Adipose_chr18$loci),]
    chr18_Adipose_genes<-c(chr18_Adipose_genes,a)
    chr18_Adipose_overloci<-rbind(chr18_Adipose_overloci,b)
  }
}
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
one2one_pig<-annotation$`Pig gene stable ID`[which(annotation$`Pig homology type`=="ortholog_one2one")]

chr1_Adipose_overloci_one2one<-chr1_Adipose_overloci[match(intersect(chr1_Adipose_genes,one2one_pig),chr1_Adipose_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Adipose_chr1$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Adipose_chr1$hum_chr<-tmp$V1
Adipose_chr1$hum_loci<-tmp$V3
FM_Adipose_hum$index<-paste0(FM_Adipose_hum$chr,"_",FM_Adipose_hum$variant_pos)
Adipose_chr1$index<-paste0("chr",Adipose_chr1$hum_chr,"_",Adipose_chr1$hum_loci)
chr1_Adipose_hum<-NULL
for(i in 1:length(same_Adipose_chr1)){
  a<-Adipose_chr1[match(same_Adipose_chr1[i],Adipose_chr1$pig_pos),]
  b<-FM_Adipose_hum[match(intersect(a$index,FM_Adipose_hum$index),FM_Adipose_hum$index),]
  chr1_Adipose_hum<-rbind(chr1_Adipose_hum,b)
}
chr1_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr1_Adipose_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr1_Adipose_hum_one2one<-chr1_Adipose_hum[match(intersect(chr1_pig2hum_one2one,chr1_Adipose_hum$gene_id),chr1_Adipose_hum$gene_id),]

chr2_Adipose_overloci_one2one<-chr2_Adipose_overloci[match(intersect(chr2_Adipose_genes,one2one_pig),chr2_Adipose_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Adipose_chr2$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Adipose_chr2$hum_chr<-tmp$V1
Adipose_chr2$hum_loci<-tmp$V3
FM_Adipose_hum$index<-paste0(FM_Adipose_hum$chr,"_",FM_Adipose_hum$variant_pos)
Adipose_chr2$index<-paste0("chr",Adipose_chr2$hum_chr,"_",Adipose_chr2$hum_loci)
chr2_Adipose_hum<-NULL
for(i in 1:length(same_Adipose_chr2)){
  a<-Adipose_chr2[match(same_Adipose_chr2[i],Adipose_chr2$pig_pos),]
  b<-FM_Adipose_hum[match(intersect(a$index,FM_Adipose_hum$index),FM_Adipose_hum$index),]
  chr2_Adipose_hum<-rbind(chr2_Adipose_hum,b)
}
chr2_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr2_Adipose_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr2_Adipose_hum_one2one<-chr2_Adipose_hum[match(intersect(chr2_pig2hum_one2one,chr2_Adipose_hum$gene_id),chr2_Adipose_hum$gene_id),]

chr3_Adipose_overloci_one2one<-chr3_Adipose_overloci[match(intersect(chr3_Adipose_genes,one2one_pig),chr3_Adipose_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Adipose_chr3$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Adipose_chr3$hum_chr<-tmp$V1
Adipose_chr3$hum_loci<-tmp$V3
FM_Adipose_hum$index<-paste0(FM_Adipose_hum$chr,"_",FM_Adipose_hum$variant_pos)
Adipose_chr3$index<-paste0("chr",Adipose_chr3$hum_chr,"_",Adipose_chr3$hum_loci)
chr3_Adipose_hum<-NULL
for(i in 1:length(same_Adipose_chr3)){
  a<-Adipose_chr3[match(same_Adipose_chr3[i],Adipose_chr3$pig_pos),]
  b<-FM_Adipose_hum[match(intersect(a$index,FM_Adipose_hum$index),FM_Adipose_hum$index),]
  chr3_Adipose_hum<-rbind(chr3_Adipose_hum,b)
}
chr3_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr3_Adipose_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr3_Adipose_hum_one2one<-chr3_Adipose_hum[match(intersect(chr3_pig2hum_one2one,chr3_Adipose_hum$gene_id),chr3_Adipose_hum$gene_id),]

chr4_Adipose_overloci_one2one<-chr4_Adipose_overloci[match(intersect(chr4_Adipose_genes,one2one_pig),chr4_Adipose_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Adipose_chr4$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Adipose_chr4$hum_chr<-tmp$V1
Adipose_chr4$hum_loci<-tmp$V3
FM_Adipose_hum$index<-paste0(FM_Adipose_hum$chr,"_",FM_Adipose_hum$variant_pos)
Adipose_chr4$index<-paste0("chr",Adipose_chr4$hum_chr,"_",Adipose_chr4$hum_loci)
chr4_Adipose_hum<-NULL
for(i in 1:length(same_Adipose_chr4)){
  a<-Adipose_chr4[match(same_Adipose_chr4[i],Adipose_chr4$pig_pos),]
  b<-FM_Adipose_hum[match(intersect(a$index,FM_Adipose_hum$index),FM_Adipose_hum$index),]
  chr4_Adipose_hum<-rbind(chr4_Adipose_hum,b)
}
chr4_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr4_Adipose_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr4_Adipose_hum_one2one<-chr4_Adipose_hum[match(intersect(chr4_pig2hum_one2one,chr4_Adipose_hum$gene_id),chr4_Adipose_hum$gene_id),]

chr5_Adipose_overloci_one2one<-chr5_Adipose_overloci[match(intersect(chr5_Adipose_genes,one2one_pig),chr5_Adipose_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Adipose_chr5$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Adipose_chr5$hum_chr<-tmp$V1
Adipose_chr5$hum_loci<-tmp$V3
FM_Adipose_hum$index<-paste0(FM_Adipose_hum$chr,"_",FM_Adipose_hum$variant_pos)
Adipose_chr5$index<-paste0("chr",Adipose_chr5$hum_chr,"_",Adipose_chr5$hum_loci)
chr5_Adipose_hum<-NULL
for(i in 1:length(same_Adipose_chr5)){
  a<-Adipose_chr5[match(same_Adipose_chr5[i],Adipose_chr5$pig_pos),]
  b<-FM_Adipose_hum[match(intersect(a$index,FM_Adipose_hum$index),FM_Adipose_hum$index),]
  chr5_Adipose_hum<-rbind(chr5_Adipose_hum,b)
}
chr5_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr5_Adipose_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr5_Adipose_hum_one2one<-chr5_Adipose_hum[match(intersect(chr5_pig2hum_one2one,chr5_Adipose_hum$gene_id),chr5_Adipose_hum$gene_id),]

chr6_Adipose_overloci_one2one<-chr6_Adipose_overloci[match(intersect(chr6_Adipose_genes,one2one_pig),chr6_Adipose_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Adipose_chr6$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Adipose_chr6$hum_chr<-tmp$V1
Adipose_chr6$hum_loci<-tmp$V3
FM_Adipose_hum$index<-paste0(FM_Adipose_hum$chr,"_",FM_Adipose_hum$variant_pos)
Adipose_chr6$index<-paste0("chr",Adipose_chr6$hum_chr,"_",Adipose_chr6$hum_loci)
chr6_Adipose_hum<-NULL
for(i in 1:length(same_Adipose_chr6)){
  a<-Adipose_chr6[match(same_Adipose_chr6[i],Adipose_chr6$pig_pos),]
  b<-FM_Adipose_hum[match(intersect(a$index,FM_Adipose_hum$index),FM_Adipose_hum$index),]
  chr6_Adipose_hum<-rbind(chr6_Adipose_hum,b)
}
chr6_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr6_Adipose_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr6_Adipose_hum_one2one<-chr6_Adipose_hum[match(intersect(chr6_pig2hum_one2one,chr6_Adipose_hum$gene_id),chr6_Adipose_hum$gene_id),]

chr7_Adipose_overloci_one2one<-chr7_Adipose_overloci[match(intersect(chr7_Adipose_genes,one2one_pig),chr7_Adipose_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Adipose_chr7$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Adipose_chr7$hum_chr<-tmp$V1
Adipose_chr7$hum_loci<-tmp$V3
FM_Adipose_hum$index<-paste0(FM_Adipose_hum$chr,"_",FM_Adipose_hum$variant_pos)
Adipose_chr7$index<-paste0("chr",Adipose_chr7$hum_chr,"_",Adipose_chr7$hum_loci)
chr7_Adipose_hum<-NULL
for(i in 1:length(same_Adipose_chr7)){
  a<-Adipose_chr7[match(same_Adipose_chr7[i],Adipose_chr7$pig_pos),]
  b<-FM_Adipose_hum[match(intersect(a$index,FM_Adipose_hum$index),FM_Adipose_hum$index),]
  chr7_Adipose_hum<-rbind(chr7_Adipose_hum,b)
}
chr7_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr7_Adipose_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr7_Adipose_hum_one2one<-chr7_Adipose_hum[match(intersect(chr7_pig2hum_one2one,chr7_Adipose_hum$gene_id),chr7_Adipose_hum$gene_id),]

chr8_Adipose_overloci_one2one<-chr8_Adipose_overloci[match(intersect(chr8_Adipose_genes,one2one_pig),chr8_Adipose_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Adipose_chr8$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Adipose_chr8$hum_chr<-tmp$V1
Adipose_chr8$hum_loci<-tmp$V3
FM_Adipose_hum$index<-paste0(FM_Adipose_hum$chr,"_",FM_Adipose_hum$variant_pos)
Adipose_chr8$index<-paste0("chr",Adipose_chr8$hum_chr,"_",Adipose_chr8$hum_loci)
chr8_Adipose_hum<-NULL
for(i in 1:length(same_Adipose_chr8)){
  a<-Adipose_chr8[match(same_Adipose_chr8[i],Adipose_chr8$pig_pos),]
  b<-FM_Adipose_hum[match(intersect(a$index,FM_Adipose_hum$index),FM_Adipose_hum$index),]
  chr8_Adipose_hum<-rbind(chr8_Adipose_hum,b)
}
chr8_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr8_Adipose_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr8_Adipose_hum_one2one<-chr8_Adipose_hum[match(intersect(chr8_pig2hum_one2one,chr8_Adipose_hum$gene_id),chr8_Adipose_hum$gene_id),]

chr9_Adipose_overloci_one2one<-chr9_Adipose_overloci[match(intersect(chr9_Adipose_genes,one2one_pig),chr9_Adipose_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Adipose_chr9$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Adipose_chr9$hum_chr<-tmp$V1
Adipose_chr9$hum_loci<-tmp$V3
FM_Adipose_hum$index<-paste0(FM_Adipose_hum$chr,"_",FM_Adipose_hum$variant_pos)
Adipose_chr9$index<-paste0("chr",Adipose_chr9$hum_chr,"_",Adipose_chr9$hum_loci)
chr9_Adipose_hum<-NULL
for(i in 1:length(same_Adipose_chr9)){
  a<-Adipose_chr9[match(same_Adipose_chr9[i],Adipose_chr9$pig_pos),]
  b<-FM_Adipose_hum[match(intersect(a$index,FM_Adipose_hum$index),FM_Adipose_hum$index),]
  chr9_Adipose_hum<-rbind(chr9_Adipose_hum,b)
}
chr9_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr9_Adipose_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr9_Adipose_hum_one2one<-chr9_Adipose_hum[match(intersect(chr9_pig2hum_one2one,chr9_Adipose_hum$gene_id),chr9_Adipose_hum$gene_id),]

chr10_Adipose_overloci_one2one<-chr10_Adipose_overloci[match(intersect(chr10_Adipose_genes,one2one_pig),chr10_Adipose_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Adipose_chr10$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Adipose_chr10$hum_chr<-tmp$V1
Adipose_chr10$hum_loci<-tmp$V3
FM_Adipose_hum$index<-paste0(FM_Adipose_hum$chr,"_",FM_Adipose_hum$variant_pos)
Adipose_chr10$index<-paste0("chr",Adipose_chr10$hum_chr,"_",Adipose_chr10$hum_loci)
chr10_Adipose_hum<-NULL
for(i in 1:length(same_Adipose_chr10)){
  a<-Adipose_chr10[match(same_Adipose_chr10[i],Adipose_chr10$pig_pos),]
  b<-FM_Adipose_hum[match(intersect(a$index,FM_Adipose_hum$index),FM_Adipose_hum$index),]
  chr10_Adipose_hum<-rbind(chr10_Adipose_hum,b)
}
chr10_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr10_Adipose_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr10_Adipose_hum_one2one<-chr10_Adipose_hum[match(intersect(chr10_pig2hum_one2one,chr10_Adipose_hum$gene_id),chr10_Adipose_hum$gene_id),]

chr11_Adipose_overloci_one2one<-chr11_Adipose_overloci[match(intersect(chr11_Adipose_genes,one2one_pig),chr11_Adipose_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Adipose_chr11$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Adipose_chr11$hum_chr<-tmp$V1
Adipose_chr11$hum_loci<-tmp$V3
FM_Adipose_hum$index<-paste0(FM_Adipose_hum$chr,"_",FM_Adipose_hum$variant_pos)
Adipose_chr11$index<-paste0("chr",Adipose_chr11$hum_chr,"_",Adipose_chr11$hum_loci)
chr11_Adipose_hum<-NULL
for(i in 1:length(same_Adipose_chr11)){
  a<-Adipose_chr11[match(same_Adipose_chr11[i],Adipose_chr11$pig_pos),]
  b<-FM_Adipose_hum[match(intersect(a$index,FM_Adipose_hum$index),FM_Adipose_hum$index),]
  chr11_Adipose_hum<-rbind(chr11_Adipose_hum,b)
}
chr11_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr11_Adipose_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr11_Adipose_hum_one2one<-chr11_Adipose_hum[match(intersect(chr11_pig2hum_one2one,chr11_Adipose_hum$gene_id),chr11_Adipose_hum$gene_id),]

chr12_Adipose_overloci_one2one<-chr12_Adipose_overloci[match(intersect(chr12_Adipose_genes,one2one_pig),chr12_Adipose_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Adipose_chr12$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Adipose_chr12$hum_chr<-tmp$V1
Adipose_chr12$hum_loci<-tmp$V3
FM_Adipose_hum$index<-paste0(FM_Adipose_hum$chr,"_",FM_Adipose_hum$variant_pos)
Adipose_chr12$index<-paste0("chr",Adipose_chr12$hum_chr,"_",Adipose_chr12$hum_loci)
chr12_Adipose_hum<-NULL
for(i in 1:length(same_Adipose_chr12)){
  a<-Adipose_chr12[match(same_Adipose_chr12[i],Adipose_chr12$pig_pos),]
  b<-FM_Adipose_hum[match(intersect(a$index,FM_Adipose_hum$index),FM_Adipose_hum$index),]
  chr12_Adipose_hum<-rbind(chr12_Adipose_hum,b)
}
chr12_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr12_Adipose_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr12_Adipose_hum_one2one<-chr12_Adipose_hum[match(intersect(chr12_pig2hum_one2one,chr12_Adipose_hum$gene_id),chr12_Adipose_hum$gene_id),]

chr13_Adipose_overloci_one2one<-chr13_Adipose_overloci[match(intersect(chr13_Adipose_genes,one2one_pig),chr13_Adipose_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Adipose_chr13$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Adipose_chr13$hum_chr<-tmp$V1
Adipose_chr13$hum_loci<-tmp$V3
FM_Adipose_hum$index<-paste0(FM_Adipose_hum$chr,"_",FM_Adipose_hum$variant_pos)
Adipose_chr13$index<-paste0("chr",Adipose_chr13$hum_chr,"_",Adipose_chr13$hum_loci)
chr13_Adipose_hum<-NULL
for(i in 1:length(same_Adipose_chr13)){
  a<-Adipose_chr13[match(same_Adipose_chr13[i],Adipose_chr13$pig_pos),]
  b<-FM_Adipose_hum[match(intersect(a$index,FM_Adipose_hum$index),FM_Adipose_hum$index),]
  chr13_Adipose_hum<-rbind(chr13_Adipose_hum,b)
}
chr13_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr13_Adipose_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr13_Adipose_hum_one2one<-chr13_Adipose_hum[match(intersect(chr13_pig2hum_one2one,chr13_Adipose_hum$gene_id),chr13_Adipose_hum$gene_id),]

chr14_Adipose_overloci_one2one<-chr14_Adipose_overloci[match(intersect(chr14_Adipose_genes,one2one_pig),chr14_Adipose_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Adipose_chr14$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Adipose_chr14$hum_chr<-tmp$V1
Adipose_chr14$hum_loci<-tmp$V3
FM_Adipose_hum$index<-paste0(FM_Adipose_hum$chr,"_",FM_Adipose_hum$variant_pos)
Adipose_chr14$index<-paste0("chr",Adipose_chr14$hum_chr,"_",Adipose_chr14$hum_loci)
chr14_Adipose_hum<-NULL
for(i in 1:length(same_Adipose_chr14)){
  a<-Adipose_chr14[match(same_Adipose_chr14[i],Adipose_chr14$pig_pos),]
  b<-FM_Adipose_hum[match(intersect(a$index,FM_Adipose_hum$index),FM_Adipose_hum$index),]
  chr14_Adipose_hum<-rbind(chr14_Adipose_hum,b)
}
chr14_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr14_Adipose_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr14_Adipose_hum_one2one<-chr14_Adipose_hum[match(intersect(chr14_pig2hum_one2one,chr14_Adipose_hum$gene_id),chr14_Adipose_hum$gene_id),]

chr15_Adipose_overloci_one2one<-chr15_Adipose_overloci[match(intersect(chr15_Adipose_genes,one2one_pig),chr15_Adipose_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Adipose_chr15$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Adipose_chr15$hum_chr<-tmp$V1
Adipose_chr15$hum_loci<-tmp$V3
FM_Adipose_hum$index<-paste0(FM_Adipose_hum$chr,"_",FM_Adipose_hum$variant_pos)
Adipose_chr15$index<-paste0("chr",Adipose_chr15$hum_chr,"_",Adipose_chr15$hum_loci)
chr15_Adipose_hum<-NULL
for(i in 1:length(same_Adipose_chr15)){
  a<-Adipose_chr15[match(same_Adipose_chr15[i],Adipose_chr15$pig_pos),]
  b<-FM_Adipose_hum[match(intersect(a$index,FM_Adipose_hum$index),FM_Adipose_hum$index),]
  chr15_Adipose_hum<-rbind(chr15_Adipose_hum,b)
}
chr15_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr15_Adipose_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr15_Adipose_hum_one2one<-chr15_Adipose_hum[match(intersect(chr15_pig2hum_one2one,chr15_Adipose_hum$gene_id),chr15_Adipose_hum$gene_id),]

chr16_Adipose_overloci_one2one<-chr16_Adipose_overloci[match(intersect(chr16_Adipose_genes,one2one_pig),chr16_Adipose_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Adipose_chr16$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Adipose_chr16$hum_chr<-tmp$V1
Adipose_chr16$hum_loci<-tmp$V3
FM_Adipose_hum$index<-paste0(FM_Adipose_hum$chr,"_",FM_Adipose_hum$variant_pos)
Adipose_chr16$index<-paste0("chr",Adipose_chr16$hum_chr,"_",Adipose_chr16$hum_loci)
chr16_Adipose_hum<-NULL
for(i in 1:length(same_Adipose_chr16)){
  a<-Adipose_chr16[match(same_Adipose_chr16[i],Adipose_chr16$pig_pos),]
  b<-FM_Adipose_hum[match(intersect(a$index,FM_Adipose_hum$index),FM_Adipose_hum$index),]
  chr16_Adipose_hum<-rbind(chr16_Adipose_hum,b)
}
chr16_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr16_Adipose_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr16_Adipose_hum_one2one<-chr16_Adipose_hum[match(intersect(chr16_pig2hum_one2one,chr16_Adipose_hum$gene_id),chr16_Adipose_hum$gene_id),]

chr17_Adipose_overloci_one2one<-chr17_Adipose_overloci[match(intersect(chr17_Adipose_genes,one2one_pig),chr17_Adipose_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Adipose_chr17$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Adipose_chr17$hum_chr<-tmp$V1
Adipose_chr17$hum_loci<-tmp$V3
FM_Adipose_hum$index<-paste0(FM_Adipose_hum$chr,"_",FM_Adipose_hum$variant_pos)
Adipose_chr17$index<-paste0("chr",Adipose_chr17$hum_chr,"_",Adipose_chr17$hum_loci)
chr17_Adipose_hum<-NULL
for(i in 1:length(same_Adipose_chr17)){
  a<-Adipose_chr17[match(same_Adipose_chr17[i],Adipose_chr17$pig_pos),]
  b<-FM_Adipose_hum[match(intersect(a$index,FM_Adipose_hum$index),FM_Adipose_hum$index),]
  chr17_Adipose_hum<-rbind(chr17_Adipose_hum,b)
}
chr17_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr17_Adipose_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr17_Adipose_hum_one2one<-chr17_Adipose_hum[match(intersect(chr17_pig2hum_one2one,chr17_Adipose_hum$gene_id),chr17_Adipose_hum$gene_id),]

chr18_Adipose_overloci_one2one<-chr18_Adipose_overloci[match(intersect(chr18_Adipose_genes,one2one_pig),chr18_Adipose_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Adipose_chr18$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Adipose_chr18$hum_chr<-tmp$V1
Adipose_chr18$hum_loci<-tmp$V3
FM_Adipose_hum$index<-paste0(FM_Adipose_hum$chr,"_",FM_Adipose_hum$variant_pos)
Adipose_chr18$index<-paste0("chr",Adipose_chr18$hum_chr,"_",Adipose_chr18$hum_loci)
chr18_Adipose_hum<-NULL
for(i in 1:length(same_Adipose_chr18)){
  a<-Adipose_chr18[match(same_Adipose_chr18[i],Adipose_chr18$pig_pos),]
  b<-FM_Adipose_hum[match(intersect(a$index,FM_Adipose_hum$index),FM_Adipose_hum$index),]
  chr18_Adipose_hum<-rbind(chr18_Adipose_hum,b)
}
chr18_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr18_Adipose_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr18_Adipose_hum_one2one<-chr18_Adipose_hum[match(intersect(chr18_pig2hum_one2one,chr18_Adipose_hum$gene_id),chr18_Adipose_hum$gene_id),]

Adipose_one2one_SNP_hum<-rbind(chr1_Adipose_hum_one2one,chr2_Adipose_hum_one2one,chr3_Adipose_hum_one2one,chr4_Adipose_hum_one2one,chr5_Adipose_hum_one2one,
                               chr6_Adipose_hum_one2one,chr7_Adipose_hum_one2one,chr8_Adipose_hum_one2one,chr9_Adipose_hum_one2one,chr10_Adipose_hum_one2one,
                               chr11_Adipose_hum_one2one,chr12_Adipose_hum_one2one,chr13_Adipose_hum_one2one,chr14_Adipose_hum_one2one,chr15_Adipose_hum_one2one,
                               chr16_Adipose_hum_one2one,chr17_Adipose_hum_one2one,chr18_Adipose_hum_one2one)

chr1_Adipose_pig_one2one<-chr1_Adipose_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr1_Adipose_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr1_Adipose_overloci_one2one$phenotype_id),1:9]
chr2_Adipose_pig_one2one<-chr2_Adipose_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr2_Adipose_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr2_Adipose_overloci_one2one$phenotype_id),1:9]
chr3_Adipose_pig_one2one<-chr3_Adipose_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr3_Adipose_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr3_Adipose_overloci_one2one$phenotype_id),1:9]
chr4_Adipose_pig_one2one<-chr4_Adipose_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr4_Adipose_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr4_Adipose_overloci_one2one$phenotype_id),1:9]
chr5_Adipose_pig_one2one<-chr5_Adipose_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr5_Adipose_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr5_Adipose_overloci_one2one$phenotype_id),1:9]
chr6_Adipose_pig_one2one<-chr6_Adipose_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr6_Adipose_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr6_Adipose_overloci_one2one$phenotype_id),1:9]
chr7_Adipose_pig_one2one<-chr7_Adipose_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr7_Adipose_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr7_Adipose_overloci_one2one$phenotype_id),1:9]
chr8_Adipose_pig_one2one<-chr8_Adipose_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr8_Adipose_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr8_Adipose_overloci_one2one$phenotype_id),1:9]
chr9_Adipose_pig_one2one<-chr9_Adipose_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr9_Adipose_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr9_Adipose_overloci_one2one$phenotype_id),1:9]
chr10_Adipose_pig_one2one<-chr10_Adipose_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr10_Adipose_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr10_Adipose_overloci_one2one$phenotype_id),1:9]
chr11_Adipose_pig_one2one<-chr11_Adipose_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr11_Adipose_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr11_Adipose_overloci_one2one$phenotype_id),1:9]
chr12_Adipose_pig_one2one<-chr12_Adipose_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr12_Adipose_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr12_Adipose_overloci_one2one$phenotype_id),1:9]
chr13_Adipose_pig_one2one<-chr13_Adipose_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr13_Adipose_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr13_Adipose_overloci_one2one$phenotype_id),1:9]
chr14_Adipose_pig_one2one<-chr14_Adipose_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr14_Adipose_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr14_Adipose_overloci_one2one$phenotype_id),1:9]
chr15_Adipose_pig_one2one<-chr15_Adipose_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr15_Adipose_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr15_Adipose_overloci_one2one$phenotype_id),1:9]
chr16_Adipose_pig_one2one<-chr16_Adipose_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr16_Adipose_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr16_Adipose_overloci_one2one$phenotype_id),1:9]
chr17_Adipose_pig_one2one<-chr17_Adipose_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr17_Adipose_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr17_Adipose_overloci_one2one$phenotype_id),1:9]
chr18_Adipose_pig_one2one<-chr18_Adipose_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr18_Adipose_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr18_Adipose_overloci_one2one$phenotype_id),1:9]

Adipose_one2one_SNP_pig<-rbind(chr1_Adipose_pig_one2one,chr2_Adipose_pig_one2one,chr3_Adipose_pig_one2one,chr4_Adipose_pig_one2one,chr5_Adipose_pig_one2one,
                               chr6_Adipose_pig_one2one,chr7_Adipose_pig_one2one,chr8_Adipose_pig_one2one,chr9_Adipose_pig_one2one,chr10_Adipose_pig_one2one,
                               chr11_Adipose_pig_one2one,chr12_Adipose_pig_one2one,chr13_Adipose_pig_one2one,chr14_Adipose_pig_one2one,chr15_Adipose_pig_one2one,
                               chr16_Adipose_pig_one2one,chr17_Adipose_pig_one2one,chr18_Adipose_pig_one2one)

Adipose_SNP_sum<-array(NA,dim=c(nrow(Adipose_one2one_SNP_hum),2))
colnames(Adipose_SNP_sum)<-c("Human","Pig")
Adipose_SNP_sum<-as.data.frame(Adipose_SNP_sum)
Adipose_SNP_sum$Human<-Adipose_one2one_SNP_hum$slope / Adipose_one2one_SNP_hum$slope_se
Adipose_SNP_sum$Pig<-Adipose_one2one_SNP_pig$slope / Adipose_one2one_SNP_pig$slope_se
cor<-cor(abs(Adipose_SNP_sum$Human),abs(Adipose_SNP_sum$Pig))
p_val<-t.test(abs(Adipose_SNP_sum$Human),abs(Adipose_SNP_sum$Pig))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/Adipose_SNP_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(Adipose_SNP_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The number of eGenes expressed in common tissues",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

save(Adipose_one2one_SNP_hum,Adipose_one2one_SNP_pig,file="/Users/baizhonghao/Downloads/human-pig GTEx/Adipose_SNP.Rdata")

#Artery_SNP_overlaploci#
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

eqtl_Artery_chr1<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Artery/Artery.cis_qtl_pairs.1.txt"))
eqtl_Artery_chr2<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Artery/Artery.cis_qtl_pairs.2.txt"))
eqtl_Artery_chr3<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Artery/Artery.cis_qtl_pairs.3.txt"))
eqtl_Artery_chr4<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Artery/Artery.cis_qtl_pairs.4.txt"))
eqtl_Artery_chr5<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Artery/Artery.cis_qtl_pairs.5.txt"))
eqtl_Artery_chr6<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Artery/Artery.cis_qtl_pairs.6.txt"))
eqtl_Artery_chr7<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Artery/Artery.cis_qtl_pairs.7.txt"))
eqtl_Artery_chr8<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Artery/Artery.cis_qtl_pairs.8.txt"))
eqtl_Artery_chr9<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Artery/Artery.cis_qtl_pairs.9.txt"))
eqtl_Artery_chr10<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Artery/Artery.cis_qtl_pairs.10.txt"))
eqtl_Artery_chr11<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Artery/Artery.cis_qtl_pairs.11.txt"))
eqtl_Artery_chr12<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Artery/Artery.cis_qtl_pairs.12.txt"))
eqtl_Artery_chr13<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Artery/Artery.cis_qtl_pairs.13.txt"))
eqtl_Artery_chr14<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Artery/Artery.cis_qtl_pairs.14.txt"))
eqtl_Artery_chr15<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Artery/Artery.cis_qtl_pairs.15.txt"))
eqtl_Artery_chr16<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Artery/Artery.cis_qtl_pairs.16.txt"))
eqtl_Artery_chr17<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Artery/Artery.cis_qtl_pairs.17.txt"))
eqtl_Artery_chr18<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Artery/Artery.cis_qtl_pairs.18.txt"))

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery_chr1$variant_id,split="_"))))
eqtl_Artery_chr1$chr<-tmp$V1
eqtl_Artery_chr1$loci<-tmp$V2

eqtl_Artery_chr1$index<-paste0(eqtl_Artery_chr1$chr,"-",eqtl_Artery_chr1$loci)
Artery_chr1$index<-paste0(Artery_chr1$chr,"-",Artery_chr1$pig_pos)
same_Artery_chr1<-intersect(Artery_chr1$pig_pos,eqtl_Artery_chr1$loci)


tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery_chr2$variant_id,split="_"))))
eqtl_Artery_chr2$chr<-tmp$V1
eqtl_Artery_chr2$loci<-tmp$V2

eqtl_Artery_chr2$index<-paste0(eqtl_Artery_chr2$chr,"-",eqtl_Artery_chr2$loci)
Artery_chr2$index<-paste0(Artery_chr2$chr,"-",Artery_chr2$pig_pos)
same_Artery_chr2<-intersect(Artery_chr2$pig_pos,eqtl_Artery_chr2$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery_chr3$variant_id,split="_"))))
eqtl_Artery_chr3$chr<-tmp$V1
eqtl_Artery_chr3$loci<-tmp$V2

eqtl_Artery_chr3$index<-paste0(eqtl_Artery_chr3$chr,"-",eqtl_Artery_chr3$loci)
Artery_chr3$index<-paste0(Artery_chr3$chr,"-",Artery_chr3$pig_pos)
same_Artery_chr3<-intersect(Artery_chr3$pig_pos,eqtl_Artery_chr3$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery_chr4$variant_id,split="_"))))
eqtl_Artery_chr4$chr<-tmp$V1
eqtl_Artery_chr4$loci<-tmp$V2

eqtl_Artery_chr4$index<-paste0(eqtl_Artery_chr4$chr,"-",eqtl_Artery_chr4$loci)
Artery_chr4$index<-paste0(Artery_chr4$chr,"-",Artery_chr4$pig_pos)
same_Artery_chr4<-intersect(Artery_chr4$pig_pos,eqtl_Artery_chr4$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery_chr5$variant_id,split="_"))))
eqtl_Artery_chr5$loci<-tmp$V2
same_Artery_chr5<-intersect(Artery_chr5$pig_pos,eqtl_Artery_chr5$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery_chr6$variant_id,split="_"))))
eqtl_Artery_chr6$loci<-tmp$V2
same_Artery_chr6<-intersect(Artery_chr6$pig_pos,eqtl_Artery_chr6$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery_chr7$variant_id,split="_"))))
eqtl_Artery_chr7$loci<-tmp$V2
same_Artery_chr7<-intersect(Artery_chr7$pig_pos,eqtl_Artery_chr7$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery_chr8$variant_id,split="_"))))
eqtl_Artery_chr8$loci<-tmp$V2
same_Artery_chr8<-intersect(Artery_chr8$pig_pos,eqtl_Artery_chr8$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery_chr9$variant_id,split="_"))))
eqtl_Artery_chr9$loci<-tmp$V2
same_Artery_chr9<-intersect(Artery_chr9$pig_pos,eqtl_Artery_chr9$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery_chr10$variant_id,split="_"))))
eqtl_Artery_chr10$loci<-tmp$V2
same_Artery_chr10<-intersect(Artery_chr10$pig_pos,eqtl_Artery_chr10$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery_chr11$variant_id,split="_"))))
eqtl_Artery_chr11$loci<-tmp$V2
same_Artery_chr11<-intersect(Artery_chr11$pig_pos,eqtl_Artery_chr11$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery_chr12$variant_id,split="_"))))
eqtl_Artery_chr12$loci<-tmp$V2
same_Artery_chr12<-intersect(Artery_chr12$pig_pos,eqtl_Artery_chr12$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery_chr13$variant_id,split="_"))))
eqtl_Artery_chr13$loci<-tmp$V2
same_Artery_chr13<-intersect(Artery_chr13$pig_pos,eqtl_Artery_chr13$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery_chr14$variant_id,split="_"))))
eqtl_Artery_chr14$loci<-tmp$V2
same_Artery_chr14<-intersect(Artery_chr14$pig_pos,eqtl_Artery_chr14$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery_chr15$variant_id,split="_"))))
eqtl_Artery_chr15$loci<-tmp$V2
same_Artery_chr15<-intersect(Artery_chr15$pig_pos,eqtl_Artery_chr15$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery_chr16$variant_id,split="_"))))
eqtl_Artery_chr16$loci<-tmp$V2
same_Artery_chr16<-intersect(Artery_chr16$pig_pos,eqtl_Artery_chr16$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery_chr17$variant_id,split="_"))))
eqtl_Artery_chr17$loci<-tmp$V2
same_Artery_chr17<-intersect(Artery_chr17$pig_pos,eqtl_Artery_chr17$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Artery_chr18$variant_id,split="_"))))
eqtl_Artery_chr18$loci<-tmp$V2
same_Artery_chr18<-intersect(Artery_chr18$pig_pos,eqtl_Artery_chr18$loci)

chr1_Artery_genes<-NULL
chr1_Artery_overloci<-NULL
if(length(same_Artery_chr1!=0)){
  for(i in 1:length(same_Artery_chr1)){
    a<-eqtl_Artery_chr1$phenotype_id[grep(same_Artery_chr1[i],eqtl_Artery_chr1$loci)]
    b<-eqtl_Artery_chr1[grep(same_Artery_chr1[i],eqtl_Artery_chr1$loci),]
    chr1_Artery_genes<-c(chr1_Artery_genes,a)
    chr1_Artery_overloci<-rbind(chr1_Artery_overloci,b)
  }
}

if(length(same_Artery_chr2!=0)){
  chr2_Artery_genes<-NULL
  chr2_Artery_overloci<-NULL
  for(i in 1:length(same_Artery_chr2)){
    a<-eqtl_Artery_chr2$phenotype_id[grep(same_Artery_chr2[i],eqtl_Artery_chr2$loci)]
    b<-eqtl_Artery_chr2[grep(same_Artery_chr2[i],eqtl_Artery_chr2$loci),]
    chr2_Artery_genes<-c(chr2_Artery_genes,a)
    chr2_Artery_overloci<-rbind(chr2_Artery_overloci,b)
  }
}

if(length(same_Artery_chr3!=0)){
  chr3_Artery_genes<-NULL
  chr3_Artery_overloci<-NULL
  for(i in 1:length(same_Artery_chr3)){
    a<-eqtl_Artery_chr3$phenotype_id[grep(same_Artery_chr3[i],eqtl_Artery_chr3$loci)]
    b<-eqtl_Artery_chr3[grep(same_Artery_chr3[i],eqtl_Artery_chr3$loci),]
    chr3_Artery_genes<-c(chr3_Artery_genes,a)
    chr3_Artery_overloci<-rbind(chr3_Artery_overloci,b)
  }
}

if(length(same_Artery_chr4!=0)){
  chr4_Artery_genes<-NULL
  chr4_Artery_overloci<-NULL
  for(i in 1:length(same_Artery_chr4)){
    a<-eqtl_Artery_chr4$phenotype_id[grep(same_Artery_chr4[i],eqtl_Artery_chr4$loci)]
    b<-eqtl_Artery_chr4[grep(same_Artery_chr4[i],eqtl_Artery_chr4$loci),]
    chr4_Artery_genes<-c(chr4_Artery_genes,a)
    chr4_Artery_overloci<-rbind(chr4_Artery_overloci,b)
  }
}

if(length(same_Artery_chr5!=0)){
  chr5_Artery_genes<-NULL
  chr5_Artery_overloci<-NULL
  for(i in 1:length(same_Artery_chr5)){
    a<-eqtl_Artery_chr5$phenotype_id[grep(same_Artery_chr5[i],eqtl_Artery_chr5$loci)]
    b<-eqtl_Artery_chr5[grep(same_Artery_chr5[i],eqtl_Artery_chr5$loci),]
    chr5_Artery_genes<-c(chr5_Artery_genes,a)
    chr5_Artery_overloci<-rbind(chr5_Artery_overloci,b)
  }
}

if(length(same_Artery_chr6!=0)){
  chr6_Artery_genes<-NULL
  chr6_Artery_overloci<-NULL
  for(i in 1:length(same_Artery_chr6)){
    a<-eqtl_Artery_chr6$phenotype_id[grep(same_Artery_chr6[i],eqtl_Artery_chr6$loci)]
    b<-eqtl_Artery_chr6[grep(same_Artery_chr6[i],eqtl_Artery_chr6$loci),]
    chr6_Artery_genes<-c(chr6_Artery_genes,a)
    chr6_Artery_overloci<-rbind(chr6_Artery_overloci,b)
  }
}

if(length(same_Artery_chr7!=0)){
  chr7_Artery_genes<-NULL
  chr7_Artery_overloci<-NULL
  for(i in 1:length(same_Artery_chr7)){
    a<-eqtl_Artery_chr7$phenotype_id[grep(same_Artery_chr7[i],eqtl_Artery_chr7$loci)]
    b<-eqtl_Artery_chr7[grep(same_Artery_chr7[i],eqtl_Artery_chr7$loci),]
    chr7_Artery_genes<-c(chr7_Artery_genes,a)
    chr7_Artery_overloci<-rbind(chr7_Artery_overloci,b)
  }
}

if(length(same_Artery_chr8!=0)){
  chr8_Artery_genes<-NULL
  chr8_Artery_overloci<-NULL
  for(i in 1:length(same_Artery_chr8)){
    a<-eqtl_Artery_chr8$phenotype_id[grep(same_Artery_chr8[i],eqtl_Artery_chr8$loci)]
    b<-eqtl_Artery_chr8[grep(same_Artery_chr8[i],eqtl_Artery_chr8$loci),]
    chr8_Artery_genes<-c(chr8_Artery_genes,a)
    chr8_Artery_overloci<-rbind(chr8_Artery_overloci,b)
  }
}

if(length(same_Artery_chr9!=0)){
  chr9_Artery_genes<-NULL
  chr9_Artery_overloci<-NULL
  for(i in 1:length(same_Artery_chr9)){
    a<-eqtl_Artery_chr9$phenotype_id[grep(same_Artery_chr9[i],eqtl_Artery_chr9$loci)]
    b<-eqtl_Artery_chr9[grep(same_Artery_chr9[i],eqtl_Artery_chr9$loci),]
    chr9_Artery_genes<-c(chr9_Artery_genes,a)
    chr9_Artery_overloci<-rbind(chr9_Artery_overloci,b)
  }
}

if(length(same_Artery_chr10!=0)){
  chr10_Artery_genes<-NULL
  chr10_Artery_overloci<-NULL
  for(i in 1:length(same_Artery_chr10)){
    a<-eqtl_Artery_chr10$phenotype_id[grep(same_Artery_chr10[i],eqtl_Artery_chr10$loci)]
    b<-eqtl_Artery_chr10[grep(same_Artery_chr10[i],eqtl_Artery_chr10$loci),]
    chr10_Artery_genes<-c(chr10_Artery_genes,a)
    chr10_Artery_overloci<-rbind(chr10_Artery_overloci,b)
  }
}

if(length(same_Artery_chr11!=0)){
  chr11_Artery_genes<-NULL
  chr11_Artery_overloci<-NULL
  for(i in 1:length(same_Artery_chr11)){
    a<-eqtl_Artery_chr11$phenotype_id[grep(same_Artery_chr11[i],eqtl_Artery_chr11$loci)]
    b<-eqtl_Artery_chr11[grep(same_Artery_chr11[i],eqtl_Artery_chr11$loci),]
    chr11_Artery_genes<-c(chr11_Artery_genes,a)
    chr11_Artery_overloci<-rbind(chr11_Artery_overloci,b)
  }
}

chr12_Artery_genes<-NULL
chr12_Artery_overloci<-NULL
if(length(same_Artery_chr12!=0)){
  for(i in 1:length(same_Artery_chr12)){
    a<-eqtl_Artery_chr12$phenotype_id[grep(same_Artery_chr12[i],eqtl_Artery_chr12$loci)]
    b<-eqtl_Artery_chr12[grep(same_Artery_chr12[i],eqtl_Artery_chr12$loci),]
    chr12_Artery_genes<-c(chr12_Artery_genes,a)
    chr12_Artery_overloci<-rbind(chr12_Artery_overloci,b)
  }
}

chr13_Artery_genes<-NULL
chr13_Artery_overloci<-NULL
if(length(same_Artery_chr13!=0)){
  for(i in 1:length(same_Artery_chr13)){
    a<-eqtl_Artery_chr13$phenotype_id[grep(same_Artery_chr13[i],eqtl_Artery_chr13$loci)]
    b<-eqtl_Artery_chr13[grep(same_Artery_chr13[i],eqtl_Artery_chr13$loci),]
    chr13_Artery_genes<-c(chr13_Artery_genes,a)
    chr13_Artery_overloci<-rbind(chr13_Artery_overloci,b)
  }
}

chr14_Artery_genes<-NULL
chr14_Artery_overloci<-NULL
if(length(same_Artery_chr14!=0)){
  for(i in 1:length(same_Artery_chr14)){
    a<-eqtl_Artery_chr14$phenotype_id[grep(same_Artery_chr14[i],eqtl_Artery_chr14$loci)]
    b<-eqtl_Artery_chr14[grep(same_Artery_chr14[i],eqtl_Artery_chr14$loci),]
    chr14_Artery_genes<-c(chr14_Artery_genes,a)
    chr14_Artery_overloci<-rbind(chr14_Artery_overloci,b)
  }
}

chr15_Artery_genes<-NULL
chr15_Artery_overloci<-NULL
if(length(same_Artery_chr15!=0)){
  for(i in 1:length(same_Artery_chr15)){
    a<-eqtl_Artery_chr15$phenotype_id[grep(same_Artery_chr15[i],eqtl_Artery_chr15$loci)]
    b<-eqtl_Artery_chr15[grep(same_Artery_chr15[i],eqtl_Artery_chr15$loci),]
    chr15_Artery_genes<-c(chr15_Artery_genes,a)
    chr15_Artery_overloci<-rbind(chr15_Artery_overloci,b)
  }
}

chr16_Artery_genes<-NULL
chr16_Artery_overloci<-NULL
if(length(same_Artery_chr16!=0)){
  for(i in 1:length(same_Artery_chr16)){
    a<-eqtl_Artery_chr16$phenotype_id[grep(same_Artery_chr16[i],eqtl_Artery_chr16$loci)]
    b<-eqtl_Artery_chr16[grep(same_Artery_chr16[i],eqtl_Artery_chr16$loci),]
    chr16_Artery_genes<-c(chr16_Artery_genes,a)
    chr16_Artery_overloci<-rbind(chr16_Artery_overloci,b)
  }
}

chr17_Artery_genes<-NULL
chr17_Artery_overloci<-NULL
if(length(same_Artery_chr17!=0)){
  for(i in 1:length(same_Artery_chr17)){
    a<-eqtl_Artery_chr17$phenotype_id[grep(same_Artery_chr17[i],eqtl_Artery_chr17$loci)]
    b<-eqtl_Artery_chr17[grep(same_Artery_chr17[i],eqtl_Artery_chr17$loci),]
    chr17_Artery_genes<-c(chr17_Artery_genes,a)
    chr17_Artery_overloci<-rbind(chr17_Artery_overloci,b)
  }
}

chr18_Artery_genes<-NULL
chr18_Artery_overloci<-NULL
if(length(same_Artery_chr18!=0)){
  for(i in 1:length(same_Artery_chr18)){
    a<-eqtl_Artery_chr18$phenotype_id[grep(same_Artery_chr18[i],eqtl_Artery_chr18$loci)]
    b<-eqtl_Artery_chr18[grep(same_Artery_chr18[i],eqtl_Artery_chr18$loci),]
    chr18_Artery_genes<-c(chr18_Artery_genes,a)
    chr18_Artery_overloci<-rbind(chr18_Artery_overloci,b)
  }
}
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
one2one_pig<-annotation$`Pig gene stable ID`[which(annotation$`Pig homology type`=="ortholog_one2one")]

chr1_Artery_overloci_one2one<-chr1_Artery_overloci[match(intersect(chr1_Artery_genes,one2one_pig),chr1_Artery_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Artery_chr1$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Artery_chr1$hum_chr<-tmp$V1
Artery_chr1$hum_loci<-tmp$V3
FM_Artery_hum$index<-paste0(FM_Artery_hum$chr,"_",FM_Artery_hum$variant_pos)
Artery_chr1$index<-paste0("chr",Artery_chr1$hum_chr,"_",Artery_chr1$hum_loci)
chr1_Artery_hum<-NULL
for(i in 1:length(same_Artery_chr1)){
  a<-Artery_chr1[match(same_Artery_chr1[i],Artery_chr1$pig_pos),]
  b<-FM_Artery_hum[match(intersect(a$index,FM_Artery_hum$index),FM_Artery_hum$index),]
  chr1_Artery_hum<-rbind(chr1_Artery_hum,b)
}
chr1_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr1_Artery_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr1_Artery_hum_one2one<-chr1_Artery_hum[match(intersect(chr1_pig2hum_one2one,chr1_Artery_hum$gene_id),chr1_Artery_hum$gene_id),]

chr2_Artery_overloci_one2one<-chr2_Artery_overloci[match(intersect(chr2_Artery_genes,one2one_pig),chr2_Artery_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Artery_chr2$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Artery_chr2$hum_chr<-tmp$V1
Artery_chr2$hum_loci<-tmp$V3
FM_Artery_hum$index<-paste0(FM_Artery_hum$chr,"_",FM_Artery_hum$variant_pos)
Artery_chr2$index<-paste0("chr",Artery_chr2$hum_chr,"_",Artery_chr2$hum_loci)
chr2_Artery_hum<-NULL
for(i in 1:length(same_Artery_chr2)){
  a<-Artery_chr2[match(same_Artery_chr2[i],Artery_chr2$pig_pos),]
  b<-FM_Artery_hum[match(intersect(a$index,FM_Artery_hum$index),FM_Artery_hum$index),]
  chr2_Artery_hum<-rbind(chr2_Artery_hum,b)
}
chr2_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr2_Artery_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr2_Artery_hum_one2one<-chr2_Artery_hum[match(intersect(chr2_pig2hum_one2one,chr2_Artery_hum$gene_id),chr2_Artery_hum$gene_id),]

chr3_Artery_overloci_one2one<-chr3_Artery_overloci[match(intersect(chr3_Artery_genes,one2one_pig),chr3_Artery_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Artery_chr3$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Artery_chr3$hum_chr<-tmp$V1
Artery_chr3$hum_loci<-tmp$V3
FM_Artery_hum$index<-paste0(FM_Artery_hum$chr,"_",FM_Artery_hum$variant_pos)
Artery_chr3$index<-paste0("chr",Artery_chr3$hum_chr,"_",Artery_chr3$hum_loci)
chr3_Artery_hum<-NULL
for(i in 1:length(same_Artery_chr3)){
  a<-Artery_chr3[match(same_Artery_chr3[i],Artery_chr3$pig_pos),]
  b<-FM_Artery_hum[match(intersect(a$index,FM_Artery_hum$index),FM_Artery_hum$index),]
  chr3_Artery_hum<-rbind(chr3_Artery_hum,b)
}
chr3_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr3_Artery_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr3_Artery_hum_one2one<-chr3_Artery_hum[match(intersect(chr3_pig2hum_one2one,chr3_Artery_hum$gene_id),chr3_Artery_hum$gene_id),]

chr4_Artery_overloci_one2one<-chr4_Artery_overloci[match(intersect(chr4_Artery_genes,one2one_pig),chr4_Artery_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Artery_chr4$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Artery_chr4$hum_chr<-tmp$V1
Artery_chr4$hum_loci<-tmp$V3
FM_Artery_hum$index<-paste0(FM_Artery_hum$chr,"_",FM_Artery_hum$variant_pos)
Artery_chr4$index<-paste0("chr",Artery_chr4$hum_chr,"_",Artery_chr4$hum_loci)
chr4_Artery_hum<-NULL
for(i in 1:length(same_Artery_chr4)){
  a<-Artery_chr4[match(same_Artery_chr4[i],Artery_chr4$pig_pos),]
  b<-FM_Artery_hum[match(intersect(a$index,FM_Artery_hum$index),FM_Artery_hum$index),]
  chr4_Artery_hum<-rbind(chr4_Artery_hum,b)
}
chr4_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr4_Artery_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr4_Artery_hum_one2one<-chr4_Artery_hum[match(intersect(chr4_pig2hum_one2one,chr4_Artery_hum$gene_id),chr4_Artery_hum$gene_id),]

chr5_Artery_overloci_one2one<-chr5_Artery_overloci[match(intersect(chr5_Artery_genes,one2one_pig),chr5_Artery_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Artery_chr5$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Artery_chr5$hum_chr<-tmp$V1
Artery_chr5$hum_loci<-tmp$V3
FM_Artery_hum$index<-paste0(FM_Artery_hum$chr,"_",FM_Artery_hum$variant_pos)
Artery_chr5$index<-paste0("chr",Artery_chr5$hum_chr,"_",Artery_chr5$hum_loci)
chr5_Artery_hum<-NULL
for(i in 1:length(same_Artery_chr5)){
  a<-Artery_chr5[match(same_Artery_chr5[i],Artery_chr5$pig_pos),]
  b<-FM_Artery_hum[match(intersect(a$index,FM_Artery_hum$index),FM_Artery_hum$index),]
  chr5_Artery_hum<-rbind(chr5_Artery_hum,b)
}
chr5_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr5_Artery_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr5_Artery_hum_one2one<-chr5_Artery_hum[match(intersect(chr5_pig2hum_one2one,chr5_Artery_hum$gene_id),chr5_Artery_hum$gene_id),]

chr6_Artery_overloci_one2one<-chr6_Artery_overloci[match(intersect(chr6_Artery_genes,one2one_pig),chr6_Artery_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Artery_chr6$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Artery_chr6$hum_chr<-tmp$V1
Artery_chr6$hum_loci<-tmp$V3
FM_Artery_hum$index<-paste0(FM_Artery_hum$chr,"_",FM_Artery_hum$variant_pos)
Artery_chr6$index<-paste0("chr",Artery_chr6$hum_chr,"_",Artery_chr6$hum_loci)
chr6_Artery_hum<-NULL
for(i in 1:length(same_Artery_chr6)){
  a<-Artery_chr6[match(same_Artery_chr6[i],Artery_chr6$pig_pos),]
  b<-FM_Artery_hum[match(intersect(a$index,FM_Artery_hum$index),FM_Artery_hum$index),]
  chr6_Artery_hum<-rbind(chr6_Artery_hum,b)
}
chr6_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr6_Artery_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr6_Artery_hum_one2one<-chr6_Artery_hum[match(intersect(chr6_pig2hum_one2one,chr6_Artery_hum$gene_id),chr6_Artery_hum$gene_id),]

chr7_Artery_overloci_one2one<-chr7_Artery_overloci[match(intersect(chr7_Artery_genes,one2one_pig),chr7_Artery_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Artery_chr7$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Artery_chr7$hum_chr<-tmp$V1
Artery_chr7$hum_loci<-tmp$V3
FM_Artery_hum$index<-paste0(FM_Artery_hum$chr,"_",FM_Artery_hum$variant_pos)
Artery_chr7$index<-paste0("chr",Artery_chr7$hum_chr,"_",Artery_chr7$hum_loci)
chr7_Artery_hum<-NULL
for(i in 1:length(same_Artery_chr7)){
  a<-Artery_chr7[match(same_Artery_chr7[i],Artery_chr7$pig_pos),]
  b<-FM_Artery_hum[match(intersect(a$index,FM_Artery_hum$index),FM_Artery_hum$index),]
  chr7_Artery_hum<-rbind(chr7_Artery_hum,b)
}
chr7_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr7_Artery_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr7_Artery_hum_one2one<-chr7_Artery_hum[match(intersect(chr7_pig2hum_one2one,chr7_Artery_hum$gene_id),chr7_Artery_hum$gene_id),]

chr8_Artery_overloci_one2one<-chr8_Artery_overloci[match(intersect(chr8_Artery_genes,one2one_pig),chr8_Artery_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Artery_chr8$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Artery_chr8$hum_chr<-tmp$V1
Artery_chr8$hum_loci<-tmp$V3
FM_Artery_hum$index<-paste0(FM_Artery_hum$chr,"_",FM_Artery_hum$variant_pos)
Artery_chr8$index<-paste0("chr",Artery_chr8$hum_chr,"_",Artery_chr8$hum_loci)
chr8_Artery_hum<-NULL
for(i in 1:length(same_Artery_chr8)){
  a<-Artery_chr8[match(same_Artery_chr8[i],Artery_chr8$pig_pos),]
  b<-FM_Artery_hum[match(intersect(a$index,FM_Artery_hum$index),FM_Artery_hum$index),]
  chr8_Artery_hum<-rbind(chr8_Artery_hum,b)
}
chr8_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr8_Artery_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr8_Artery_hum_one2one<-chr8_Artery_hum[match(intersect(chr8_pig2hum_one2one,chr8_Artery_hum$gene_id),chr8_Artery_hum$gene_id),]

chr9_Artery_overloci_one2one<-chr9_Artery_overloci[match(intersect(chr9_Artery_genes,one2one_pig),chr9_Artery_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Artery_chr9$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Artery_chr9$hum_chr<-tmp$V1
Artery_chr9$hum_loci<-tmp$V3
FM_Artery_hum$index<-paste0(FM_Artery_hum$chr,"_",FM_Artery_hum$variant_pos)
Artery_chr9$index<-paste0("chr",Artery_chr9$hum_chr,"_",Artery_chr9$hum_loci)
chr9_Artery_hum<-NULL
for(i in 1:length(same_Artery_chr9)){
  a<-Artery_chr9[match(same_Artery_chr9[i],Artery_chr9$pig_pos),]
  b<-FM_Artery_hum[match(intersect(a$index,FM_Artery_hum$index),FM_Artery_hum$index),]
  chr9_Artery_hum<-rbind(chr9_Artery_hum,b)
}
chr9_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr9_Artery_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr9_Artery_hum_one2one<-chr9_Artery_hum[match(intersect(chr9_pig2hum_one2one,chr9_Artery_hum$gene_id),chr9_Artery_hum$gene_id),]

chr10_Artery_overloci_one2one<-chr10_Artery_overloci[match(intersect(chr10_Artery_genes,one2one_pig),chr10_Artery_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Artery_chr10$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Artery_chr10$hum_chr<-tmp$V1
Artery_chr10$hum_loci<-tmp$V3
FM_Artery_hum$index<-paste0(FM_Artery_hum$chr,"_",FM_Artery_hum$variant_pos)
Artery_chr10$index<-paste0("chr",Artery_chr10$hum_chr,"_",Artery_chr10$hum_loci)
chr10_Artery_hum<-NULL
for(i in 1:length(same_Artery_chr10)){
  a<-Artery_chr10[match(same_Artery_chr10[i],Artery_chr10$pig_pos),]
  b<-FM_Artery_hum[match(intersect(a$index,FM_Artery_hum$index),FM_Artery_hum$index),]
  chr10_Artery_hum<-rbind(chr10_Artery_hum,b)
}
chr10_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr10_Artery_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr10_Artery_hum_one2one<-chr10_Artery_hum[match(intersect(chr10_pig2hum_one2one,chr10_Artery_hum$gene_id),chr10_Artery_hum$gene_id),]

chr11_Artery_overloci_one2one<-chr11_Artery_overloci[match(intersect(chr11_Artery_genes,one2one_pig),chr11_Artery_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Artery_chr11$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Artery_chr11$hum_chr<-tmp$V1
Artery_chr11$hum_loci<-tmp$V3
FM_Artery_hum$index<-paste0(FM_Artery_hum$chr,"_",FM_Artery_hum$variant_pos)
Artery_chr11$index<-paste0("chr",Artery_chr11$hum_chr,"_",Artery_chr11$hum_loci)
chr11_Artery_hum<-NULL
for(i in 1:length(same_Artery_chr11)){
  a<-Artery_chr11[match(same_Artery_chr11[i],Artery_chr11$pig_pos),]
  b<-FM_Artery_hum[match(intersect(a$index,FM_Artery_hum$index),FM_Artery_hum$index),]
  chr11_Artery_hum<-rbind(chr11_Artery_hum,b)
}
chr11_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr11_Artery_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr11_Artery_hum_one2one<-chr11_Artery_hum[match(intersect(chr11_pig2hum_one2one,chr11_Artery_hum$gene_id),chr11_Artery_hum$gene_id),]

chr12_Artery_overloci_one2one<-chr12_Artery_overloci[match(intersect(chr12_Artery_genes,one2one_pig),chr12_Artery_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Artery_chr12$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Artery_chr12$hum_chr<-tmp$V1
Artery_chr12$hum_loci<-tmp$V3
FM_Artery_hum$index<-paste0(FM_Artery_hum$chr,"_",FM_Artery_hum$variant_pos)
Artery_chr12$index<-paste0("chr",Artery_chr12$hum_chr,"_",Artery_chr12$hum_loci)
chr12_Artery_hum<-NULL
for(i in 1:length(same_Artery_chr12)){
  a<-Artery_chr12[match(same_Artery_chr12[i],Artery_chr12$pig_pos),]
  b<-FM_Artery_hum[match(intersect(a$index,FM_Artery_hum$index),FM_Artery_hum$index),]
  chr12_Artery_hum<-rbind(chr12_Artery_hum,b)
}
chr12_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr12_Artery_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr12_Artery_hum_one2one<-chr12_Artery_hum[match(intersect(chr12_pig2hum_one2one,chr12_Artery_hum$gene_id),chr12_Artery_hum$gene_id),]

chr13_Artery_overloci_one2one<-chr13_Artery_overloci[match(intersect(chr13_Artery_genes,one2one_pig),chr13_Artery_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Artery_chr13$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Artery_chr13$hum_chr<-tmp$V1
Artery_chr13$hum_loci<-tmp$V3
FM_Artery_hum$index<-paste0(FM_Artery_hum$chr,"_",FM_Artery_hum$variant_pos)
Artery_chr13$index<-paste0("chr",Artery_chr13$hum_chr,"_",Artery_chr13$hum_loci)
chr13_Artery_hum<-NULL
for(i in 1:length(same_Artery_chr13)){
  a<-Artery_chr13[match(same_Artery_chr13[i],Artery_chr13$pig_pos),]
  b<-FM_Artery_hum[match(intersect(a$index,FM_Artery_hum$index),FM_Artery_hum$index),]
  chr13_Artery_hum<-rbind(chr13_Artery_hum,b)
}
chr13_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr13_Artery_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr13_Artery_hum_one2one<-chr13_Artery_hum[match(intersect(chr13_pig2hum_one2one,chr13_Artery_hum$gene_id),chr13_Artery_hum$gene_id),]

chr14_Artery_overloci_one2one<-chr14_Artery_overloci[match(intersect(chr14_Artery_genes,one2one_pig),chr14_Artery_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Artery_chr14$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Artery_chr14$hum_chr<-tmp$V1
Artery_chr14$hum_loci<-tmp$V3
FM_Artery_hum$index<-paste0(FM_Artery_hum$chr,"_",FM_Artery_hum$variant_pos)
Artery_chr14$index<-paste0("chr",Artery_chr14$hum_chr,"_",Artery_chr14$hum_loci)
chr14_Artery_hum<-NULL
for(i in 1:length(same_Artery_chr14)){
  a<-Artery_chr14[match(same_Artery_chr14[i],Artery_chr14$pig_pos),]
  b<-FM_Artery_hum[match(intersect(a$index,FM_Artery_hum$index),FM_Artery_hum$index),]
  chr14_Artery_hum<-rbind(chr14_Artery_hum,b)
}
chr14_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr14_Artery_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr14_Artery_hum_one2one<-chr14_Artery_hum[match(intersect(chr14_pig2hum_one2one,chr14_Artery_hum$gene_id),chr14_Artery_hum$gene_id),]

chr15_Artery_overloci_one2one<-chr15_Artery_overloci[match(intersect(chr15_Artery_genes,one2one_pig),chr15_Artery_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Artery_chr15$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Artery_chr15$hum_chr<-tmp$V1
Artery_chr15$hum_loci<-tmp$V3
FM_Artery_hum$index<-paste0(FM_Artery_hum$chr,"_",FM_Artery_hum$variant_pos)
Artery_chr15$index<-paste0("chr",Artery_chr15$hum_chr,"_",Artery_chr15$hum_loci)
chr15_Artery_hum<-NULL
for(i in 1:length(same_Artery_chr15)){
  a<-Artery_chr15[match(same_Artery_chr15[i],Artery_chr15$pig_pos),]
  b<-FM_Artery_hum[match(intersect(a$index,FM_Artery_hum$index),FM_Artery_hum$index),]
  chr15_Artery_hum<-rbind(chr15_Artery_hum,b)
}
chr15_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr15_Artery_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr15_Artery_hum_one2one<-chr15_Artery_hum[match(intersect(chr15_pig2hum_one2one,chr15_Artery_hum$gene_id),chr15_Artery_hum$gene_id),]

chr16_Artery_overloci_one2one<-chr16_Artery_overloci[match(intersect(chr16_Artery_genes,one2one_pig),chr16_Artery_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Artery_chr16$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Artery_chr16$hum_chr<-tmp$V1
Artery_chr16$hum_loci<-tmp$V3
FM_Artery_hum$index<-paste0(FM_Artery_hum$chr,"_",FM_Artery_hum$variant_pos)
Artery_chr16$index<-paste0("chr",Artery_chr16$hum_chr,"_",Artery_chr16$hum_loci)
chr16_Artery_hum<-NULL
for(i in 1:length(same_Artery_chr16)){
  a<-Artery_chr16[match(same_Artery_chr16[i],Artery_chr16$pig_pos),]
  b<-FM_Artery_hum[match(intersect(a$index,FM_Artery_hum$index),FM_Artery_hum$index),]
  chr16_Artery_hum<-rbind(chr16_Artery_hum,b)
}
chr16_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr16_Artery_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr16_Artery_hum_one2one<-chr16_Artery_hum[match(intersect(chr16_pig2hum_one2one,chr16_Artery_hum$gene_id),chr16_Artery_hum$gene_id),]

chr17_Artery_overloci_one2one<-chr17_Artery_overloci[match(intersect(chr17_Artery_genes,one2one_pig),chr17_Artery_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Artery_chr17$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Artery_chr17$hum_chr<-tmp$V1
Artery_chr17$hum_loci<-tmp$V3
FM_Artery_hum$index<-paste0(FM_Artery_hum$chr,"_",FM_Artery_hum$variant_pos)
Artery_chr17$index<-paste0("chr",Artery_chr17$hum_chr,"_",Artery_chr17$hum_loci)
chr17_Artery_hum<-NULL
for(i in 1:length(same_Artery_chr17)){
  a<-Artery_chr17[match(same_Artery_chr17[i],Artery_chr17$pig_pos),]
  b<-FM_Artery_hum[match(intersect(a$index,FM_Artery_hum$index),FM_Artery_hum$index),]
  chr17_Artery_hum<-rbind(chr17_Artery_hum,b)
}
chr17_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr17_Artery_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr17_Artery_hum_one2one<-chr17_Artery_hum[match(intersect(chr17_pig2hum_one2one,chr17_Artery_hum$gene_id),chr17_Artery_hum$gene_id),]

chr18_Artery_overloci_one2one<-chr18_Artery_overloci[match(intersect(chr18_Artery_genes,one2one_pig),chr18_Artery_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Artery_chr18$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Artery_chr18$hum_chr<-tmp$V1
Artery_chr18$hum_loci<-tmp$V3
FM_Artery_hum$index<-paste0(FM_Artery_hum$chr,"_",FM_Artery_hum$variant_pos)
Artery_chr18$index<-paste0("chr",Artery_chr18$hum_chr,"_",Artery_chr18$hum_loci)
chr18_Artery_hum<-NULL
for(i in 1:length(same_Artery_chr18)){
  a<-Artery_chr18[match(same_Artery_chr18[i],Artery_chr18$pig_pos),]
  b<-FM_Artery_hum[match(intersect(a$index,FM_Artery_hum$index),FM_Artery_hum$index),]
  chr18_Artery_hum<-rbind(chr18_Artery_hum,b)
}
chr18_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr18_Artery_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr18_Artery_hum_one2one<-chr18_Artery_hum[match(intersect(chr18_pig2hum_one2one,chr18_Artery_hum$gene_id),chr18_Artery_hum$gene_id),]

Artery_one2one_SNP_hum<-rbind(chr1_Artery_hum_one2one,chr2_Artery_hum_one2one,chr3_Artery_hum_one2one,chr4_Artery_hum_one2one,chr5_Artery_hum_one2one,
                              chr6_Artery_hum_one2one,chr7_Artery_hum_one2one,chr8_Artery_hum_one2one,chr9_Artery_hum_one2one,chr10_Artery_hum_one2one,
                              chr11_Artery_hum_one2one,chr12_Artery_hum_one2one,chr13_Artery_hum_one2one,chr14_Artery_hum_one2one,chr15_Artery_hum_one2one,
                              chr16_Artery_hum_one2one,chr17_Artery_hum_one2one,chr18_Artery_hum_one2one)

chr1_Artery_pig_one2one<-chr1_Artery_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr1_Artery_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr1_Artery_overloci_one2one$phenotype_id),1:9]
chr2_Artery_pig_one2one<-chr2_Artery_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr2_Artery_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr2_Artery_overloci_one2one$phenotype_id),1:9]
chr3_Artery_pig_one2one<-chr3_Artery_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr3_Artery_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr3_Artery_overloci_one2one$phenotype_id),1:9]
chr4_Artery_pig_one2one<-chr4_Artery_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr4_Artery_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr4_Artery_overloci_one2one$phenotype_id),1:9]
chr5_Artery_pig_one2one<-chr5_Artery_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr5_Artery_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr5_Artery_overloci_one2one$phenotype_id),1:9]
chr6_Artery_pig_one2one<-chr6_Artery_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr6_Artery_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr6_Artery_overloci_one2one$phenotype_id),1:9]
chr7_Artery_pig_one2one<-chr7_Artery_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr7_Artery_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr7_Artery_overloci_one2one$phenotype_id),1:9]
chr8_Artery_pig_one2one<-chr8_Artery_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr8_Artery_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr8_Artery_overloci_one2one$phenotype_id),1:9]
chr9_Artery_pig_one2one<-chr9_Artery_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr9_Artery_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr9_Artery_overloci_one2one$phenotype_id),1:9]
chr10_Artery_pig_one2one<-chr10_Artery_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr10_Artery_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr10_Artery_overloci_one2one$phenotype_id),1:9]
chr11_Artery_pig_one2one<-chr11_Artery_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr11_Artery_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr11_Artery_overloci_one2one$phenotype_id),1:9]
chr12_Artery_pig_one2one<-chr12_Artery_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr12_Artery_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr12_Artery_overloci_one2one$phenotype_id),1:9]
chr13_Artery_pig_one2one<-chr13_Artery_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr13_Artery_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr13_Artery_overloci_one2one$phenotype_id),1:9]
chr14_Artery_pig_one2one<-chr14_Artery_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr14_Artery_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr14_Artery_overloci_one2one$phenotype_id),1:9]
chr15_Artery_pig_one2one<-chr15_Artery_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr15_Artery_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr15_Artery_overloci_one2one$phenotype_id),1:9]
chr16_Artery_pig_one2one<-chr16_Artery_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr16_Artery_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr16_Artery_overloci_one2one$phenotype_id),1:9]
chr17_Artery_pig_one2one<-chr17_Artery_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr17_Artery_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr17_Artery_overloci_one2one$phenotype_id),1:9]
chr18_Artery_pig_one2one<-chr18_Artery_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr18_Artery_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr18_Artery_overloci_one2one$phenotype_id),1:9]

Artery_one2one_SNP_pig<-rbind(chr1_Artery_pig_one2one,chr2_Artery_pig_one2one,chr3_Artery_pig_one2one,chr4_Artery_pig_one2one,chr5_Artery_pig_one2one,
                              chr6_Artery_pig_one2one,chr7_Artery_pig_one2one,chr8_Artery_pig_one2one,chr9_Artery_pig_one2one,chr10_Artery_pig_one2one,
                              chr11_Artery_pig_one2one,chr12_Artery_pig_one2one,chr13_Artery_pig_one2one,chr14_Artery_pig_one2one,chr15_Artery_pig_one2one,
                              chr16_Artery_pig_one2one,chr17_Artery_pig_one2one,chr18_Artery_pig_one2one)

Artery_SNP_sum<-array(NA,dim=c(nrow(Artery_one2one_SNP_hum),2))
colnames(Artery_SNP_sum)<-c("Human","Pig")
Artery_SNP_sum<-as.data.frame(Artery_SNP_sum)
Artery_SNP_sum$Human<-Artery_one2one_SNP_hum$slope / Artery_one2one_SNP_hum$slope_se
Artery_SNP_sum$Pig<-Artery_one2one_SNP_pig$slope / Artery_one2one_SNP_pig$slope_se
cor<-cor(abs(Artery_SNP_sum$Human),abs(Artery_SNP_sum$Pig))
p_val<-t.test(abs(Artery_SNP_sum$Human),abs(Artery_SNP_sum$Pig))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/Artery_SNP_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(Artery_SNP_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The number of eGenes expressed in common tissues",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

save(Artery_one2one_SNP_hum,Artery_one2one_SNP_pig,file="/Users/baizhonghao/Downloads/human-pig GTEx/Artery_SNP.Rdata")

#Blood_SNP_overlaploci#
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

eqtl_Blood_chr1<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Blood/Blood.cis_qtl_pairs.1.txt"))
eqtl_Blood_chr2<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Blood/Blood.cis_qtl_pairs.2.txt"))
eqtl_Blood_chr3<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Blood/Blood.cis_qtl_pairs.3.txt"))
eqtl_Blood_chr4<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Blood/Blood.cis_qtl_pairs.4.txt"))
eqtl_Blood_chr5<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Blood/Blood.cis_qtl_pairs.5.txt"))
eqtl_Blood_chr6<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Blood/Blood.cis_qtl_pairs.6.txt"))
eqtl_Blood_chr7<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Blood/Blood.cis_qtl_pairs.7.txt"))
eqtl_Blood_chr8<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Blood/Blood.cis_qtl_pairs.8.txt"))
eqtl_Blood_chr9<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Blood/Blood.cis_qtl_pairs.9.txt"))
eqtl_Blood_chr10<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Blood/Blood.cis_qtl_pairs.10.txt"))
eqtl_Blood_chr11<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Blood/Blood.cis_qtl_pairs.11.txt"))
eqtl_Blood_chr12<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Blood/Blood.cis_qtl_pairs.12.txt"))
eqtl_Blood_chr13<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Blood/Blood.cis_qtl_pairs.13.txt"))
eqtl_Blood_chr14<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Blood/Blood.cis_qtl_pairs.14.txt"))
eqtl_Blood_chr15<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Blood/Blood.cis_qtl_pairs.15.txt"))
eqtl_Blood_chr16<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Blood/Blood.cis_qtl_pairs.16.txt"))
eqtl_Blood_chr17<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Blood/Blood.cis_qtl_pairs.17.txt"))
eqtl_Blood_chr18<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Blood/Blood.cis_qtl_pairs.18.txt"))

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood_chr1$variant_id,split="_"))))
eqtl_Blood_chr1$chr<-tmp$V1
eqtl_Blood_chr1$loci<-tmp$V2

eqtl_Blood_chr1$index<-paste0(eqtl_Blood_chr1$chr,"-",eqtl_Blood_chr1$loci)
Blood_chr1$index<-paste0(Blood_chr1$chr,"-",Blood_chr1$pig_pos)
same_Blood_chr1<-intersect(Blood_chr1$pig_pos,eqtl_Blood_chr1$loci)


tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood_chr2$variant_id,split="_"))))
eqtl_Blood_chr2$chr<-tmp$V1
eqtl_Blood_chr2$loci<-tmp$V2

eqtl_Blood_chr2$index<-paste0(eqtl_Blood_chr2$chr,"-",eqtl_Blood_chr2$loci)
Blood_chr2$index<-paste0(Blood_chr2$chr,"-",Blood_chr2$pig_pos)
same_Blood_chr2<-intersect(Blood_chr2$pig_pos,eqtl_Blood_chr2$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood_chr3$variant_id,split="_"))))
eqtl_Blood_chr3$chr<-tmp$V1
eqtl_Blood_chr3$loci<-tmp$V2

eqtl_Blood_chr3$index<-paste0(eqtl_Blood_chr3$chr,"-",eqtl_Blood_chr3$loci)
Blood_chr3$index<-paste0(Blood_chr3$chr,"-",Blood_chr3$pig_pos)
same_Blood_chr3<-intersect(Blood_chr3$pig_pos,eqtl_Blood_chr3$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood_chr4$variant_id,split="_"))))
eqtl_Blood_chr4$chr<-tmp$V1
eqtl_Blood_chr4$loci<-tmp$V2

eqtl_Blood_chr4$index<-paste0(eqtl_Blood_chr4$chr,"-",eqtl_Blood_chr4$loci)
Blood_chr4$index<-paste0(Blood_chr4$chr,"-",Blood_chr4$pig_pos)
same_Blood_chr4<-intersect(Blood_chr4$pig_pos,eqtl_Blood_chr4$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood_chr5$variant_id,split="_"))))
eqtl_Blood_chr5$loci<-tmp$V2
same_Blood_chr5<-intersect(Blood_chr5$pig_pos,eqtl_Blood_chr5$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood_chr6$variant_id,split="_"))))
eqtl_Blood_chr6$loci<-tmp$V2
same_Blood_chr6<-intersect(Blood_chr6$pig_pos,eqtl_Blood_chr6$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood_chr7$variant_id,split="_"))))
eqtl_Blood_chr7$loci<-tmp$V2
same_Blood_chr7<-intersect(Blood_chr7$pig_pos,eqtl_Blood_chr7$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood_chr8$variant_id,split="_"))))
eqtl_Blood_chr8$loci<-tmp$V2
same_Blood_chr8<-intersect(Blood_chr8$pig_pos,eqtl_Blood_chr8$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood_chr9$variant_id,split="_"))))
eqtl_Blood_chr9$loci<-tmp$V2
same_Blood_chr9<-intersect(Blood_chr9$pig_pos,eqtl_Blood_chr9$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood_chr10$variant_id,split="_"))))
eqtl_Blood_chr10$loci<-tmp$V2
same_Blood_chr10<-intersect(Blood_chr10$pig_pos,eqtl_Blood_chr10$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood_chr11$variant_id,split="_"))))
eqtl_Blood_chr11$loci<-tmp$V2
same_Blood_chr11<-intersect(Blood_chr11$pig_pos,eqtl_Blood_chr11$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood_chr12$variant_id,split="_"))))
eqtl_Blood_chr12$loci<-tmp$V2
same_Blood_chr12<-intersect(Blood_chr12$pig_pos,eqtl_Blood_chr12$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood_chr13$variant_id,split="_"))))
eqtl_Blood_chr13$loci<-tmp$V2
same_Blood_chr13<-intersect(Blood_chr13$pig_pos,eqtl_Blood_chr13$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood_chr14$variant_id,split="_"))))
eqtl_Blood_chr14$loci<-tmp$V2
same_Blood_chr14<-intersect(Blood_chr14$pig_pos,eqtl_Blood_chr14$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood_chr15$variant_id,split="_"))))
eqtl_Blood_chr15$loci<-tmp$V2
same_Blood_chr15<-intersect(Blood_chr15$pig_pos,eqtl_Blood_chr15$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood_chr16$variant_id,split="_"))))
eqtl_Blood_chr16$loci<-tmp$V2
same_Blood_chr16<-intersect(Blood_chr16$pig_pos,eqtl_Blood_chr16$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood_chr17$variant_id,split="_"))))
eqtl_Blood_chr17$loci<-tmp$V2
same_Blood_chr17<-intersect(Blood_chr17$pig_pos,eqtl_Blood_chr17$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Blood_chr18$variant_id,split="_"))))
eqtl_Blood_chr18$loci<-tmp$V2
same_Blood_chr18<-intersect(Blood_chr18$pig_pos,eqtl_Blood_chr18$loci)

chr1_Blood_genes<-NULL
chr1_Blood_overloci<-NULL
if(length(same_Blood_chr1!=0)){
  for(i in 1:length(same_Blood_chr1)){
    a<-eqtl_Blood_chr1$phenotype_id[grep(same_Blood_chr1[i],eqtl_Blood_chr1$loci)]
    b<-eqtl_Blood_chr1[grep(same_Blood_chr1[i],eqtl_Blood_chr1$loci),]
    chr1_Blood_genes<-c(chr1_Blood_genes,a)
    chr1_Blood_overloci<-rbind(chr1_Blood_overloci,b)
  }
}

chr2_Blood_genes<-NULL
chr2_Blood_overloci<-NULL
if(length(same_Blood_chr2!=0)){
  for(i in 1:length(same_Blood_chr2)){
    a<-eqtl_Blood_chr2$phenotype_id[grep(same_Blood_chr2[i],eqtl_Blood_chr2$loci)]
    b<-eqtl_Blood_chr2[grep(same_Blood_chr2[i],eqtl_Blood_chr2$loci),]
    chr2_Blood_genes<-c(chr2_Blood_genes,a)
    chr2_Blood_overloci<-rbind(chr2_Blood_overloci,b)
  }
}

chr3_Blood_genes<-NULL
chr3_Blood_overloci<-NULL
if(length(same_Blood_chr3!=0)){
  for(i in 1:length(same_Blood_chr3)){
    a<-eqtl_Blood_chr3$phenotype_id[grep(same_Blood_chr3[i],eqtl_Blood_chr3$loci)]
    b<-eqtl_Blood_chr3[grep(same_Blood_chr3[i],eqtl_Blood_chr3$loci),]
    chr3_Blood_genes<-c(chr3_Blood_genes,a)
    chr3_Blood_overloci<-rbind(chr3_Blood_overloci,b)
  }
}

chr4_Blood_genes<-NULL
chr4_Blood_overloci<-NULL
if(length(same_Blood_chr4!=0)){
  for(i in 1:length(same_Blood_chr4)){
    a<-eqtl_Blood_chr4$phenotype_id[grep(same_Blood_chr4[i],eqtl_Blood_chr4$loci)]
    b<-eqtl_Blood_chr4[grep(same_Blood_chr4[i],eqtl_Blood_chr4$loci),]
    chr4_Blood_genes<-c(chr4_Blood_genes,a)
    chr4_Blood_overloci<-rbind(chr4_Blood_overloci,b)
  }
}

chr5_Blood_genes<-NULL
chr5_Blood_overloci<-NULL
if(length(same_Blood_chr5!=0)){
  for(i in 1:length(same_Blood_chr5)){
    a<-eqtl_Blood_chr5$phenotype_id[grep(same_Blood_chr5[i],eqtl_Blood_chr5$loci)]
    b<-eqtl_Blood_chr5[grep(same_Blood_chr5[i],eqtl_Blood_chr5$loci),]
    chr5_Blood_genes<-c(chr5_Blood_genes,a)
    chr5_Blood_overloci<-rbind(chr5_Blood_overloci,b)
  }
}

chr6_Blood_genes<-NULL
chr6_Blood_overloci<-NULL
if(length(same_Blood_chr6!=0)){
  for(i in 1:length(same_Blood_chr6)){
    a<-eqtl_Blood_chr6$phenotype_id[grep(same_Blood_chr6[i],eqtl_Blood_chr6$loci)]
    b<-eqtl_Blood_chr6[grep(same_Blood_chr6[i],eqtl_Blood_chr6$loci),]
    chr6_Blood_genes<-c(chr6_Blood_genes,a)
    chr6_Blood_overloci<-rbind(chr6_Blood_overloci,b)
  }
}

chr7_Blood_genes<-NULL
chr7_Blood_overloci<-NULL
if(length(same_Blood_chr7!=0)){
  for(i in 1:length(same_Blood_chr7)){
    a<-eqtl_Blood_chr7$phenotype_id[grep(same_Blood_chr7[i],eqtl_Blood_chr7$loci)]
    b<-eqtl_Blood_chr7[grep(same_Blood_chr7[i],eqtl_Blood_chr7$loci),]
    chr7_Blood_genes<-c(chr7_Blood_genes,a)
    chr7_Blood_overloci<-rbind(chr7_Blood_overloci,b)
  }
}

chr8_Blood_genes<-NULL
chr8_Blood_overloci<-NULL
if(length(same_Blood_chr8!=0)){
  for(i in 1:length(same_Blood_chr8)){
    a<-eqtl_Blood_chr8$phenotype_id[grep(same_Blood_chr8[i],eqtl_Blood_chr8$loci)]
    b<-eqtl_Blood_chr8[grep(same_Blood_chr8[i],eqtl_Blood_chr8$loci),]
    chr8_Blood_genes<-c(chr8_Blood_genes,a)
    chr8_Blood_overloci<-rbind(chr8_Blood_overloci,b)
  }
}

if(length(same_Blood_chr9!=0)){
  chr9_Blood_genes<-NULL
  chr9_Blood_overloci<-NULL
  for(i in 1:length(same_Blood_chr9)){
    a<-eqtl_Blood_chr9$phenotype_id[grep(same_Blood_chr9[i],eqtl_Blood_chr9$loci)]
    b<-eqtl_Blood_chr9[grep(same_Blood_chr9[i],eqtl_Blood_chr9$loci),]
    chr9_Blood_genes<-c(chr9_Blood_genes,a)
    chr9_Blood_overloci<-rbind(chr9_Blood_overloci,b)
  }
}

chr10_Blood_genes<-NULL
chr10_Blood_overloci<-NULL
if(length(same_Blood_chr10!=0)){
  for(i in 1:length(same_Blood_chr10)){
    a<-eqtl_Blood_chr10$phenotype_id[grep(same_Blood_chr10[i],eqtl_Blood_chr10$loci)]
    b<-eqtl_Blood_chr10[grep(same_Blood_chr10[i],eqtl_Blood_chr10$loci),]
    chr10_Blood_genes<-c(chr10_Blood_genes,a)
    chr10_Blood_overloci<-rbind(chr10_Blood_overloci,b)
  }
}

chr11_Blood_genes<-NULL
chr11_Blood_overloci<-NULL
if(length(same_Blood_chr11!=0)){
  for(i in 1:length(same_Blood_chr11)){
    a<-eqtl_Blood_chr11$phenotype_id[grep(same_Blood_chr11[i],eqtl_Blood_chr11$loci)]
    b<-eqtl_Blood_chr11[grep(same_Blood_chr11[i],eqtl_Blood_chr11$loci),]
    chr11_Blood_genes<-c(chr11_Blood_genes,a)
    chr11_Blood_overloci<-rbind(chr11_Blood_overloci,b)
  }
}

chr12_Blood_genes<-NULL
chr12_Blood_overloci<-NULL
if(length(same_Blood_chr12!=0)){
  for(i in 1:length(same_Blood_chr12)){
    a<-eqtl_Blood_chr12$phenotype_id[grep(same_Blood_chr12[i],eqtl_Blood_chr12$loci)]
    b<-eqtl_Blood_chr12[grep(same_Blood_chr12[i],eqtl_Blood_chr12$loci),]
    chr12_Blood_genes<-c(chr12_Blood_genes,a)
    chr12_Blood_overloci<-rbind(chr12_Blood_overloci,b)
  }
}

chr13_Blood_genes<-NULL
chr13_Blood_overloci<-NULL
if(length(same_Blood_chr13!=0)){
  for(i in 1:length(same_Blood_chr13)){
    a<-eqtl_Blood_chr13$phenotype_id[grep(same_Blood_chr13[i],eqtl_Blood_chr13$loci)]
    b<-eqtl_Blood_chr13[grep(same_Blood_chr13[i],eqtl_Blood_chr13$loci),]
    chr13_Blood_genes<-c(chr13_Blood_genes,a)
    chr13_Blood_overloci<-rbind(chr13_Blood_overloci,b)
  }
}

chr14_Blood_genes<-NULL
chr14_Blood_overloci<-NULL
if(length(same_Blood_chr14!=0)){
  for(i in 1:length(same_Blood_chr14)){
    a<-eqtl_Blood_chr14$phenotype_id[grep(same_Blood_chr14[i],eqtl_Blood_chr14$loci)]
    b<-eqtl_Blood_chr14[grep(same_Blood_chr14[i],eqtl_Blood_chr14$loci),]
    chr14_Blood_genes<-c(chr14_Blood_genes,a)
    chr14_Blood_overloci<-rbind(chr14_Blood_overloci,b)
  }
}

chr15_Blood_genes<-NULL
chr15_Blood_overloci<-NULL
if(length(same_Blood_chr15!=0)){
  for(i in 1:length(same_Blood_chr15)){
    a<-eqtl_Blood_chr15$phenotype_id[grep(same_Blood_chr15[i],eqtl_Blood_chr15$loci)]
    b<-eqtl_Blood_chr15[grep(same_Blood_chr15[i],eqtl_Blood_chr15$loci),]
    chr15_Blood_genes<-c(chr15_Blood_genes,a)
    chr15_Blood_overloci<-rbind(chr15_Blood_overloci,b)
  }
}

chr16_Blood_genes<-NULL
chr16_Blood_overloci<-NULL
if(length(same_Blood_chr16!=0)){
  for(i in 1:length(same_Blood_chr16)){
    a<-eqtl_Blood_chr16$phenotype_id[grep(same_Blood_chr16[i],eqtl_Blood_chr16$loci)]
    b<-eqtl_Blood_chr16[grep(same_Blood_chr16[i],eqtl_Blood_chr16$loci),]
    chr16_Blood_genes<-c(chr16_Blood_genes,a)
    chr16_Blood_overloci<-rbind(chr16_Blood_overloci,b)
  }
}

chr17_Blood_genes<-NULL
chr17_Blood_overloci<-NULL
if(length(same_Blood_chr17!=0)){
  for(i in 1:length(same_Blood_chr17)){
    a<-eqtl_Blood_chr17$phenotype_id[grep(same_Blood_chr17[i],eqtl_Blood_chr17$loci)]
    b<-eqtl_Blood_chr17[grep(same_Blood_chr17[i],eqtl_Blood_chr17$loci),]
    chr17_Blood_genes<-c(chr17_Blood_genes,a)
    chr17_Blood_overloci<-rbind(chr17_Blood_overloci,b)
  }
}

chr18_Blood_genes<-NULL
chr18_Blood_overloci<-NULL
if(length(same_Blood_chr18!=0)){
  for(i in 1:length(same_Blood_chr18)){
    a<-eqtl_Blood_chr18$phenotype_id[grep(same_Blood_chr18[i],eqtl_Blood_chr18$loci)]
    b<-eqtl_Blood_chr18[grep(same_Blood_chr18[i],eqtl_Blood_chr18$loci),]
    chr18_Blood_genes<-c(chr18_Blood_genes,a)
    chr18_Blood_overloci<-rbind(chr18_Blood_overloci,b)
  }
}
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
one2one_pig<-annotation$`Pig gene stable ID`[which(annotation$`Pig homology type`=="ortholog_one2one")]

chr1_Blood_overloci_one2one<-chr1_Blood_overloci[match(intersect(chr1_Blood_genes,one2one_pig),chr1_Blood_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Blood_chr1$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Blood_chr1$hum_chr<-tmp$V1
Blood_chr1$hum_loci<-tmp$V3
FM_Blood_hum$index<-paste0(FM_Blood_hum$chr,"_",FM_Blood_hum$variant_pos)
Blood_chr1$index<-paste0("chr",Blood_chr1$hum_chr,"_",Blood_chr1$hum_loci)
chr1_Blood_hum<-NULL
for(i in 1:length(same_Blood_chr1)){
  a<-Blood_chr1[match(same_Blood_chr1[i],Blood_chr1$pig_pos),]
  b<-FM_Blood_hum[match(intersect(a$index,FM_Blood_hum$index),FM_Blood_hum$index),]
  chr1_Blood_hum<-rbind(chr1_Blood_hum,b)
}
chr1_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr1_Blood_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr1_Blood_hum_one2one<-chr1_Blood_hum[match(intersect(chr1_pig2hum_one2one,chr1_Blood_hum$gene_id),chr1_Blood_hum$gene_id),]

chr2_Blood_overloci_one2one<-chr2_Blood_overloci[match(intersect(chr2_Blood_genes,one2one_pig),chr2_Blood_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Blood_chr2$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Blood_chr2$hum_chr<-tmp$V1
Blood_chr2$hum_loci<-tmp$V3
FM_Blood_hum$index<-paste0(FM_Blood_hum$chr,"_",FM_Blood_hum$variant_pos)
Blood_chr2$index<-paste0("chr",Blood_chr2$hum_chr,"_",Blood_chr2$hum_loci)
chr2_Blood_hum<-NULL
for(i in 1:length(same_Blood_chr2)){
  a<-Blood_chr2[match(same_Blood_chr2[i],Blood_chr2$pig_pos),]
  b<-FM_Blood_hum[match(intersect(a$index,FM_Blood_hum$index),FM_Blood_hum$index),]
  chr2_Blood_hum<-rbind(chr2_Blood_hum,b)
}
chr2_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr2_Blood_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr2_Blood_hum_one2one<-chr2_Blood_hum[match(intersect(chr2_pig2hum_one2one,chr2_Blood_hum$gene_id),chr2_Blood_hum$gene_id),]

chr3_Blood_overloci_one2one<-chr3_Blood_overloci[match(intersect(chr3_Blood_genes,one2one_pig),chr3_Blood_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Blood_chr3$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Blood_chr3$hum_chr<-tmp$V1
Blood_chr3$hum_loci<-tmp$V3
FM_Blood_hum$index<-paste0(FM_Blood_hum$chr,"_",FM_Blood_hum$variant_pos)
Blood_chr3$index<-paste0("chr",Blood_chr3$hum_chr,"_",Blood_chr3$hum_loci)
chr3_Blood_hum<-NULL
for(i in 1:length(same_Blood_chr3)){
  a<-Blood_chr3[match(same_Blood_chr3[i],Blood_chr3$pig_pos),]
  b<-FM_Blood_hum[match(intersect(a$index,FM_Blood_hum$index),FM_Blood_hum$index),]
  chr3_Blood_hum<-rbind(chr3_Blood_hum,b)
}
chr3_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr3_Blood_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr3_Blood_hum_one2one<-chr3_Blood_hum[match(intersect(chr3_pig2hum_one2one,chr3_Blood_hum$gene_id),chr3_Blood_hum$gene_id),]

chr4_Blood_overloci_one2one<-chr4_Blood_overloci[match(intersect(chr4_Blood_genes,one2one_pig),chr4_Blood_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Blood_chr4$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Blood_chr4$hum_chr<-tmp$V1
Blood_chr4$hum_loci<-tmp$V3
FM_Blood_hum$index<-paste0(FM_Blood_hum$chr,"_",FM_Blood_hum$variant_pos)
Blood_chr4$index<-paste0("chr",Blood_chr4$hum_chr,"_",Blood_chr4$hum_loci)
chr4_Blood_hum<-NULL
for(i in 1:length(same_Blood_chr4)){
  a<-Blood_chr4[match(same_Blood_chr4[i],Blood_chr4$pig_pos),]
  b<-FM_Blood_hum[match(intersect(a$index,FM_Blood_hum$index),FM_Blood_hum$index),]
  chr4_Blood_hum<-rbind(chr4_Blood_hum,b)
}
chr4_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr4_Blood_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr4_Blood_hum_one2one<-chr4_Blood_hum[match(intersect(chr4_pig2hum_one2one,chr4_Blood_hum$gene_id),chr4_Blood_hum$gene_id),]

chr5_Blood_overloci_one2one<-chr5_Blood_overloci[match(intersect(chr5_Blood_genes,one2one_pig),chr5_Blood_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Blood_chr5$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Blood_chr5$hum_chr<-tmp$V1
Blood_chr5$hum_loci<-tmp$V3
FM_Blood_hum$index<-paste0(FM_Blood_hum$chr,"_",FM_Blood_hum$variant_pos)
Blood_chr5$index<-paste0("chr",Blood_chr5$hum_chr,"_",Blood_chr5$hum_loci)
chr5_Blood_hum<-NULL
for(i in 1:length(same_Blood_chr5)){
  a<-Blood_chr5[match(same_Blood_chr5[i],Blood_chr5$pig_pos),]
  b<-FM_Blood_hum[match(intersect(a$index,FM_Blood_hum$index),FM_Blood_hum$index),]
  chr5_Blood_hum<-rbind(chr5_Blood_hum,b)
}
chr5_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr5_Blood_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr5_Blood_hum_one2one<-chr5_Blood_hum[match(intersect(chr5_pig2hum_one2one,chr5_Blood_hum$gene_id),chr5_Blood_hum$gene_id),]

chr6_Blood_overloci_one2one<-chr6_Blood_overloci[match(intersect(chr6_Blood_genes,one2one_pig),chr6_Blood_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Blood_chr6$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Blood_chr6$hum_chr<-tmp$V1
Blood_chr6$hum_loci<-tmp$V3
FM_Blood_hum$index<-paste0(FM_Blood_hum$chr,"_",FM_Blood_hum$variant_pos)
Blood_chr6$index<-paste0("chr",Blood_chr6$hum_chr,"_",Blood_chr6$hum_loci)
chr6_Blood_hum<-NULL
for(i in 1:length(same_Blood_chr6)){
  a<-Blood_chr6[match(same_Blood_chr6[i],Blood_chr6$pig_pos),]
  b<-FM_Blood_hum[match(intersect(a$index,FM_Blood_hum$index),FM_Blood_hum$index),]
  chr6_Blood_hum<-rbind(chr6_Blood_hum,b)
}
chr6_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr6_Blood_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr6_Blood_hum_one2one<-chr6_Blood_hum[match(intersect(chr6_pig2hum_one2one,chr6_Blood_hum$gene_id),chr6_Blood_hum$gene_id),]

chr7_Blood_overloci_one2one<-chr7_Blood_overloci[match(intersect(chr7_Blood_genes,one2one_pig),chr7_Blood_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Blood_chr7$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Blood_chr7$hum_chr<-tmp$V1
Blood_chr7$hum_loci<-tmp$V3
FM_Blood_hum$index<-paste0(FM_Blood_hum$chr,"_",FM_Blood_hum$variant_pos)
Blood_chr7$index<-paste0("chr",Blood_chr7$hum_chr,"_",Blood_chr7$hum_loci)
chr7_Blood_hum<-NULL
for(i in 1:length(same_Blood_chr7)){
  a<-Blood_chr7[match(same_Blood_chr7[i],Blood_chr7$pig_pos),]
  b<-FM_Blood_hum[match(intersect(a$index,FM_Blood_hum$index),FM_Blood_hum$index),]
  chr7_Blood_hum<-rbind(chr7_Blood_hum,b)
}
chr7_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr7_Blood_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr7_Blood_hum_one2one<-chr7_Blood_hum[match(intersect(chr7_pig2hum_one2one,chr7_Blood_hum$gene_id),chr7_Blood_hum$gene_id),]

chr8_Blood_overloci_one2one<-chr8_Blood_overloci[match(intersect(chr8_Blood_genes,one2one_pig),chr8_Blood_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Blood_chr8$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Blood_chr8$hum_chr<-tmp$V1
Blood_chr8$hum_loci<-tmp$V3
FM_Blood_hum$index<-paste0(FM_Blood_hum$chr,"_",FM_Blood_hum$variant_pos)
Blood_chr8$index<-paste0("chr",Blood_chr8$hum_chr,"_",Blood_chr8$hum_loci)
chr8_Blood_hum<-NULL
for(i in 1:length(same_Blood_chr8)){
  a<-Blood_chr8[match(same_Blood_chr8[i],Blood_chr8$pig_pos),]
  b<-FM_Blood_hum[match(intersect(a$index,FM_Blood_hum$index),FM_Blood_hum$index),]
  chr8_Blood_hum<-rbind(chr8_Blood_hum,b)
}
chr8_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr8_Blood_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr8_Blood_hum_one2one<-chr8_Blood_hum[match(intersect(chr8_pig2hum_one2one,chr8_Blood_hum$gene_id),chr8_Blood_hum$gene_id),]

chr9_Blood_overloci_one2one<-chr9_Blood_overloci[match(intersect(chr9_Blood_genes,one2one_pig),chr9_Blood_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Blood_chr9$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Blood_chr9$hum_chr<-tmp$V1
Blood_chr9$hum_loci<-tmp$V3
FM_Blood_hum$index<-paste0(FM_Blood_hum$chr,"_",FM_Blood_hum$variant_pos)
Blood_chr9$index<-paste0("chr",Blood_chr9$hum_chr,"_",Blood_chr9$hum_loci)
chr9_Blood_hum<-NULL
for(i in 1:length(same_Blood_chr9)){
  a<-Blood_chr9[match(same_Blood_chr9[i],Blood_chr9$pig_pos),]
  b<-FM_Blood_hum[match(intersect(a$index,FM_Blood_hum$index),FM_Blood_hum$index),]
  chr9_Blood_hum<-rbind(chr9_Blood_hum,b)
}
chr9_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr9_Blood_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr9_Blood_hum_one2one<-chr9_Blood_hum[match(intersect(chr9_pig2hum_one2one,chr9_Blood_hum$gene_id),chr9_Blood_hum$gene_id),]

chr10_Blood_overloci_one2one<-chr10_Blood_overloci[match(intersect(chr10_Blood_genes,one2one_pig),chr10_Blood_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Blood_chr10$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Blood_chr10$hum_chr<-tmp$V1
Blood_chr10$hum_loci<-tmp$V3
FM_Blood_hum$index<-paste0(FM_Blood_hum$chr,"_",FM_Blood_hum$variant_pos)
Blood_chr10$index<-paste0("chr",Blood_chr10$hum_chr,"_",Blood_chr10$hum_loci)
chr10_Blood_hum<-NULL
for(i in 1:length(same_Blood_chr10)){
  a<-Blood_chr10[match(same_Blood_chr10[i],Blood_chr10$pig_pos),]
  b<-FM_Blood_hum[match(intersect(a$index,FM_Blood_hum$index),FM_Blood_hum$index),]
  chr10_Blood_hum<-rbind(chr10_Blood_hum,b)
}
chr10_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr10_Blood_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr10_Blood_hum_one2one<-chr10_Blood_hum[match(intersect(chr10_pig2hum_one2one,chr10_Blood_hum$gene_id),chr10_Blood_hum$gene_id),]

chr11_Blood_overloci_one2one<-chr11_Blood_overloci[match(intersect(chr11_Blood_genes,one2one_pig),chr11_Blood_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Blood_chr11$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Blood_chr11$hum_chr<-tmp$V1
Blood_chr11$hum_loci<-tmp$V3
FM_Blood_hum$index<-paste0(FM_Blood_hum$chr,"_",FM_Blood_hum$variant_pos)
Blood_chr11$index<-paste0("chr",Blood_chr11$hum_chr,"_",Blood_chr11$hum_loci)
chr11_Blood_hum<-NULL
for(i in 1:length(same_Blood_chr11)){
  a<-Blood_chr11[match(same_Blood_chr11[i],Blood_chr11$pig_pos),]
  b<-FM_Blood_hum[match(intersect(a$index,FM_Blood_hum$index),FM_Blood_hum$index),]
  chr11_Blood_hum<-rbind(chr11_Blood_hum,b)
}
chr11_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr11_Blood_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr11_Blood_hum_one2one<-chr11_Blood_hum[match(intersect(chr11_pig2hum_one2one,chr11_Blood_hum$gene_id),chr11_Blood_hum$gene_id),]

chr12_Blood_overloci_one2one<-chr12_Blood_overloci[match(intersect(chr12_Blood_genes,one2one_pig),chr12_Blood_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Blood_chr12$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Blood_chr12$hum_chr<-tmp$V1
Blood_chr12$hum_loci<-tmp$V3
FM_Blood_hum$index<-paste0(FM_Blood_hum$chr,"_",FM_Blood_hum$variant_pos)
Blood_chr12$index<-paste0("chr",Blood_chr12$hum_chr,"_",Blood_chr12$hum_loci)
chr12_Blood_hum<-NULL
for(i in 1:length(same_Blood_chr12)){
  a<-Blood_chr12[match(same_Blood_chr12[i],Blood_chr12$pig_pos),]
  b<-FM_Blood_hum[match(intersect(a$index,FM_Blood_hum$index),FM_Blood_hum$index),]
  chr12_Blood_hum<-rbind(chr12_Blood_hum,b)
}
chr12_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr12_Blood_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr12_Blood_hum_one2one<-chr12_Blood_hum[match(intersect(chr12_pig2hum_one2one,chr12_Blood_hum$gene_id),chr12_Blood_hum$gene_id),]

chr13_Blood_overloci_one2one<-chr13_Blood_overloci[match(intersect(chr13_Blood_genes,one2one_pig),chr13_Blood_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Blood_chr13$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Blood_chr13$hum_chr<-tmp$V1
Blood_chr13$hum_loci<-tmp$V3
FM_Blood_hum$index<-paste0(FM_Blood_hum$chr,"_",FM_Blood_hum$variant_pos)
Blood_chr13$index<-paste0("chr",Blood_chr13$hum_chr,"_",Blood_chr13$hum_loci)
chr13_Blood_hum<-NULL
for(i in 1:length(same_Blood_chr13)){
  a<-Blood_chr13[match(same_Blood_chr13[i],Blood_chr13$pig_pos),]
  b<-FM_Blood_hum[match(intersect(a$index,FM_Blood_hum$index),FM_Blood_hum$index),]
  chr13_Blood_hum<-rbind(chr13_Blood_hum,b)
}
chr13_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr13_Blood_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr13_Blood_hum_one2one<-chr13_Blood_hum[match(intersect(chr13_pig2hum_one2one,chr13_Blood_hum$gene_id),chr13_Blood_hum$gene_id),]

chr14_Blood_overloci_one2one<-chr14_Blood_overloci[match(intersect(chr14_Blood_genes,one2one_pig),chr14_Blood_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Blood_chr14$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Blood_chr14$hum_chr<-tmp$V1
Blood_chr14$hum_loci<-tmp$V3
FM_Blood_hum$index<-paste0(FM_Blood_hum$chr,"_",FM_Blood_hum$variant_pos)
Blood_chr14$index<-paste0("chr",Blood_chr14$hum_chr,"_",Blood_chr14$hum_loci)
chr14_Blood_hum<-NULL
for(i in 1:length(same_Blood_chr14)){
  a<-Blood_chr14[match(same_Blood_chr14[i],Blood_chr14$pig_pos),]
  b<-FM_Blood_hum[match(intersect(a$index,FM_Blood_hum$index),FM_Blood_hum$index),]
  chr14_Blood_hum<-rbind(chr14_Blood_hum,b)
}
chr14_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr14_Blood_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr14_Blood_hum_one2one<-chr14_Blood_hum[match(intersect(chr14_pig2hum_one2one,chr14_Blood_hum$gene_id),chr14_Blood_hum$gene_id),]

chr15_Blood_overloci_one2one<-chr15_Blood_overloci[match(intersect(chr15_Blood_genes,one2one_pig),chr15_Blood_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Blood_chr15$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Blood_chr15$hum_chr<-tmp$V1
Blood_chr15$hum_loci<-tmp$V3
FM_Blood_hum$index<-paste0(FM_Blood_hum$chr,"_",FM_Blood_hum$variant_pos)
Blood_chr15$index<-paste0("chr",Blood_chr15$hum_chr,"_",Blood_chr15$hum_loci)
chr15_Blood_hum<-NULL
for(i in 1:length(same_Blood_chr15)){
  a<-Blood_chr15[match(same_Blood_chr15[i],Blood_chr15$pig_pos),]
  b<-FM_Blood_hum[match(intersect(a$index,FM_Blood_hum$index),FM_Blood_hum$index),]
  chr15_Blood_hum<-rbind(chr15_Blood_hum,b)
}
chr15_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr15_Blood_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr15_Blood_hum_one2one<-chr15_Blood_hum[match(intersect(chr15_pig2hum_one2one,chr15_Blood_hum$gene_id),chr15_Blood_hum$gene_id),]

chr16_Blood_overloci_one2one<-chr16_Blood_overloci[match(intersect(chr16_Blood_genes,one2one_pig),chr16_Blood_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Blood_chr16$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Blood_chr16$hum_chr<-tmp$V1
Blood_chr16$hum_loci<-tmp$V3
FM_Blood_hum$index<-paste0(FM_Blood_hum$chr,"_",FM_Blood_hum$variant_pos)
Blood_chr16$index<-paste0("chr",Blood_chr16$hum_chr,"_",Blood_chr16$hum_loci)
chr16_Blood_hum<-NULL
for(i in 1:length(same_Blood_chr16)){
  a<-Blood_chr16[match(same_Blood_chr16[i],Blood_chr16$pig_pos),]
  b<-FM_Blood_hum[match(intersect(a$index,FM_Blood_hum$index),FM_Blood_hum$index),]
  chr16_Blood_hum<-rbind(chr16_Blood_hum,b)
}
chr16_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr16_Blood_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr16_Blood_hum_one2one<-chr16_Blood_hum[match(intersect(chr16_pig2hum_one2one,chr16_Blood_hum$gene_id),chr16_Blood_hum$gene_id),]

chr17_Blood_overloci_one2one<-chr17_Blood_overloci[match(intersect(chr17_Blood_genes,one2one_pig),chr17_Blood_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Blood_chr17$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Blood_chr17$hum_chr<-tmp$V1
Blood_chr17$hum_loci<-tmp$V3
FM_Blood_hum$index<-paste0(FM_Blood_hum$chr,"_",FM_Blood_hum$variant_pos)
Blood_chr17$index<-paste0("chr",Blood_chr17$hum_chr,"_",Blood_chr17$hum_loci)
chr17_Blood_hum<-NULL
for(i in 1:length(same_Blood_chr17)){
  a<-Blood_chr17[match(same_Blood_chr17[i],Blood_chr17$pig_pos),]
  b<-FM_Blood_hum[match(intersect(a$index,FM_Blood_hum$index),FM_Blood_hum$index),]
  chr17_Blood_hum<-rbind(chr17_Blood_hum,b)
}
chr17_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr17_Blood_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr17_Blood_hum_one2one<-chr17_Blood_hum[match(intersect(chr17_pig2hum_one2one,chr17_Blood_hum$gene_id),chr17_Blood_hum$gene_id),]

chr18_Blood_overloci_one2one<-chr18_Blood_overloci[match(intersect(chr18_Blood_genes,one2one_pig),chr18_Blood_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Blood_chr18$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Blood_chr18$hum_chr<-tmp$V1
Blood_chr18$hum_loci<-tmp$V3
FM_Blood_hum$index<-paste0(FM_Blood_hum$chr,"_",FM_Blood_hum$variant_pos)
Blood_chr18$index<-paste0("chr",Blood_chr18$hum_chr,"_",Blood_chr18$hum_loci)
chr18_Blood_hum<-NULL
for(i in 1:length(same_Blood_chr18)){
  a<-Blood_chr18[match(same_Blood_chr18[i],Blood_chr18$pig_pos),]
  b<-FM_Blood_hum[match(intersect(a$index,FM_Blood_hum$index),FM_Blood_hum$index),]
  chr18_Blood_hum<-rbind(chr18_Blood_hum,b)
}
chr18_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr18_Blood_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr18_Blood_hum_one2one<-chr18_Blood_hum[match(intersect(chr18_pig2hum_one2one,chr18_Blood_hum$gene_id),chr18_Blood_hum$gene_id),]

Blood_one2one_SNP_hum<-rbind(chr1_Blood_hum_one2one,chr2_Blood_hum_one2one,chr3_Blood_hum_one2one,chr4_Blood_hum_one2one,chr5_Blood_hum_one2one,
                             chr6_Blood_hum_one2one,chr7_Blood_hum_one2one,chr8_Blood_hum_one2one,chr9_Blood_hum_one2one,chr10_Blood_hum_one2one,
                             chr11_Blood_hum_one2one,chr12_Blood_hum_one2one,chr13_Blood_hum_one2one,chr14_Blood_hum_one2one,chr15_Blood_hum_one2one,
                             chr16_Blood_hum_one2one,chr17_Blood_hum_one2one,chr18_Blood_hum_one2one)

chr1_Blood_pig_one2one<-chr1_Blood_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr1_Blood_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr1_Blood_overloci_one2one$phenotype_id),1:9]
chr2_Blood_pig_one2one<-chr2_Blood_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr2_Blood_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr2_Blood_overloci_one2one$phenotype_id),1:9]
chr3_Blood_pig_one2one<-chr3_Blood_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr3_Blood_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr3_Blood_overloci_one2one$phenotype_id),1:9]
chr4_Blood_pig_one2one<-chr4_Blood_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr4_Blood_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr4_Blood_overloci_one2one$phenotype_id),1:9]
chr5_Blood_pig_one2one<-chr5_Blood_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr5_Blood_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr5_Blood_overloci_one2one$phenotype_id),1:9]
chr6_Blood_pig_one2one<-chr6_Blood_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr6_Blood_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr6_Blood_overloci_one2one$phenotype_id),1:9]
chr7_Blood_pig_one2one<-chr7_Blood_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr7_Blood_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr7_Blood_overloci_one2one$phenotype_id),1:9]
chr8_Blood_pig_one2one<-chr8_Blood_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr8_Blood_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr8_Blood_overloci_one2one$phenotype_id),1:9]
chr9_Blood_pig_one2one<-chr9_Blood_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr9_Blood_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr9_Blood_overloci_one2one$phenotype_id),1:9]
chr10_Blood_pig_one2one<-chr10_Blood_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr10_Blood_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr10_Blood_overloci_one2one$phenotype_id),1:9]
chr11_Blood_pig_one2one<-chr11_Blood_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr11_Blood_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr11_Blood_overloci_one2one$phenotype_id),1:9]
chr12_Blood_pig_one2one<-chr12_Blood_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr12_Blood_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr12_Blood_overloci_one2one$phenotype_id),1:9]
chr13_Blood_pig_one2one<-chr13_Blood_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr13_Blood_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr13_Blood_overloci_one2one$phenotype_id),1:9]
chr14_Blood_pig_one2one<-chr14_Blood_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr14_Blood_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr14_Blood_overloci_one2one$phenotype_id),1:9]
chr15_Blood_pig_one2one<-chr15_Blood_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr15_Blood_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr15_Blood_overloci_one2one$phenotype_id),1:9]
chr16_Blood_pig_one2one<-chr16_Blood_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr16_Blood_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr16_Blood_overloci_one2one$phenotype_id),1:9]
chr17_Blood_pig_one2one<-chr17_Blood_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr17_Blood_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr17_Blood_overloci_one2one$phenotype_id),1:9]
chr18_Blood_pig_one2one<-chr18_Blood_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr18_Blood_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr18_Blood_overloci_one2one$phenotype_id),1:9]

Blood_one2one_SNP_pig<-rbind(chr1_Blood_pig_one2one,chr2_Blood_pig_one2one,chr3_Blood_pig_one2one,chr4_Blood_pig_one2one,chr5_Blood_pig_one2one,
                             chr6_Blood_pig_one2one,chr7_Blood_pig_one2one,chr8_Blood_pig_one2one,chr9_Blood_pig_one2one,chr10_Blood_pig_one2one,
                             chr11_Blood_pig_one2one,chr12_Blood_pig_one2one,chr13_Blood_pig_one2one,chr14_Blood_pig_one2one,chr15_Blood_pig_one2one,
                             chr16_Blood_pig_one2one,chr17_Blood_pig_one2one,chr18_Blood_pig_one2one)

Blood_SNP_sum<-array(NA,dim=c(nrow(Blood_one2one_SNP_hum),2))
colnames(Blood_SNP_sum)<-c("Human","Pig")
Blood_SNP_sum<-as.data.frame(Blood_SNP_sum)
Blood_SNP_sum$Human<-Blood_one2one_SNP_hum$slope / Blood_one2one_SNP_hum$slope_se
Blood_SNP_sum$Pig<-Blood_one2one_SNP_pig$slope / Blood_one2one_SNP_pig$slope_se
cor<-cor(abs(Blood_SNP_sum$Human),abs(Blood_SNP_sum$Pig))
p_val<-t.test(abs(Blood_SNP_sum$Human),abs(Blood_SNP_sum$Pig))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/Blood_SNP_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(Blood_SNP_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The number of eGenes expressed in common tissues",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

save(Blood_one2one_SNP_hum,Blood_one2one_SNP_pig,file="/Users/baizhonghao/Downloads/human-pig GTEx/Blood_SNP.Rdata")

#Colon_SNP_overlaploci#
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

eqtl_Colon_chr1<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Colon/Colon.cis_qtl_pairs.1.txt"))
eqtl_Colon_chr2<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Colon/Colon.cis_qtl_pairs.2.txt"))
eqtl_Colon_chr3<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Colon/Colon.cis_qtl_pairs.3.txt"))
eqtl_Colon_chr4<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Colon/Colon.cis_qtl_pairs.4.txt"))
eqtl_Colon_chr5<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Colon/Colon.cis_qtl_pairs.5.txt"))
eqtl_Colon_chr6<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Colon/Colon.cis_qtl_pairs.6.txt"))
eqtl_Colon_chr7<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Colon/Colon.cis_qtl_pairs.7.txt"))
eqtl_Colon_chr8<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Colon/Colon.cis_qtl_pairs.8.txt"))
eqtl_Colon_chr9<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Colon/Colon.cis_qtl_pairs.9.txt"))
eqtl_Colon_chr10<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Colon/Colon.cis_qtl_pairs.10.txt"))
eqtl_Colon_chr11<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Colon/Colon.cis_qtl_pairs.11.txt"))
eqtl_Colon_chr12<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Colon/Colon.cis_qtl_pairs.12.txt"))
eqtl_Colon_chr13<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Colon/Colon.cis_qtl_pairs.13.txt"))
eqtl_Colon_chr14<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Colon/Colon.cis_qtl_pairs.14.txt"))
eqtl_Colon_chr15<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Colon/Colon.cis_qtl_pairs.15.txt"))
eqtl_Colon_chr16<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Colon/Colon.cis_qtl_pairs.16.txt"))
eqtl_Colon_chr17<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Colon/Colon.cis_qtl_pairs.17.txt"))
eqtl_Colon_chr18<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Colon/Colon.cis_qtl_pairs.18.txt"))

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon_chr1$variant_id,split="_"))))
eqtl_Colon_chr1$chr<-tmp$V1
eqtl_Colon_chr1$loci<-tmp$V2

eqtl_Colon_chr1$index<-paste0(eqtl_Colon_chr1$chr,"-",eqtl_Colon_chr1$loci)
Colon_chr1$index<-paste0(Colon_chr1$chr,"-",Colon_chr1$pig_pos)
same_Colon_chr1<-intersect(Colon_chr1$pig_pos,eqtl_Colon_chr1$loci)


tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon_chr2$variant_id,split="_"))))
eqtl_Colon_chr2$chr<-tmp$V1
eqtl_Colon_chr2$loci<-tmp$V2

eqtl_Colon_chr2$index<-paste0(eqtl_Colon_chr2$chr,"-",eqtl_Colon_chr2$loci)
Colon_chr2$index<-paste0(Colon_chr2$chr,"-",Colon_chr2$pig_pos)
same_Colon_chr2<-intersect(Colon_chr2$pig_pos,eqtl_Colon_chr2$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon_chr3$variant_id,split="_"))))
eqtl_Colon_chr3$chr<-tmp$V1
eqtl_Colon_chr3$loci<-tmp$V2

eqtl_Colon_chr3$index<-paste0(eqtl_Colon_chr3$chr,"-",eqtl_Colon_chr3$loci)
Colon_chr3$index<-paste0(Colon_chr3$chr,"-",Colon_chr3$pig_pos)
same_Colon_chr3<-intersect(Colon_chr3$pig_pos,eqtl_Colon_chr3$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon_chr4$variant_id,split="_"))))
eqtl_Colon_chr4$chr<-tmp$V1
eqtl_Colon_chr4$loci<-tmp$V2

eqtl_Colon_chr4$index<-paste0(eqtl_Colon_chr4$chr,"-",eqtl_Colon_chr4$loci)
Colon_chr4$index<-paste0(Colon_chr4$chr,"-",Colon_chr4$pig_pos)
same_Colon_chr4<-intersect(Colon_chr4$pig_pos,eqtl_Colon_chr4$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon_chr5$variant_id,split="_"))))
eqtl_Colon_chr5$loci<-tmp$V2
same_Colon_chr5<-intersect(Colon_chr5$pig_pos,eqtl_Colon_chr5$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon_chr6$variant_id,split="_"))))
eqtl_Colon_chr6$loci<-tmp$V2
same_Colon_chr6<-intersect(Colon_chr6$pig_pos,eqtl_Colon_chr6$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon_chr7$variant_id,split="_"))))
eqtl_Colon_chr7$loci<-tmp$V2
same_Colon_chr7<-intersect(Colon_chr7$pig_pos,eqtl_Colon_chr7$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon_chr8$variant_id,split="_"))))
eqtl_Colon_chr8$loci<-tmp$V2
same_Colon_chr8<-intersect(Colon_chr8$pig_pos,eqtl_Colon_chr8$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon_chr9$variant_id,split="_"))))
eqtl_Colon_chr9$loci<-tmp$V2
same_Colon_chr9<-intersect(Colon_chr9$pig_pos,eqtl_Colon_chr9$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon_chr10$variant_id,split="_"))))
eqtl_Colon_chr10$loci<-tmp$V2
same_Colon_chr10<-intersect(Colon_chr10$pig_pos,eqtl_Colon_chr10$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon_chr11$variant_id,split="_"))))
eqtl_Colon_chr11$loci<-tmp$V2
same_Colon_chr11<-intersect(Colon_chr11$pig_pos,eqtl_Colon_chr11$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon_chr12$variant_id,split="_"))))
eqtl_Colon_chr12$loci<-tmp$V2
same_Colon_chr12<-intersect(Colon_chr12$pig_pos,eqtl_Colon_chr12$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon_chr13$variant_id,split="_"))))
eqtl_Colon_chr13$loci<-tmp$V2
same_Colon_chr13<-intersect(Colon_chr13$pig_pos,eqtl_Colon_chr13$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon_chr14$variant_id,split="_"))))
eqtl_Colon_chr14$loci<-tmp$V2
same_Colon_chr14<-intersect(Colon_chr14$pig_pos,eqtl_Colon_chr14$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon_chr15$variant_id,split="_"))))
eqtl_Colon_chr15$loci<-tmp$V2
same_Colon_chr15<-intersect(Colon_chr15$pig_pos,eqtl_Colon_chr15$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon_chr16$variant_id,split="_"))))
eqtl_Colon_chr16$loci<-tmp$V2
same_Colon_chr16<-intersect(Colon_chr16$pig_pos,eqtl_Colon_chr16$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon_chr17$variant_id,split="_"))))
eqtl_Colon_chr17$loci<-tmp$V2
same_Colon_chr17<-intersect(Colon_chr17$pig_pos,eqtl_Colon_chr17$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Colon_chr18$variant_id,split="_"))))
eqtl_Colon_chr18$loci<-tmp$V2
same_Colon_chr18<-intersect(Colon_chr18$pig_pos,eqtl_Colon_chr18$loci)

chr1_Colon_genes<-NULL
chr1_Colon_overloci<-NULL
if(length(same_Colon_chr1!=0)){
  for(i in 1:length(same_Colon_chr1)){
    a<-eqtl_Colon_chr1$phenotype_id[grep(same_Colon_chr1[i],eqtl_Colon_chr1$loci)]
    b<-eqtl_Colon_chr1[grep(same_Colon_chr1[i],eqtl_Colon_chr1$loci),]
    chr1_Colon_genes<-c(chr1_Colon_genes,a)
    chr1_Colon_overloci<-rbind(chr1_Colon_overloci,b)
  }
}

chr2_Colon_genes<-NULL
chr2_Colon_overloci<-NULL
if(length(same_Colon_chr2!=0)){
  for(i in 1:length(same_Colon_chr2)){
    a<-eqtl_Colon_chr2$phenotype_id[grep(same_Colon_chr2[i],eqtl_Colon_chr2$loci)]
    b<-eqtl_Colon_chr2[grep(same_Colon_chr2[i],eqtl_Colon_chr2$loci),]
    chr2_Colon_genes<-c(chr2_Colon_genes,a)
    chr2_Colon_overloci<-rbind(chr2_Colon_overloci,b)
  }
}

chr3_Colon_genes<-NULL
chr3_Colon_overloci<-NULL
if(length(same_Colon_chr3!=0)){
  for(i in 1:length(same_Colon_chr3)){
    a<-eqtl_Colon_chr3$phenotype_id[grep(same_Colon_chr3[i],eqtl_Colon_chr3$loci)]
    b<-eqtl_Colon_chr3[grep(same_Colon_chr3[i],eqtl_Colon_chr3$loci),]
    chr3_Colon_genes<-c(chr3_Colon_genes,a)
    chr3_Colon_overloci<-rbind(chr3_Colon_overloci,b)
  }
}

chr4_Colon_genes<-NULL
chr4_Colon_overloci<-NULL
if(length(same_Colon_chr4!=0)){
  for(i in 1:length(same_Colon_chr4)){
    a<-eqtl_Colon_chr4$phenotype_id[grep(same_Colon_chr4[i],eqtl_Colon_chr4$loci)]
    b<-eqtl_Colon_chr4[grep(same_Colon_chr4[i],eqtl_Colon_chr4$loci),]
    chr4_Colon_genes<-c(chr4_Colon_genes,a)
    chr4_Colon_overloci<-rbind(chr4_Colon_overloci,b)
  }
}

chr5_Colon_genes<-NULL
chr5_Colon_overloci<-NULL
if(length(same_Colon_chr5!=0)){
  for(i in 1:length(same_Colon_chr5)){
    a<-eqtl_Colon_chr5$phenotype_id[grep(same_Colon_chr5[i],eqtl_Colon_chr5$loci)]
    b<-eqtl_Colon_chr5[grep(same_Colon_chr5[i],eqtl_Colon_chr5$loci),]
    chr5_Colon_genes<-c(chr5_Colon_genes,a)
    chr5_Colon_overloci<-rbind(chr5_Colon_overloci,b)
  }
}

chr6_Colon_genes<-NULL
chr6_Colon_overloci<-NULL
if(length(same_Colon_chr6!=0)){
  for(i in 1:length(same_Colon_chr6)){
    a<-eqtl_Colon_chr6$phenotype_id[grep(same_Colon_chr6[i],eqtl_Colon_chr6$loci)]
    b<-eqtl_Colon_chr6[grep(same_Colon_chr6[i],eqtl_Colon_chr6$loci),]
    chr6_Colon_genes<-c(chr6_Colon_genes,a)
    chr6_Colon_overloci<-rbind(chr6_Colon_overloci,b)
  }
}

chr7_Colon_genes<-NULL
chr7_Colon_overloci<-NULL
if(length(same_Colon_chr7!=0)){
  for(i in 1:length(same_Colon_chr7)){
    a<-eqtl_Colon_chr7$phenotype_id[grep(same_Colon_chr7[i],eqtl_Colon_chr7$loci)]
    b<-eqtl_Colon_chr7[grep(same_Colon_chr7[i],eqtl_Colon_chr7$loci),]
    chr7_Colon_genes<-c(chr7_Colon_genes,a)
    chr7_Colon_overloci<-rbind(chr7_Colon_overloci,b)
  }
}

chr8_Colon_genes<-NULL
chr8_Colon_overloci<-NULL
if(length(same_Colon_chr8!=0)){
  for(i in 1:length(same_Colon_chr8)){
    a<-eqtl_Colon_chr8$phenotype_id[grep(same_Colon_chr8[i],eqtl_Colon_chr8$loci)]
    b<-eqtl_Colon_chr8[grep(same_Colon_chr8[i],eqtl_Colon_chr8$loci),]
    chr8_Colon_genes<-c(chr8_Colon_genes,a)
    chr8_Colon_overloci<-rbind(chr8_Colon_overloci,b)
  }
}

if(length(same_Colon_chr9!=0)){
  chr9_Colon_genes<-NULL
  chr9_Colon_overloci<-NULL
  for(i in 1:length(same_Colon_chr9)){
    a<-eqtl_Colon_chr9$phenotype_id[grep(same_Colon_chr9[i],eqtl_Colon_chr9$loci)]
    b<-eqtl_Colon_chr9[grep(same_Colon_chr9[i],eqtl_Colon_chr9$loci),]
    chr9_Colon_genes<-c(chr9_Colon_genes,a)
    chr9_Colon_overloci<-rbind(chr9_Colon_overloci,b)
  }
}

chr10_Colon_genes<-NULL
chr10_Colon_overloci<-NULL
if(length(same_Colon_chr10!=0)){
  for(i in 1:length(same_Colon_chr10)){
    a<-eqtl_Colon_chr10$phenotype_id[grep(same_Colon_chr10[i],eqtl_Colon_chr10$loci)]
    b<-eqtl_Colon_chr10[grep(same_Colon_chr10[i],eqtl_Colon_chr10$loci),]
    chr10_Colon_genes<-c(chr10_Colon_genes,a)
    chr10_Colon_overloci<-rbind(chr10_Colon_overloci,b)
  }
}

chr11_Colon_genes<-NULL
chr11_Colon_overloci<-NULL
if(length(same_Colon_chr11!=0)){
  for(i in 1:length(same_Colon_chr11)){
    a<-eqtl_Colon_chr11$phenotype_id[grep(same_Colon_chr11[i],eqtl_Colon_chr11$loci)]
    b<-eqtl_Colon_chr11[grep(same_Colon_chr11[i],eqtl_Colon_chr11$loci),]
    chr11_Colon_genes<-c(chr11_Colon_genes,a)
    chr11_Colon_overloci<-rbind(chr11_Colon_overloci,b)
  }
}

chr12_Colon_genes<-NULL
chr12_Colon_overloci<-NULL
if(length(same_Colon_chr12!=0)){
  for(i in 1:length(same_Colon_chr12)){
    a<-eqtl_Colon_chr12$phenotype_id[grep(same_Colon_chr12[i],eqtl_Colon_chr12$loci)]
    b<-eqtl_Colon_chr12[grep(same_Colon_chr12[i],eqtl_Colon_chr12$loci),]
    chr12_Colon_genes<-c(chr12_Colon_genes,a)
    chr12_Colon_overloci<-rbind(chr12_Colon_overloci,b)
  }
}

chr13_Colon_genes<-NULL
chr13_Colon_overloci<-NULL
if(length(same_Colon_chr13!=0)){
  for(i in 1:length(same_Colon_chr13)){
    a<-eqtl_Colon_chr13$phenotype_id[grep(same_Colon_chr13[i],eqtl_Colon_chr13$loci)]
    b<-eqtl_Colon_chr13[grep(same_Colon_chr13[i],eqtl_Colon_chr13$loci),]
    chr13_Colon_genes<-c(chr13_Colon_genes,a)
    chr13_Colon_overloci<-rbind(chr13_Colon_overloci,b)
  }
}

chr14_Colon_genes<-NULL
chr14_Colon_overloci<-NULL
if(length(same_Colon_chr14!=0)){
  for(i in 1:length(same_Colon_chr14)){
    a<-eqtl_Colon_chr14$phenotype_id[grep(same_Colon_chr14[i],eqtl_Colon_chr14$loci)]
    b<-eqtl_Colon_chr14[grep(same_Colon_chr14[i],eqtl_Colon_chr14$loci),]
    chr14_Colon_genes<-c(chr14_Colon_genes,a)
    chr14_Colon_overloci<-rbind(chr14_Colon_overloci,b)
  }
}

chr15_Colon_genes<-NULL
chr15_Colon_overloci<-NULL
if(length(same_Colon_chr15!=0)){
  for(i in 1:length(same_Colon_chr15)){
    a<-eqtl_Colon_chr15$phenotype_id[grep(same_Colon_chr15[i],eqtl_Colon_chr15$loci)]
    b<-eqtl_Colon_chr15[grep(same_Colon_chr15[i],eqtl_Colon_chr15$loci),]
    chr15_Colon_genes<-c(chr15_Colon_genes,a)
    chr15_Colon_overloci<-rbind(chr15_Colon_overloci,b)
  }
}

chr16_Colon_genes<-NULL
chr16_Colon_overloci<-NULL
if(length(same_Colon_chr16!=0)){
  for(i in 1:length(same_Colon_chr16)){
    a<-eqtl_Colon_chr16$phenotype_id[grep(same_Colon_chr16[i],eqtl_Colon_chr16$loci)]
    b<-eqtl_Colon_chr16[grep(same_Colon_chr16[i],eqtl_Colon_chr16$loci),]
    chr16_Colon_genes<-c(chr16_Colon_genes,a)
    chr16_Colon_overloci<-rbind(chr16_Colon_overloci,b)
  }
}

chr17_Colon_genes<-NULL
chr17_Colon_overloci<-NULL
if(length(same_Colon_chr17!=0)){
  for(i in 1:length(same_Colon_chr17)){
    a<-eqtl_Colon_chr17$phenotype_id[grep(same_Colon_chr17[i],eqtl_Colon_chr17$loci)]
    b<-eqtl_Colon_chr17[grep(same_Colon_chr17[i],eqtl_Colon_chr17$loci),]
    chr17_Colon_genes<-c(chr17_Colon_genes,a)
    chr17_Colon_overloci<-rbind(chr17_Colon_overloci,b)
  }
}

chr18_Colon_genes<-NULL
chr18_Colon_overloci<-NULL
if(length(same_Colon_chr18!=0)){
  for(i in 1:length(same_Colon_chr18)){
    a<-eqtl_Colon_chr18$phenotype_id[grep(same_Colon_chr18[i],eqtl_Colon_chr18$loci)]
    b<-eqtl_Colon_chr18[grep(same_Colon_chr18[i],eqtl_Colon_chr18$loci),]
    chr18_Colon_genes<-c(chr18_Colon_genes,a)
    chr18_Colon_overloci<-rbind(chr18_Colon_overloci,b)
  }
}
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
one2one_pig<-annotation$`Pig gene stable ID`[which(annotation$`Pig homology type`=="ortholog_one2one")]

chr1_Colon_overloci_one2one<-chr1_Colon_overloci[match(intersect(chr1_Colon_genes,one2one_pig),chr1_Colon_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Colon_chr1$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Colon_chr1$hum_chr<-tmp$V1
Colon_chr1$hum_loci<-tmp$V3
FM_Colon_hum$index<-paste0(FM_Colon_hum$chr,"_",FM_Colon_hum$variant_pos)
Colon_chr1$index<-paste0("chr",Colon_chr1$hum_chr,"_",Colon_chr1$hum_loci)
chr1_Colon_hum<-NULL
for(i in 1:length(same_Colon_chr1)){
  a<-Colon_chr1[match(same_Colon_chr1[i],Colon_chr1$pig_pos),]
  b<-FM_Colon_hum[match(intersect(a$index,FM_Colon_hum$index),FM_Colon_hum$index),]
  chr1_Colon_hum<-rbind(chr1_Colon_hum,b)
}
chr1_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr1_Colon_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr1_Colon_hum_one2one<-chr1_Colon_hum[match(intersect(chr1_pig2hum_one2one,chr1_Colon_hum$gene_id),chr1_Colon_hum$gene_id),]

chr2_Colon_overloci_one2one<-chr2_Colon_overloci[match(intersect(chr2_Colon_genes,one2one_pig),chr2_Colon_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Colon_chr2$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Colon_chr2$hum_chr<-tmp$V1
Colon_chr2$hum_loci<-tmp$V3
FM_Colon_hum$index<-paste0(FM_Colon_hum$chr,"_",FM_Colon_hum$variant_pos)
Colon_chr2$index<-paste0("chr",Colon_chr2$hum_chr,"_",Colon_chr2$hum_loci)
chr2_Colon_hum<-NULL
for(i in 1:length(same_Colon_chr2)){
  a<-Colon_chr2[match(same_Colon_chr2[i],Colon_chr2$pig_pos),]
  b<-FM_Colon_hum[match(intersect(a$index,FM_Colon_hum$index),FM_Colon_hum$index),]
  chr2_Colon_hum<-rbind(chr2_Colon_hum,b)
}
chr2_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr2_Colon_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr2_Colon_hum_one2one<-chr2_Colon_hum[match(intersect(chr2_pig2hum_one2one,chr2_Colon_hum$gene_id),chr2_Colon_hum$gene_id),]

chr3_Colon_overloci_one2one<-chr3_Colon_overloci[match(intersect(chr3_Colon_genes,one2one_pig),chr3_Colon_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Colon_chr3$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Colon_chr3$hum_chr<-tmp$V1
Colon_chr3$hum_loci<-tmp$V3
FM_Colon_hum$index<-paste0(FM_Colon_hum$chr,"_",FM_Colon_hum$variant_pos)
Colon_chr3$index<-paste0("chr",Colon_chr3$hum_chr,"_",Colon_chr3$hum_loci)
chr3_Colon_hum<-NULL
for(i in 1:length(same_Colon_chr3)){
  a<-Colon_chr3[match(same_Colon_chr3[i],Colon_chr3$pig_pos),]
  b<-FM_Colon_hum[match(intersect(a$index,FM_Colon_hum$index),FM_Colon_hum$index),]
  chr3_Colon_hum<-rbind(chr3_Colon_hum,b)
}
chr3_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr3_Colon_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr3_Colon_hum_one2one<-chr3_Colon_hum[match(intersect(chr3_pig2hum_one2one,chr3_Colon_hum$gene_id),chr3_Colon_hum$gene_id),]

chr4_Colon_overloci_one2one<-chr4_Colon_overloci[match(intersect(chr4_Colon_genes,one2one_pig),chr4_Colon_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Colon_chr4$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Colon_chr4$hum_chr<-tmp$V1
Colon_chr4$hum_loci<-tmp$V3
FM_Colon_hum$index<-paste0(FM_Colon_hum$chr,"_",FM_Colon_hum$variant_pos)
Colon_chr4$index<-paste0("chr",Colon_chr4$hum_chr,"_",Colon_chr4$hum_loci)
chr4_Colon_hum<-NULL
for(i in 1:length(same_Colon_chr4)){
  a<-Colon_chr4[match(same_Colon_chr4[i],Colon_chr4$pig_pos),]
  b<-FM_Colon_hum[match(intersect(a$index,FM_Colon_hum$index),FM_Colon_hum$index),]
  chr4_Colon_hum<-rbind(chr4_Colon_hum,b)
}
chr4_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr4_Colon_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr4_Colon_hum_one2one<-chr4_Colon_hum[match(intersect(chr4_pig2hum_one2one,chr4_Colon_hum$gene_id),chr4_Colon_hum$gene_id),]

chr5_Colon_overloci_one2one<-chr5_Colon_overloci[match(intersect(chr5_Colon_genes,one2one_pig),chr5_Colon_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Colon_chr5$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Colon_chr5$hum_chr<-tmp$V1
Colon_chr5$hum_loci<-tmp$V3
FM_Colon_hum$index<-paste0(FM_Colon_hum$chr,"_",FM_Colon_hum$variant_pos)
Colon_chr5$index<-paste0("chr",Colon_chr5$hum_chr,"_",Colon_chr5$hum_loci)
chr5_Colon_hum<-NULL
for(i in 1:length(same_Colon_chr5)){
  a<-Colon_chr5[match(same_Colon_chr5[i],Colon_chr5$pig_pos),]
  b<-FM_Colon_hum[match(intersect(a$index,FM_Colon_hum$index),FM_Colon_hum$index),]
  chr5_Colon_hum<-rbind(chr5_Colon_hum,b)
}
chr5_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr5_Colon_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr5_Colon_hum_one2one<-chr5_Colon_hum[match(intersect(chr5_pig2hum_one2one,chr5_Colon_hum$gene_id),chr5_Colon_hum$gene_id),]

chr6_Colon_overloci_one2one<-chr6_Colon_overloci[match(intersect(chr6_Colon_genes,one2one_pig),chr6_Colon_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Colon_chr6$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Colon_chr6$hum_chr<-tmp$V1
Colon_chr6$hum_loci<-tmp$V3
FM_Colon_hum$index<-paste0(FM_Colon_hum$chr,"_",FM_Colon_hum$variant_pos)
Colon_chr6$index<-paste0("chr",Colon_chr6$hum_chr,"_",Colon_chr6$hum_loci)
chr6_Colon_hum<-NULL
for(i in 1:length(same_Colon_chr6)){
  a<-Colon_chr6[match(same_Colon_chr6[i],Colon_chr6$pig_pos),]
  b<-FM_Colon_hum[match(intersect(a$index,FM_Colon_hum$index),FM_Colon_hum$index),]
  chr6_Colon_hum<-rbind(chr6_Colon_hum,b)
}
chr6_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr6_Colon_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr6_Colon_hum_one2one<-chr6_Colon_hum[match(intersect(chr6_pig2hum_one2one,chr6_Colon_hum$gene_id),chr6_Colon_hum$gene_id),]

chr7_Colon_overloci_one2one<-chr7_Colon_overloci[match(intersect(chr7_Colon_genes,one2one_pig),chr7_Colon_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Colon_chr7$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Colon_chr7$hum_chr<-tmp$V1
Colon_chr7$hum_loci<-tmp$V3
FM_Colon_hum$index<-paste0(FM_Colon_hum$chr,"_",FM_Colon_hum$variant_pos)
Colon_chr7$index<-paste0("chr",Colon_chr7$hum_chr,"_",Colon_chr7$hum_loci)
chr7_Colon_hum<-NULL
for(i in 1:length(same_Colon_chr7)){
  a<-Colon_chr7[match(same_Colon_chr7[i],Colon_chr7$pig_pos),]
  b<-FM_Colon_hum[match(intersect(a$index,FM_Colon_hum$index),FM_Colon_hum$index),]
  chr7_Colon_hum<-rbind(chr7_Colon_hum,b)
}
chr7_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr7_Colon_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr7_Colon_hum_one2one<-chr7_Colon_hum[match(intersect(chr7_pig2hum_one2one,chr7_Colon_hum$gene_id),chr7_Colon_hum$gene_id),]

chr8_Colon_overloci_one2one<-chr8_Colon_overloci[match(intersect(chr8_Colon_genes,one2one_pig),chr8_Colon_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Colon_chr8$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Colon_chr8$hum_chr<-tmp$V1
Colon_chr8$hum_loci<-tmp$V3
FM_Colon_hum$index<-paste0(FM_Colon_hum$chr,"_",FM_Colon_hum$variant_pos)
Colon_chr8$index<-paste0("chr",Colon_chr8$hum_chr,"_",Colon_chr8$hum_loci)
chr8_Colon_hum<-NULL
for(i in 1:length(same_Colon_chr8)){
  a<-Colon_chr8[match(same_Colon_chr8[i],Colon_chr8$pig_pos),]
  b<-FM_Colon_hum[match(intersect(a$index,FM_Colon_hum$index),FM_Colon_hum$index),]
  chr8_Colon_hum<-rbind(chr8_Colon_hum,b)
}
chr8_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr8_Colon_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr8_Colon_hum_one2one<-chr8_Colon_hum[match(intersect(chr8_pig2hum_one2one,chr8_Colon_hum$gene_id),chr8_Colon_hum$gene_id),]

chr9_Colon_overloci_one2one<-chr9_Colon_overloci[match(intersect(chr9_Colon_genes,one2one_pig),chr9_Colon_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Colon_chr9$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Colon_chr9$hum_chr<-tmp$V1
Colon_chr9$hum_loci<-tmp$V3
FM_Colon_hum$index<-paste0(FM_Colon_hum$chr,"_",FM_Colon_hum$variant_pos)
Colon_chr9$index<-paste0("chr",Colon_chr9$hum_chr,"_",Colon_chr9$hum_loci)
chr9_Colon_hum<-NULL
for(i in 1:length(same_Colon_chr9)){
  a<-Colon_chr9[match(same_Colon_chr9[i],Colon_chr9$pig_pos),]
  b<-FM_Colon_hum[match(intersect(a$index,FM_Colon_hum$index),FM_Colon_hum$index),]
  chr9_Colon_hum<-rbind(chr9_Colon_hum,b)
}
chr9_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr9_Colon_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr9_Colon_hum_one2one<-chr9_Colon_hum[match(intersect(chr9_pig2hum_one2one,chr9_Colon_hum$gene_id),chr9_Colon_hum$gene_id),]

chr10_Colon_overloci_one2one<-chr10_Colon_overloci[match(intersect(chr10_Colon_genes,one2one_pig),chr10_Colon_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Colon_chr10$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Colon_chr10$hum_chr<-tmp$V1
Colon_chr10$hum_loci<-tmp$V3
FM_Colon_hum$index<-paste0(FM_Colon_hum$chr,"_",FM_Colon_hum$variant_pos)
Colon_chr10$index<-paste0("chr",Colon_chr10$hum_chr,"_",Colon_chr10$hum_loci)
chr10_Colon_hum<-NULL
for(i in 1:length(same_Colon_chr10)){
  a<-Colon_chr10[match(same_Colon_chr10[i],Colon_chr10$pig_pos),]
  b<-FM_Colon_hum[match(intersect(a$index,FM_Colon_hum$index),FM_Colon_hum$index),]
  chr10_Colon_hum<-rbind(chr10_Colon_hum,b)
}
chr10_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr10_Colon_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr10_Colon_hum_one2one<-chr10_Colon_hum[match(intersect(chr10_pig2hum_one2one,chr10_Colon_hum$gene_id),chr10_Colon_hum$gene_id),]

chr11_Colon_overloci_one2one<-chr11_Colon_overloci[match(intersect(chr11_Colon_genes,one2one_pig),chr11_Colon_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Colon_chr11$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Colon_chr11$hum_chr<-tmp$V1
Colon_chr11$hum_loci<-tmp$V3
FM_Colon_hum$index<-paste0(FM_Colon_hum$chr,"_",FM_Colon_hum$variant_pos)
Colon_chr11$index<-paste0("chr",Colon_chr11$hum_chr,"_",Colon_chr11$hum_loci)
chr11_Colon_hum<-NULL
for(i in 1:length(same_Colon_chr11)){
  a<-Colon_chr11[match(same_Colon_chr11[i],Colon_chr11$pig_pos),]
  b<-FM_Colon_hum[match(intersect(a$index,FM_Colon_hum$index),FM_Colon_hum$index),]
  chr11_Colon_hum<-rbind(chr11_Colon_hum,b)
}
chr11_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr11_Colon_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr11_Colon_hum_one2one<-chr11_Colon_hum[match(intersect(chr11_pig2hum_one2one,chr11_Colon_hum$gene_id),chr11_Colon_hum$gene_id),]

chr12_Colon_overloci_one2one<-chr12_Colon_overloci[match(intersect(chr12_Colon_genes,one2one_pig),chr12_Colon_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Colon_chr12$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Colon_chr12$hum_chr<-tmp$V1
Colon_chr12$hum_loci<-tmp$V3
FM_Colon_hum$index<-paste0(FM_Colon_hum$chr,"_",FM_Colon_hum$variant_pos)
Colon_chr12$index<-paste0("chr",Colon_chr12$hum_chr,"_",Colon_chr12$hum_loci)
chr12_Colon_hum<-NULL
for(i in 1:length(same_Colon_chr12)){
  a<-Colon_chr12[match(same_Colon_chr12[i],Colon_chr12$pig_pos),]
  b<-FM_Colon_hum[match(intersect(a$index,FM_Colon_hum$index),FM_Colon_hum$index),]
  chr12_Colon_hum<-rbind(chr12_Colon_hum,b)
}
chr12_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr12_Colon_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr12_Colon_hum_one2one<-chr12_Colon_hum[match(intersect(chr12_pig2hum_one2one,chr12_Colon_hum$gene_id),chr12_Colon_hum$gene_id),]

chr13_Colon_overloci_one2one<-chr13_Colon_overloci[match(intersect(chr13_Colon_genes,one2one_pig),chr13_Colon_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Colon_chr13$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Colon_chr13$hum_chr<-tmp$V1
Colon_chr13$hum_loci<-tmp$V3
FM_Colon_hum$index<-paste0(FM_Colon_hum$chr,"_",FM_Colon_hum$variant_pos)
Colon_chr13$index<-paste0("chr",Colon_chr13$hum_chr,"_",Colon_chr13$hum_loci)
chr13_Colon_hum<-NULL
for(i in 1:length(same_Colon_chr13)){
  a<-Colon_chr13[match(same_Colon_chr13[i],Colon_chr13$pig_pos),]
  b<-FM_Colon_hum[match(intersect(a$index,FM_Colon_hum$index),FM_Colon_hum$index),]
  chr13_Colon_hum<-rbind(chr13_Colon_hum,b)
}
chr13_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr13_Colon_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr13_Colon_hum_one2one<-chr13_Colon_hum[match(intersect(chr13_pig2hum_one2one,chr13_Colon_hum$gene_id),chr13_Colon_hum$gene_id),]

chr14_Colon_overloci_one2one<-chr14_Colon_overloci[match(intersect(chr14_Colon_genes,one2one_pig),chr14_Colon_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Colon_chr14$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Colon_chr14$hum_chr<-tmp$V1
Colon_chr14$hum_loci<-tmp$V3
FM_Colon_hum$index<-paste0(FM_Colon_hum$chr,"_",FM_Colon_hum$variant_pos)
Colon_chr14$index<-paste0("chr",Colon_chr14$hum_chr,"_",Colon_chr14$hum_loci)
chr14_Colon_hum<-NULL
for(i in 1:length(same_Colon_chr14)){
  a<-Colon_chr14[match(same_Colon_chr14[i],Colon_chr14$pig_pos),]
  b<-FM_Colon_hum[match(intersect(a$index,FM_Colon_hum$index),FM_Colon_hum$index),]
  chr14_Colon_hum<-rbind(chr14_Colon_hum,b)
}
chr14_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr14_Colon_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr14_Colon_hum_one2one<-chr14_Colon_hum[match(intersect(chr14_pig2hum_one2one,chr14_Colon_hum$gene_id),chr14_Colon_hum$gene_id),]

chr15_Colon_overloci_one2one<-chr15_Colon_overloci[match(intersect(chr15_Colon_genes,one2one_pig),chr15_Colon_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Colon_chr15$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Colon_chr15$hum_chr<-tmp$V1
Colon_chr15$hum_loci<-tmp$V3
FM_Colon_hum$index<-paste0(FM_Colon_hum$chr,"_",FM_Colon_hum$variant_pos)
Colon_chr15$index<-paste0("chr",Colon_chr15$hum_chr,"_",Colon_chr15$hum_loci)
chr15_Colon_hum<-NULL
for(i in 1:length(same_Colon_chr15)){
  a<-Colon_chr15[match(same_Colon_chr15[i],Colon_chr15$pig_pos),]
  b<-FM_Colon_hum[match(intersect(a$index,FM_Colon_hum$index),FM_Colon_hum$index),]
  chr15_Colon_hum<-rbind(chr15_Colon_hum,b)
}
chr15_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr15_Colon_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr15_Colon_hum_one2one<-chr15_Colon_hum[match(intersect(chr15_pig2hum_one2one,chr15_Colon_hum$gene_id),chr15_Colon_hum$gene_id),]

chr16_Colon_overloci_one2one<-chr16_Colon_overloci[match(intersect(chr16_Colon_genes,one2one_pig),chr16_Colon_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Colon_chr16$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Colon_chr16$hum_chr<-tmp$V1
Colon_chr16$hum_loci<-tmp$V3
FM_Colon_hum$index<-paste0(FM_Colon_hum$chr,"_",FM_Colon_hum$variant_pos)
Colon_chr16$index<-paste0("chr",Colon_chr16$hum_chr,"_",Colon_chr16$hum_loci)
chr16_Colon_hum<-NULL
for(i in 1:length(same_Colon_chr16)){
  a<-Colon_chr16[match(same_Colon_chr16[i],Colon_chr16$pig_pos),]
  b<-FM_Colon_hum[match(intersect(a$index,FM_Colon_hum$index),FM_Colon_hum$index),]
  chr16_Colon_hum<-rbind(chr16_Colon_hum,b)
}
chr16_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr16_Colon_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr16_Colon_hum_one2one<-chr16_Colon_hum[match(intersect(chr16_pig2hum_one2one,chr16_Colon_hum$gene_id),chr16_Colon_hum$gene_id),]

chr17_Colon_overloci_one2one<-chr17_Colon_overloci[match(intersect(chr17_Colon_genes,one2one_pig),chr17_Colon_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Colon_chr17$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Colon_chr17$hum_chr<-tmp$V1
Colon_chr17$hum_loci<-tmp$V3
FM_Colon_hum$index<-paste0(FM_Colon_hum$chr,"_",FM_Colon_hum$variant_pos)
Colon_chr17$index<-paste0("chr",Colon_chr17$hum_chr,"_",Colon_chr17$hum_loci)
chr17_Colon_hum<-NULL
for(i in 1:length(same_Colon_chr17)){
  a<-Colon_chr17[match(same_Colon_chr17[i],Colon_chr17$pig_pos),]
  b<-FM_Colon_hum[match(intersect(a$index,FM_Colon_hum$index),FM_Colon_hum$index),]
  chr17_Colon_hum<-rbind(chr17_Colon_hum,b)
}
chr17_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr17_Colon_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr17_Colon_hum_one2one<-chr17_Colon_hum[match(intersect(chr17_pig2hum_one2one,chr17_Colon_hum$gene_id),chr17_Colon_hum$gene_id),]

chr18_Colon_overloci_one2one<-chr18_Colon_overloci[match(intersect(chr18_Colon_genes,one2one_pig),chr18_Colon_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Colon_chr18$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Colon_chr18$hum_chr<-tmp$V1
Colon_chr18$hum_loci<-tmp$V3
FM_Colon_hum$index<-paste0(FM_Colon_hum$chr,"_",FM_Colon_hum$variant_pos)
Colon_chr18$index<-paste0("chr",Colon_chr18$hum_chr,"_",Colon_chr18$hum_loci)
chr18_Colon_hum<-NULL
for(i in 1:length(same_Colon_chr18)){
  a<-Colon_chr18[match(same_Colon_chr18[i],Colon_chr18$pig_pos),]
  b<-FM_Colon_hum[match(intersect(a$index,FM_Colon_hum$index),FM_Colon_hum$index),]
  chr18_Colon_hum<-rbind(chr18_Colon_hum,b)
}
chr18_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr18_Colon_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr18_Colon_hum_one2one<-chr18_Colon_hum[match(intersect(chr18_pig2hum_one2one,chr18_Colon_hum$gene_id),chr18_Colon_hum$gene_id),]

Colon_one2one_SNP_hum<-rbind(chr1_Colon_hum_one2one,chr2_Colon_hum_one2one,chr3_Colon_hum_one2one,chr4_Colon_hum_one2one,chr5_Colon_hum_one2one,
                             chr6_Colon_hum_one2one,chr7_Colon_hum_one2one,chr8_Colon_hum_one2one,chr9_Colon_hum_one2one,chr10_Colon_hum_one2one,
                             chr11_Colon_hum_one2one,chr12_Colon_hum_one2one,chr13_Colon_hum_one2one,chr14_Colon_hum_one2one,chr15_Colon_hum_one2one,
                             chr16_Colon_hum_one2one,chr17_Colon_hum_one2one,chr18_Colon_hum_one2one)

chr1_Colon_pig_one2one<-chr1_Colon_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr1_Colon_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr1_Colon_overloci_one2one$phenotype_id),1:9]
chr2_Colon_pig_one2one<-chr2_Colon_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr2_Colon_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr2_Colon_overloci_one2one$phenotype_id),1:9]
chr3_Colon_pig_one2one<-chr3_Colon_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr3_Colon_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr3_Colon_overloci_one2one$phenotype_id),1:9]
chr4_Colon_pig_one2one<-chr4_Colon_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr4_Colon_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr4_Colon_overloci_one2one$phenotype_id),1:9]
chr5_Colon_pig_one2one<-chr5_Colon_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr5_Colon_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr5_Colon_overloci_one2one$phenotype_id),1:9]
chr6_Colon_pig_one2one<-chr6_Colon_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr6_Colon_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr6_Colon_overloci_one2one$phenotype_id),1:9]
chr7_Colon_pig_one2one<-chr7_Colon_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr7_Colon_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr7_Colon_overloci_one2one$phenotype_id),1:9]
chr8_Colon_pig_one2one<-chr8_Colon_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr8_Colon_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr8_Colon_overloci_one2one$phenotype_id),1:9]
chr9_Colon_pig_one2one<-chr9_Colon_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr9_Colon_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr9_Colon_overloci_one2one$phenotype_id),1:9]
chr10_Colon_pig_one2one<-chr10_Colon_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr10_Colon_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr10_Colon_overloci_one2one$phenotype_id),1:9]
chr11_Colon_pig_one2one<-chr11_Colon_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr11_Colon_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr11_Colon_overloci_one2one$phenotype_id),1:9]
chr12_Colon_pig_one2one<-chr12_Colon_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr12_Colon_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr12_Colon_overloci_one2one$phenotype_id),1:9]
chr13_Colon_pig_one2one<-chr13_Colon_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr13_Colon_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr13_Colon_overloci_one2one$phenotype_id),1:9]
chr14_Colon_pig_one2one<-chr14_Colon_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr14_Colon_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr14_Colon_overloci_one2one$phenotype_id),1:9]
chr15_Colon_pig_one2one<-chr15_Colon_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr15_Colon_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr15_Colon_overloci_one2one$phenotype_id),1:9]
chr16_Colon_pig_one2one<-chr16_Colon_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr16_Colon_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr16_Colon_overloci_one2one$phenotype_id),1:9]
chr17_Colon_pig_one2one<-chr17_Colon_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr17_Colon_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr17_Colon_overloci_one2one$phenotype_id),1:9]
chr18_Colon_pig_one2one<-chr18_Colon_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr18_Colon_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr18_Colon_overloci_one2one$phenotype_id),1:9]

Colon_one2one_SNP_pig<-rbind(chr1_Colon_pig_one2one,chr2_Colon_pig_one2one,chr3_Colon_pig_one2one,chr4_Colon_pig_one2one,chr5_Colon_pig_one2one,
                             chr6_Colon_pig_one2one,chr7_Colon_pig_one2one,chr8_Colon_pig_one2one,chr9_Colon_pig_one2one,chr10_Colon_pig_one2one,
                             chr11_Colon_pig_one2one,chr12_Colon_pig_one2one,chr13_Colon_pig_one2one,chr14_Colon_pig_one2one,chr15_Colon_pig_one2one,
                             chr16_Colon_pig_one2one,chr17_Colon_pig_one2one,chr18_Colon_pig_one2one)

Colon_SNP_sum<-array(NA,dim=c(nrow(Colon_one2one_SNP_hum),2))
colnames(Colon_SNP_sum)<-c("Human","Pig")
Colon_SNP_sum<-as.data.frame(Colon_SNP_sum)
Colon_SNP_sum$Human<-Colon_one2one_SNP_hum$slope / Colon_one2one_SNP_hum$slope_se
Colon_SNP_sum$Pig<-Colon_one2one_SNP_pig$slope / Colon_one2one_SNP_pig$slope_se
cor<-cor(abs(Colon_SNP_sum$Human),abs(Colon_SNP_sum$Pig))
p_val<-t.test(abs(Colon_SNP_sum$Human),abs(Colon_SNP_sum$Pig))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/Colon_SNP_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(Colon_SNP_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The number of eGenes expressed in common tissues",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

save(Colon_one2one_SNP_hum,Colon_one2one_SNP_pig,file="/Users/baizhonghao/Downloads/human-pig GTEx/Colon_SNP.Rdata")

#Frontal_cortex_SNP_overlaploci#
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

eqtl_Frontal_cortex_chr1<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Frontal_cortex/Frontal_cortex.cis_qtl_pairs.1.txt"))
eqtl_Frontal_cortex_chr2<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Frontal_cortex/Frontal_cortex.cis_qtl_pairs.2.txt"))
eqtl_Frontal_cortex_chr3<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Frontal_cortex/Frontal_cortex.cis_qtl_pairs.3.txt"))
eqtl_Frontal_cortex_chr4<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Frontal_cortex/Frontal_cortex.cis_qtl_pairs.4.txt"))
eqtl_Frontal_cortex_chr5<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Frontal_cortex/Frontal_cortex.cis_qtl_pairs.5.txt"))
eqtl_Frontal_cortex_chr6<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Frontal_cortex/Frontal_cortex.cis_qtl_pairs.6.txt"))
eqtl_Frontal_cortex_chr7<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Frontal_cortex/Frontal_cortex.cis_qtl_pairs.7.txt"))
eqtl_Frontal_cortex_chr8<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Frontal_cortex/Frontal_cortex.cis_qtl_pairs.8.txt"))
eqtl_Frontal_cortex_chr9<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Frontal_cortex/Frontal_cortex.cis_qtl_pairs.9.txt"))
eqtl_Frontal_cortex_chr10<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Frontal_cortex/Frontal_cortex.cis_qtl_pairs.10.txt"))
eqtl_Frontal_cortex_chr11<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Frontal_cortex/Frontal_cortex.cis_qtl_pairs.11.txt"))
eqtl_Frontal_cortex_chr12<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Frontal_cortex/Frontal_cortex.cis_qtl_pairs.12.txt"))
eqtl_Frontal_cortex_chr13<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Frontal_cortex/Frontal_cortex.cis_qtl_pairs.13.txt"))
eqtl_Frontal_cortex_chr14<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Frontal_cortex/Frontal_cortex.cis_qtl_pairs.14.txt"))
eqtl_Frontal_cortex_chr15<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Frontal_cortex/Frontal_cortex.cis_qtl_pairs.15.txt"))
eqtl_Frontal_cortex_chr16<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Frontal_cortex/Frontal_cortex.cis_qtl_pairs.16.txt"))
eqtl_Frontal_cortex_chr17<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Frontal_cortex/Frontal_cortex.cis_qtl_pairs.17.txt"))
eqtl_Frontal_cortex_chr18<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Frontal_cortex/Frontal_cortex.cis_qtl_pairs.18.txt"))

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex_chr1$variant_id,split="_"))))
eqtl_Frontal_cortex_chr1$chr<-tmp$V1
eqtl_Frontal_cortex_chr1$loci<-tmp$V2

eqtl_Frontal_cortex_chr1$index<-paste0(eqtl_Frontal_cortex_chr1$chr,"-",eqtl_Frontal_cortex_chr1$loci)
Frontal_cortex_chr1$index<-paste0(Frontal_cortex_chr1$chr,"-",Frontal_cortex_chr1$pig_pos)
same_Frontal_cortex_chr1<-intersect(Frontal_cortex_chr1$pig_pos,eqtl_Frontal_cortex_chr1$loci)


tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex_chr2$variant_id,split="_"))))
eqtl_Frontal_cortex_chr2$chr<-tmp$V1
eqtl_Frontal_cortex_chr2$loci<-tmp$V2

eqtl_Frontal_cortex_chr2$index<-paste0(eqtl_Frontal_cortex_chr2$chr,"-",eqtl_Frontal_cortex_chr2$loci)
Frontal_cortex_chr2$index<-paste0(Frontal_cortex_chr2$chr,"-",Frontal_cortex_chr2$pig_pos)
same_Frontal_cortex_chr2<-intersect(Frontal_cortex_chr2$pig_pos,eqtl_Frontal_cortex_chr2$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex_chr3$variant_id,split="_"))))
eqtl_Frontal_cortex_chr3$chr<-tmp$V1
eqtl_Frontal_cortex_chr3$loci<-tmp$V2

eqtl_Frontal_cortex_chr3$index<-paste0(eqtl_Frontal_cortex_chr3$chr,"-",eqtl_Frontal_cortex_chr3$loci)
Frontal_cortex_chr3$index<-paste0(Frontal_cortex_chr3$chr,"-",Frontal_cortex_chr3$pig_pos)
same_Frontal_cortex_chr3<-intersect(Frontal_cortex_chr3$pig_pos,eqtl_Frontal_cortex_chr3$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex_chr4$variant_id,split="_"))))
eqtl_Frontal_cortex_chr4$chr<-tmp$V1
eqtl_Frontal_cortex_chr4$loci<-tmp$V2

eqtl_Frontal_cortex_chr4$index<-paste0(eqtl_Frontal_cortex_chr4$chr,"-",eqtl_Frontal_cortex_chr4$loci)
Frontal_cortex_chr4$index<-paste0(Frontal_cortex_chr4$chr,"-",Frontal_cortex_chr4$pig_pos)
same_Frontal_cortex_chr4<-intersect(Frontal_cortex_chr4$pig_pos,eqtl_Frontal_cortex_chr4$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex_chr5$variant_id,split="_"))))
eqtl_Frontal_cortex_chr5$loci<-tmp$V2
same_Frontal_cortex_chr5<-intersect(Frontal_cortex_chr5$pig_pos,eqtl_Frontal_cortex_chr5$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex_chr6$variant_id,split="_"))))
eqtl_Frontal_cortex_chr6$loci<-tmp$V2
same_Frontal_cortex_chr6<-intersect(Frontal_cortex_chr6$pig_pos,eqtl_Frontal_cortex_chr6$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex_chr7$variant_id,split="_"))))
eqtl_Frontal_cortex_chr7$loci<-tmp$V2
same_Frontal_cortex_chr7<-intersect(Frontal_cortex_chr7$pig_pos,eqtl_Frontal_cortex_chr7$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex_chr8$variant_id,split="_"))))
eqtl_Frontal_cortex_chr8$loci<-tmp$V2
same_Frontal_cortex_chr8<-intersect(Frontal_cortex_chr8$pig_pos,eqtl_Frontal_cortex_chr8$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex_chr9$variant_id,split="_"))))
eqtl_Frontal_cortex_chr9$loci<-tmp$V2
same_Frontal_cortex_chr9<-intersect(Frontal_cortex_chr9$pig_pos,eqtl_Frontal_cortex_chr9$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex_chr10$variant_id,split="_"))))
eqtl_Frontal_cortex_chr10$loci<-tmp$V2
same_Frontal_cortex_chr10<-intersect(Frontal_cortex_chr10$pig_pos,eqtl_Frontal_cortex_chr10$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex_chr11$variant_id,split="_"))))
eqtl_Frontal_cortex_chr11$loci<-tmp$V2
same_Frontal_cortex_chr11<-intersect(Frontal_cortex_chr11$pig_pos,eqtl_Frontal_cortex_chr11$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex_chr12$variant_id,split="_"))))
eqtl_Frontal_cortex_chr12$loci<-tmp$V2
same_Frontal_cortex_chr12<-intersect(Frontal_cortex_chr12$pig_pos,eqtl_Frontal_cortex_chr12$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex_chr13$variant_id,split="_"))))
eqtl_Frontal_cortex_chr13$loci<-tmp$V2
same_Frontal_cortex_chr13<-intersect(Frontal_cortex_chr13$pig_pos,eqtl_Frontal_cortex_chr13$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex_chr14$variant_id,split="_"))))
eqtl_Frontal_cortex_chr14$loci<-tmp$V2
same_Frontal_cortex_chr14<-intersect(Frontal_cortex_chr14$pig_pos,eqtl_Frontal_cortex_chr14$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex_chr15$variant_id,split="_"))))
eqtl_Frontal_cortex_chr15$loci<-tmp$V2
same_Frontal_cortex_chr15<-intersect(Frontal_cortex_chr15$pig_pos,eqtl_Frontal_cortex_chr15$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex_chr16$variant_id,split="_"))))
eqtl_Frontal_cortex_chr16$loci<-tmp$V2
same_Frontal_cortex_chr16<-intersect(Frontal_cortex_chr16$pig_pos,eqtl_Frontal_cortex_chr16$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex_chr17$variant_id,split="_"))))
eqtl_Frontal_cortex_chr17$loci<-tmp$V2
same_Frontal_cortex_chr17<-intersect(Frontal_cortex_chr17$pig_pos,eqtl_Frontal_cortex_chr17$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Frontal_cortex_chr18$variant_id,split="_"))))
eqtl_Frontal_cortex_chr18$loci<-tmp$V2
same_Frontal_cortex_chr18<-intersect(Frontal_cortex_chr18$pig_pos,eqtl_Frontal_cortex_chr18$loci)

chr1_Frontal_cortex_genes<-NULL
chr1_Frontal_cortex_overloci<-NULL
if(length(same_Frontal_cortex_chr1!=0)){
  for(i in 1:length(same_Frontal_cortex_chr1)){
    a<-eqtl_Frontal_cortex_chr1$phenotype_id[grep(same_Frontal_cortex_chr1[i],eqtl_Frontal_cortex_chr1$loci)]
    b<-eqtl_Frontal_cortex_chr1[grep(same_Frontal_cortex_chr1[i],eqtl_Frontal_cortex_chr1$loci),]
    chr1_Frontal_cortex_genes<-c(chr1_Frontal_cortex_genes,a)
    chr1_Frontal_cortex_overloci<-rbind(chr1_Frontal_cortex_overloci,b)
  }
}

chr2_Frontal_cortex_genes<-NULL
chr2_Frontal_cortex_overloci<-NULL
if(length(same_Frontal_cortex_chr2!=0)){
  for(i in 1:length(same_Frontal_cortex_chr2)){
    a<-eqtl_Frontal_cortex_chr2$phenotype_id[grep(same_Frontal_cortex_chr2[i],eqtl_Frontal_cortex_chr2$loci)]
    b<-eqtl_Frontal_cortex_chr2[grep(same_Frontal_cortex_chr2[i],eqtl_Frontal_cortex_chr2$loci),]
    chr2_Frontal_cortex_genes<-c(chr2_Frontal_cortex_genes,a)
    chr2_Frontal_cortex_overloci<-rbind(chr2_Frontal_cortex_overloci,b)
  }
}

chr3_Frontal_cortex_genes<-NULL
chr3_Frontal_cortex_overloci<-NULL
if(length(same_Frontal_cortex_chr3!=0)){
  for(i in 1:length(same_Frontal_cortex_chr3)){
    a<-eqtl_Frontal_cortex_chr3$phenotype_id[grep(same_Frontal_cortex_chr3[i],eqtl_Frontal_cortex_chr3$loci)]
    b<-eqtl_Frontal_cortex_chr3[grep(same_Frontal_cortex_chr3[i],eqtl_Frontal_cortex_chr3$loci),]
    chr3_Frontal_cortex_genes<-c(chr3_Frontal_cortex_genes,a)
    chr3_Frontal_cortex_overloci<-rbind(chr3_Frontal_cortex_overloci,b)
  }
}

chr4_Frontal_cortex_genes<-NULL
chr4_Frontal_cortex_overloci<-NULL
if(length(same_Frontal_cortex_chr4!=0)){
  for(i in 1:length(same_Frontal_cortex_chr4)){
    a<-eqtl_Frontal_cortex_chr4$phenotype_id[grep(same_Frontal_cortex_chr4[i],eqtl_Frontal_cortex_chr4$loci)]
    b<-eqtl_Frontal_cortex_chr4[grep(same_Frontal_cortex_chr4[i],eqtl_Frontal_cortex_chr4$loci),]
    chr4_Frontal_cortex_genes<-c(chr4_Frontal_cortex_genes,a)
    chr4_Frontal_cortex_overloci<-rbind(chr4_Frontal_cortex_overloci,b)
  }
}

chr5_Frontal_cortex_genes<-NULL
chr5_Frontal_cortex_overloci<-NULL
if(length(same_Frontal_cortex_chr5!=0)){
  for(i in 1:length(same_Frontal_cortex_chr5)){
    a<-eqtl_Frontal_cortex_chr5$phenotype_id[grep(same_Frontal_cortex_chr5[i],eqtl_Frontal_cortex_chr5$loci)]
    b<-eqtl_Frontal_cortex_chr5[grep(same_Frontal_cortex_chr5[i],eqtl_Frontal_cortex_chr5$loci),]
    chr5_Frontal_cortex_genes<-c(chr5_Frontal_cortex_genes,a)
    chr5_Frontal_cortex_overloci<-rbind(chr5_Frontal_cortex_overloci,b)
  }
}

chr6_Frontal_cortex_genes<-NULL
chr6_Frontal_cortex_overloci<-NULL
if(length(same_Frontal_cortex_chr6!=0)){
  for(i in 1:length(same_Frontal_cortex_chr6)){
    a<-eqtl_Frontal_cortex_chr6$phenotype_id[grep(same_Frontal_cortex_chr6[i],eqtl_Frontal_cortex_chr6$loci)]
    b<-eqtl_Frontal_cortex_chr6[grep(same_Frontal_cortex_chr6[i],eqtl_Frontal_cortex_chr6$loci),]
    chr6_Frontal_cortex_genes<-c(chr6_Frontal_cortex_genes,a)
    chr6_Frontal_cortex_overloci<-rbind(chr6_Frontal_cortex_overloci,b)
  }
}

chr7_Frontal_cortex_genes<-NULL
chr7_Frontal_cortex_overloci<-NULL
if(length(same_Frontal_cortex_chr7!=0)){
  for(i in 1:length(same_Frontal_cortex_chr7)){
    a<-eqtl_Frontal_cortex_chr7$phenotype_id[grep(same_Frontal_cortex_chr7[i],eqtl_Frontal_cortex_chr7$loci)]
    b<-eqtl_Frontal_cortex_chr7[grep(same_Frontal_cortex_chr7[i],eqtl_Frontal_cortex_chr7$loci),]
    chr7_Frontal_cortex_genes<-c(chr7_Frontal_cortex_genes,a)
    chr7_Frontal_cortex_overloci<-rbind(chr7_Frontal_cortex_overloci,b)
  }
}

chr8_Frontal_cortex_genes<-NULL
chr8_Frontal_cortex_overloci<-NULL
if(length(same_Frontal_cortex_chr8!=0)){
  for(i in 1:length(same_Frontal_cortex_chr8)){
    a<-eqtl_Frontal_cortex_chr8$phenotype_id[grep(same_Frontal_cortex_chr8[i],eqtl_Frontal_cortex_chr8$loci)]
    b<-eqtl_Frontal_cortex_chr8[grep(same_Frontal_cortex_chr8[i],eqtl_Frontal_cortex_chr8$loci),]
    chr8_Frontal_cortex_genes<-c(chr8_Frontal_cortex_genes,a)
    chr8_Frontal_cortex_overloci<-rbind(chr8_Frontal_cortex_overloci,b)
  }
}

if(length(same_Frontal_cortex_chr9!=0)){
  chr9_Frontal_cortex_genes<-NULL
  chr9_Frontal_cortex_overloci<-NULL
  for(i in 1:length(same_Frontal_cortex_chr9)){
    a<-eqtl_Frontal_cortex_chr9$phenotype_id[grep(same_Frontal_cortex_chr9[i],eqtl_Frontal_cortex_chr9$loci)]
    b<-eqtl_Frontal_cortex_chr9[grep(same_Frontal_cortex_chr9[i],eqtl_Frontal_cortex_chr9$loci),]
    chr9_Frontal_cortex_genes<-c(chr9_Frontal_cortex_genes,a)
    chr9_Frontal_cortex_overloci<-rbind(chr9_Frontal_cortex_overloci,b)
  }
}

chr10_Frontal_cortex_genes<-NULL
chr10_Frontal_cortex_overloci<-NULL
if(length(same_Frontal_cortex_chr10!=0)){
  for(i in 1:length(same_Frontal_cortex_chr10)){
    a<-eqtl_Frontal_cortex_chr10$phenotype_id[grep(same_Frontal_cortex_chr10[i],eqtl_Frontal_cortex_chr10$loci)]
    b<-eqtl_Frontal_cortex_chr10[grep(same_Frontal_cortex_chr10[i],eqtl_Frontal_cortex_chr10$loci),]
    chr10_Frontal_cortex_genes<-c(chr10_Frontal_cortex_genes,a)
    chr10_Frontal_cortex_overloci<-rbind(chr10_Frontal_cortex_overloci,b)
  }
}

chr11_Frontal_cortex_genes<-NULL
chr11_Frontal_cortex_overloci<-NULL
if(length(same_Frontal_cortex_chr11!=0)){
  for(i in 1:length(same_Frontal_cortex_chr11)){
    a<-eqtl_Frontal_cortex_chr11$phenotype_id[grep(same_Frontal_cortex_chr11[i],eqtl_Frontal_cortex_chr11$loci)]
    b<-eqtl_Frontal_cortex_chr11[grep(same_Frontal_cortex_chr11[i],eqtl_Frontal_cortex_chr11$loci),]
    chr11_Frontal_cortex_genes<-c(chr11_Frontal_cortex_genes,a)
    chr11_Frontal_cortex_overloci<-rbind(chr11_Frontal_cortex_overloci,b)
  }
}

chr12_Frontal_cortex_genes<-NULL
chr12_Frontal_cortex_overloci<-NULL
if(length(same_Frontal_cortex_chr12!=0)){
  for(i in 1:length(same_Frontal_cortex_chr12)){
    a<-eqtl_Frontal_cortex_chr12$phenotype_id[grep(same_Frontal_cortex_chr12[i],eqtl_Frontal_cortex_chr12$loci)]
    b<-eqtl_Frontal_cortex_chr12[grep(same_Frontal_cortex_chr12[i],eqtl_Frontal_cortex_chr12$loci),]
    chr12_Frontal_cortex_genes<-c(chr12_Frontal_cortex_genes,a)
    chr12_Frontal_cortex_overloci<-rbind(chr12_Frontal_cortex_overloci,b)
  }
}

chr13_Frontal_cortex_genes<-NULL
chr13_Frontal_cortex_overloci<-NULL
if(length(same_Frontal_cortex_chr13!=0)){
  for(i in 1:length(same_Frontal_cortex_chr13)){
    a<-eqtl_Frontal_cortex_chr13$phenotype_id[grep(same_Frontal_cortex_chr13[i],eqtl_Frontal_cortex_chr13$loci)]
    b<-eqtl_Frontal_cortex_chr13[grep(same_Frontal_cortex_chr13[i],eqtl_Frontal_cortex_chr13$loci),]
    chr13_Frontal_cortex_genes<-c(chr13_Frontal_cortex_genes,a)
    chr13_Frontal_cortex_overloci<-rbind(chr13_Frontal_cortex_overloci,b)
  }
}

chr14_Frontal_cortex_genes<-NULL
chr14_Frontal_cortex_overloci<-NULL
if(length(same_Frontal_cortex_chr14!=0)){
  for(i in 1:length(same_Frontal_cortex_chr14)){
    a<-eqtl_Frontal_cortex_chr14$phenotype_id[grep(same_Frontal_cortex_chr14[i],eqtl_Frontal_cortex_chr14$loci)]
    b<-eqtl_Frontal_cortex_chr14[grep(same_Frontal_cortex_chr14[i],eqtl_Frontal_cortex_chr14$loci),]
    chr14_Frontal_cortex_genes<-c(chr14_Frontal_cortex_genes,a)
    chr14_Frontal_cortex_overloci<-rbind(chr14_Frontal_cortex_overloci,b)
  }
}

chr15_Frontal_cortex_genes<-NULL
chr15_Frontal_cortex_overloci<-NULL
if(length(same_Frontal_cortex_chr15!=0)){
  for(i in 1:length(same_Frontal_cortex_chr15)){
    a<-eqtl_Frontal_cortex_chr15$phenotype_id[grep(same_Frontal_cortex_chr15[i],eqtl_Frontal_cortex_chr15$loci)]
    b<-eqtl_Frontal_cortex_chr15[grep(same_Frontal_cortex_chr15[i],eqtl_Frontal_cortex_chr15$loci),]
    chr15_Frontal_cortex_genes<-c(chr15_Frontal_cortex_genes,a)
    chr15_Frontal_cortex_overloci<-rbind(chr15_Frontal_cortex_overloci,b)
  }
}

chr16_Frontal_cortex_genes<-NULL
chr16_Frontal_cortex_overloci<-NULL
if(length(same_Frontal_cortex_chr16!=0)){
  for(i in 1:length(same_Frontal_cortex_chr16)){
    a<-eqtl_Frontal_cortex_chr16$phenotype_id[grep(same_Frontal_cortex_chr16[i],eqtl_Frontal_cortex_chr16$loci)]
    b<-eqtl_Frontal_cortex_chr16[grep(same_Frontal_cortex_chr16[i],eqtl_Frontal_cortex_chr16$loci),]
    chr16_Frontal_cortex_genes<-c(chr16_Frontal_cortex_genes,a)
    chr16_Frontal_cortex_overloci<-rbind(chr16_Frontal_cortex_overloci,b)
  }
}

chr17_Frontal_cortex_genes<-NULL
chr17_Frontal_cortex_overloci<-NULL
if(length(same_Frontal_cortex_chr17!=0)){
  for(i in 1:length(same_Frontal_cortex_chr17)){
    a<-eqtl_Frontal_cortex_chr17$phenotype_id[grep(same_Frontal_cortex_chr17[i],eqtl_Frontal_cortex_chr17$loci)]
    b<-eqtl_Frontal_cortex_chr17[grep(same_Frontal_cortex_chr17[i],eqtl_Frontal_cortex_chr17$loci),]
    chr17_Frontal_cortex_genes<-c(chr17_Frontal_cortex_genes,a)
    chr17_Frontal_cortex_overloci<-rbind(chr17_Frontal_cortex_overloci,b)
  }
}

chr18_Frontal_cortex_genes<-NULL
chr18_Frontal_cortex_overloci<-NULL
if(length(same_Frontal_cortex_chr18!=0)){
  for(i in 1:length(same_Frontal_cortex_chr18)){
    a<-eqtl_Frontal_cortex_chr18$phenotype_id[grep(same_Frontal_cortex_chr18[i],eqtl_Frontal_cortex_chr18$loci)]
    b<-eqtl_Frontal_cortex_chr18[grep(same_Frontal_cortex_chr18[i],eqtl_Frontal_cortex_chr18$loci),]
    chr18_Frontal_cortex_genes<-c(chr18_Frontal_cortex_genes,a)
    chr18_Frontal_cortex_overloci<-rbind(chr18_Frontal_cortex_overloci,b)
  }
}
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
one2one_pig<-annotation$`Pig gene stable ID`[which(annotation$`Pig homology type`=="ortholog_one2one")]

chr1_Frontal_cortex_overloci_one2one<-chr1_Frontal_cortex_overloci[match(intersect(chr1_Frontal_cortex_genes,one2one_pig),chr1_Frontal_cortex_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Frontal_cortex_chr1$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Frontal_cortex_chr1$hum_chr<-tmp$V1
Frontal_cortex_chr1$hum_loci<-tmp$V3
FM_Frontal_cortex_hum$index<-paste0(FM_Frontal_cortex_hum$chr,"_",FM_Frontal_cortex_hum$variant_pos)
Frontal_cortex_chr1$index<-paste0("chr",Frontal_cortex_chr1$hum_chr,"_",Frontal_cortex_chr1$hum_loci)
chr1_Frontal_cortex_hum<-NULL
for(i in 1:length(same_Frontal_cortex_chr1)){
  a<-Frontal_cortex_chr1[match(same_Frontal_cortex_chr1[i],Frontal_cortex_chr1$pig_pos),]
  b<-FM_Frontal_cortex_hum[match(intersect(a$index,FM_Frontal_cortex_hum$index),FM_Frontal_cortex_hum$index),]
  chr1_Frontal_cortex_hum<-rbind(chr1_Frontal_cortex_hum,b)
}
chr1_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr1_Frontal_cortex_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr1_Frontal_cortex_hum_one2one<-chr1_Frontal_cortex_hum[match(intersect(chr1_pig2hum_one2one,chr1_Frontal_cortex_hum$gene_id),chr1_Frontal_cortex_hum$gene_id),]

chr2_Frontal_cortex_overloci_one2one<-chr2_Frontal_cortex_overloci[match(intersect(chr2_Frontal_cortex_genes,one2one_pig),chr2_Frontal_cortex_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Frontal_cortex_chr2$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Frontal_cortex_chr2$hum_chr<-tmp$V1
Frontal_cortex_chr2$hum_loci<-tmp$V3
FM_Frontal_cortex_hum$index<-paste0(FM_Frontal_cortex_hum$chr,"_",FM_Frontal_cortex_hum$variant_pos)
Frontal_cortex_chr2$index<-paste0("chr",Frontal_cortex_chr2$hum_chr,"_",Frontal_cortex_chr2$hum_loci)
chr2_Frontal_cortex_hum<-NULL
for(i in 1:length(same_Frontal_cortex_chr2)){
  a<-Frontal_cortex_chr2[match(same_Frontal_cortex_chr2[i],Frontal_cortex_chr2$pig_pos),]
  b<-FM_Frontal_cortex_hum[match(intersect(a$index,FM_Frontal_cortex_hum$index),FM_Frontal_cortex_hum$index),]
  chr2_Frontal_cortex_hum<-rbind(chr2_Frontal_cortex_hum,b)
}
chr2_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr2_Frontal_cortex_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr2_Frontal_cortex_hum_one2one<-chr2_Frontal_cortex_hum[match(intersect(chr2_pig2hum_one2one,chr2_Frontal_cortex_hum$gene_id),chr2_Frontal_cortex_hum$gene_id),]

chr3_Frontal_cortex_overloci_one2one<-chr3_Frontal_cortex_overloci[match(intersect(chr3_Frontal_cortex_genes,one2one_pig),chr3_Frontal_cortex_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Frontal_cortex_chr3$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Frontal_cortex_chr3$hum_chr<-tmp$V1
Frontal_cortex_chr3$hum_loci<-tmp$V3
FM_Frontal_cortex_hum$index<-paste0(FM_Frontal_cortex_hum$chr,"_",FM_Frontal_cortex_hum$variant_pos)
Frontal_cortex_chr3$index<-paste0("chr",Frontal_cortex_chr3$hum_chr,"_",Frontal_cortex_chr3$hum_loci)
chr3_Frontal_cortex_hum<-NULL
for(i in 1:length(same_Frontal_cortex_chr3)){
  a<-Frontal_cortex_chr3[match(same_Frontal_cortex_chr3[i],Frontal_cortex_chr3$pig_pos),]
  b<-FM_Frontal_cortex_hum[match(intersect(a$index,FM_Frontal_cortex_hum$index),FM_Frontal_cortex_hum$index),]
  chr3_Frontal_cortex_hum<-rbind(chr3_Frontal_cortex_hum,b)
}
chr3_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr3_Frontal_cortex_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr3_Frontal_cortex_hum_one2one<-chr3_Frontal_cortex_hum[match(intersect(chr3_pig2hum_one2one,chr3_Frontal_cortex_hum$gene_id),chr3_Frontal_cortex_hum$gene_id),]

chr4_Frontal_cortex_overloci_one2one<-chr4_Frontal_cortex_overloci[match(intersect(chr4_Frontal_cortex_genes,one2one_pig),chr4_Frontal_cortex_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Frontal_cortex_chr4$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Frontal_cortex_chr4$hum_chr<-tmp$V1
Frontal_cortex_chr4$hum_loci<-tmp$V3
FM_Frontal_cortex_hum$index<-paste0(FM_Frontal_cortex_hum$chr,"_",FM_Frontal_cortex_hum$variant_pos)
Frontal_cortex_chr4$index<-paste0("chr",Frontal_cortex_chr4$hum_chr,"_",Frontal_cortex_chr4$hum_loci)
chr4_Frontal_cortex_hum<-NULL
for(i in 1:length(same_Frontal_cortex_chr4)){
  a<-Frontal_cortex_chr4[match(same_Frontal_cortex_chr4[i],Frontal_cortex_chr4$pig_pos),]
  b<-FM_Frontal_cortex_hum[match(intersect(a$index,FM_Frontal_cortex_hum$index),FM_Frontal_cortex_hum$index),]
  chr4_Frontal_cortex_hum<-rbind(chr4_Frontal_cortex_hum,b)
}
chr4_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr4_Frontal_cortex_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr4_Frontal_cortex_hum_one2one<-chr4_Frontal_cortex_hum[match(intersect(chr4_pig2hum_one2one,chr4_Frontal_cortex_hum$gene_id),chr4_Frontal_cortex_hum$gene_id),]

chr5_Frontal_cortex_overloci_one2one<-chr5_Frontal_cortex_overloci[match(intersect(chr5_Frontal_cortex_genes,one2one_pig),chr5_Frontal_cortex_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Frontal_cortex_chr5$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Frontal_cortex_chr5$hum_chr<-tmp$V1
Frontal_cortex_chr5$hum_loci<-tmp$V3
FM_Frontal_cortex_hum$index<-paste0(FM_Frontal_cortex_hum$chr,"_",FM_Frontal_cortex_hum$variant_pos)
Frontal_cortex_chr5$index<-paste0("chr",Frontal_cortex_chr5$hum_chr,"_",Frontal_cortex_chr5$hum_loci)
chr5_Frontal_cortex_hum<-NULL
for(i in 1:length(same_Frontal_cortex_chr5)){
  a<-Frontal_cortex_chr5[match(same_Frontal_cortex_chr5[i],Frontal_cortex_chr5$pig_pos),]
  b<-FM_Frontal_cortex_hum[match(intersect(a$index,FM_Frontal_cortex_hum$index),FM_Frontal_cortex_hum$index),]
  chr5_Frontal_cortex_hum<-rbind(chr5_Frontal_cortex_hum,b)
}
chr5_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr5_Frontal_cortex_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr5_Frontal_cortex_hum_one2one<-chr5_Frontal_cortex_hum[match(intersect(chr5_pig2hum_one2one,chr5_Frontal_cortex_hum$gene_id),chr5_Frontal_cortex_hum$gene_id),]

chr6_Frontal_cortex_overloci_one2one<-chr6_Frontal_cortex_overloci[match(intersect(chr6_Frontal_cortex_genes,one2one_pig),chr6_Frontal_cortex_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Frontal_cortex_chr6$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Frontal_cortex_chr6$hum_chr<-tmp$V1
Frontal_cortex_chr6$hum_loci<-tmp$V3
FM_Frontal_cortex_hum$index<-paste0(FM_Frontal_cortex_hum$chr,"_",FM_Frontal_cortex_hum$variant_pos)
Frontal_cortex_chr6$index<-paste0("chr",Frontal_cortex_chr6$hum_chr,"_",Frontal_cortex_chr6$hum_loci)
chr6_Frontal_cortex_hum<-NULL
for(i in 1:length(same_Frontal_cortex_chr6)){
  a<-Frontal_cortex_chr6[match(same_Frontal_cortex_chr6[i],Frontal_cortex_chr6$pig_pos),]
  b<-FM_Frontal_cortex_hum[match(intersect(a$index,FM_Frontal_cortex_hum$index),FM_Frontal_cortex_hum$index),]
  chr6_Frontal_cortex_hum<-rbind(chr6_Frontal_cortex_hum,b)
}
chr6_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr6_Frontal_cortex_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr6_Frontal_cortex_hum_one2one<-chr6_Frontal_cortex_hum[match(intersect(chr6_pig2hum_one2one,chr6_Frontal_cortex_hum$gene_id),chr6_Frontal_cortex_hum$gene_id),]

chr7_Frontal_cortex_overloci_one2one<-chr7_Frontal_cortex_overloci[match(intersect(chr7_Frontal_cortex_genes,one2one_pig),chr7_Frontal_cortex_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Frontal_cortex_chr7$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Frontal_cortex_chr7$hum_chr<-tmp$V1
Frontal_cortex_chr7$hum_loci<-tmp$V3
FM_Frontal_cortex_hum$index<-paste0(FM_Frontal_cortex_hum$chr,"_",FM_Frontal_cortex_hum$variant_pos)
Frontal_cortex_chr7$index<-paste0("chr",Frontal_cortex_chr7$hum_chr,"_",Frontal_cortex_chr7$hum_loci)
chr7_Frontal_cortex_hum<-NULL
for(i in 1:length(same_Frontal_cortex_chr7)){
  a<-Frontal_cortex_chr7[match(same_Frontal_cortex_chr7[i],Frontal_cortex_chr7$pig_pos),]
  b<-FM_Frontal_cortex_hum[match(intersect(a$index,FM_Frontal_cortex_hum$index),FM_Frontal_cortex_hum$index),]
  chr7_Frontal_cortex_hum<-rbind(chr7_Frontal_cortex_hum,b)
}
chr7_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr7_Frontal_cortex_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr7_Frontal_cortex_hum_one2one<-chr7_Frontal_cortex_hum[match(intersect(chr7_pig2hum_one2one,chr7_Frontal_cortex_hum$gene_id),chr7_Frontal_cortex_hum$gene_id),]

chr8_Frontal_cortex_overloci_one2one<-chr8_Frontal_cortex_overloci[match(intersect(chr8_Frontal_cortex_genes,one2one_pig),chr8_Frontal_cortex_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Frontal_cortex_chr8$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Frontal_cortex_chr8$hum_chr<-tmp$V1
Frontal_cortex_chr8$hum_loci<-tmp$V3
FM_Frontal_cortex_hum$index<-paste0(FM_Frontal_cortex_hum$chr,"_",FM_Frontal_cortex_hum$variant_pos)
Frontal_cortex_chr8$index<-paste0("chr",Frontal_cortex_chr8$hum_chr,"_",Frontal_cortex_chr8$hum_loci)
chr8_Frontal_cortex_hum<-NULL
for(i in 1:length(same_Frontal_cortex_chr8)){
  a<-Frontal_cortex_chr8[match(same_Frontal_cortex_chr8[i],Frontal_cortex_chr8$pig_pos),]
  b<-FM_Frontal_cortex_hum[match(intersect(a$index,FM_Frontal_cortex_hum$index),FM_Frontal_cortex_hum$index),]
  chr8_Frontal_cortex_hum<-rbind(chr8_Frontal_cortex_hum,b)
}
chr8_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr8_Frontal_cortex_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr8_Frontal_cortex_hum_one2one<-chr8_Frontal_cortex_hum[match(intersect(chr8_pig2hum_one2one,chr8_Frontal_cortex_hum$gene_id),chr8_Frontal_cortex_hum$gene_id),]

chr9_Frontal_cortex_overloci_one2one<-chr9_Frontal_cortex_overloci[match(intersect(chr9_Frontal_cortex_genes,one2one_pig),chr9_Frontal_cortex_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Frontal_cortex_chr9$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Frontal_cortex_chr9$hum_chr<-tmp$V1
Frontal_cortex_chr9$hum_loci<-tmp$V3
FM_Frontal_cortex_hum$index<-paste0(FM_Frontal_cortex_hum$chr,"_",FM_Frontal_cortex_hum$variant_pos)
Frontal_cortex_chr9$index<-paste0("chr",Frontal_cortex_chr9$hum_chr,"_",Frontal_cortex_chr9$hum_loci)
chr9_Frontal_cortex_hum<-NULL
for(i in 1:length(same_Frontal_cortex_chr9)){
  a<-Frontal_cortex_chr9[match(same_Frontal_cortex_chr9[i],Frontal_cortex_chr9$pig_pos),]
  b<-FM_Frontal_cortex_hum[match(intersect(a$index,FM_Frontal_cortex_hum$index),FM_Frontal_cortex_hum$index),]
  chr9_Frontal_cortex_hum<-rbind(chr9_Frontal_cortex_hum,b)
}
chr9_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr9_Frontal_cortex_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr9_Frontal_cortex_hum_one2one<-chr9_Frontal_cortex_hum[match(intersect(chr9_pig2hum_one2one,chr9_Frontal_cortex_hum$gene_id),chr9_Frontal_cortex_hum$gene_id),]

chr10_Frontal_cortex_overloci_one2one<-chr10_Frontal_cortex_overloci[match(intersect(chr10_Frontal_cortex_genes,one2one_pig),chr10_Frontal_cortex_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Frontal_cortex_chr10$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Frontal_cortex_chr10$hum_chr<-tmp$V1
Frontal_cortex_chr10$hum_loci<-tmp$V3
FM_Frontal_cortex_hum$index<-paste0(FM_Frontal_cortex_hum$chr,"_",FM_Frontal_cortex_hum$variant_pos)
Frontal_cortex_chr10$index<-paste0("chr",Frontal_cortex_chr10$hum_chr,"_",Frontal_cortex_chr10$hum_loci)
chr10_Frontal_cortex_hum<-NULL
for(i in 1:length(same_Frontal_cortex_chr10)){
  a<-Frontal_cortex_chr10[match(same_Frontal_cortex_chr10[i],Frontal_cortex_chr10$pig_pos),]
  b<-FM_Frontal_cortex_hum[match(intersect(a$index,FM_Frontal_cortex_hum$index),FM_Frontal_cortex_hum$index),]
  chr10_Frontal_cortex_hum<-rbind(chr10_Frontal_cortex_hum,b)
}
chr10_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr10_Frontal_cortex_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr10_Frontal_cortex_hum_one2one<-chr10_Frontal_cortex_hum[match(intersect(chr10_pig2hum_one2one,chr10_Frontal_cortex_hum$gene_id),chr10_Frontal_cortex_hum$gene_id),]

chr11_Frontal_cortex_overloci_one2one<-chr11_Frontal_cortex_overloci[match(intersect(chr11_Frontal_cortex_genes,one2one_pig),chr11_Frontal_cortex_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Frontal_cortex_chr11$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Frontal_cortex_chr11$hum_chr<-tmp$V1
Frontal_cortex_chr11$hum_loci<-tmp$V3
FM_Frontal_cortex_hum$index<-paste0(FM_Frontal_cortex_hum$chr,"_",FM_Frontal_cortex_hum$variant_pos)
Frontal_cortex_chr11$index<-paste0("chr",Frontal_cortex_chr11$hum_chr,"_",Frontal_cortex_chr11$hum_loci)
chr11_Frontal_cortex_hum<-NULL
for(i in 1:length(same_Frontal_cortex_chr11)){
  a<-Frontal_cortex_chr11[match(same_Frontal_cortex_chr11[i],Frontal_cortex_chr11$pig_pos),]
  b<-FM_Frontal_cortex_hum[match(intersect(a$index,FM_Frontal_cortex_hum$index),FM_Frontal_cortex_hum$index),]
  chr11_Frontal_cortex_hum<-rbind(chr11_Frontal_cortex_hum,b)
}
chr11_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr11_Frontal_cortex_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr11_Frontal_cortex_hum_one2one<-chr11_Frontal_cortex_hum[match(intersect(chr11_pig2hum_one2one,chr11_Frontal_cortex_hum$gene_id),chr11_Frontal_cortex_hum$gene_id),]

chr12_Frontal_cortex_overloci_one2one<-chr12_Frontal_cortex_overloci[match(intersect(chr12_Frontal_cortex_genes,one2one_pig),chr12_Frontal_cortex_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Frontal_cortex_chr12$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Frontal_cortex_chr12$hum_chr<-tmp$V1
Frontal_cortex_chr12$hum_loci<-tmp$V3
FM_Frontal_cortex_hum$index<-paste0(FM_Frontal_cortex_hum$chr,"_",FM_Frontal_cortex_hum$variant_pos)
Frontal_cortex_chr12$index<-paste0("chr",Frontal_cortex_chr12$hum_chr,"_",Frontal_cortex_chr12$hum_loci)
chr12_Frontal_cortex_hum<-NULL
for(i in 1:length(same_Frontal_cortex_chr12)){
  a<-Frontal_cortex_chr12[match(same_Frontal_cortex_chr12[i],Frontal_cortex_chr12$pig_pos),]
  b<-FM_Frontal_cortex_hum[match(intersect(a$index,FM_Frontal_cortex_hum$index),FM_Frontal_cortex_hum$index),]
  chr12_Frontal_cortex_hum<-rbind(chr12_Frontal_cortex_hum,b)
}
chr12_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr12_Frontal_cortex_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr12_Frontal_cortex_hum_one2one<-chr12_Frontal_cortex_hum[match(intersect(chr12_pig2hum_one2one,chr12_Frontal_cortex_hum$gene_id),chr12_Frontal_cortex_hum$gene_id),]

chr13_Frontal_cortex_overloci_one2one<-chr13_Frontal_cortex_overloci[match(intersect(chr13_Frontal_cortex_genes,one2one_pig),chr13_Frontal_cortex_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Frontal_cortex_chr13$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Frontal_cortex_chr13$hum_chr<-tmp$V1
Frontal_cortex_chr13$hum_loci<-tmp$V3
FM_Frontal_cortex_hum$index<-paste0(FM_Frontal_cortex_hum$chr,"_",FM_Frontal_cortex_hum$variant_pos)
Frontal_cortex_chr13$index<-paste0("chr",Frontal_cortex_chr13$hum_chr,"_",Frontal_cortex_chr13$hum_loci)
chr13_Frontal_cortex_hum<-NULL
for(i in 1:length(same_Frontal_cortex_chr13)){
  a<-Frontal_cortex_chr13[match(same_Frontal_cortex_chr13[i],Frontal_cortex_chr13$pig_pos),]
  b<-FM_Frontal_cortex_hum[match(intersect(a$index,FM_Frontal_cortex_hum$index),FM_Frontal_cortex_hum$index),]
  chr13_Frontal_cortex_hum<-rbind(chr13_Frontal_cortex_hum,b)
}
chr13_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr13_Frontal_cortex_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr13_Frontal_cortex_hum_one2one<-chr13_Frontal_cortex_hum[match(intersect(chr13_pig2hum_one2one,chr13_Frontal_cortex_hum$gene_id),chr13_Frontal_cortex_hum$gene_id),]

chr14_Frontal_cortex_overloci_one2one<-chr14_Frontal_cortex_overloci[match(intersect(chr14_Frontal_cortex_genes,one2one_pig),chr14_Frontal_cortex_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Frontal_cortex_chr14$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Frontal_cortex_chr14$hum_chr<-tmp$V1
Frontal_cortex_chr14$hum_loci<-tmp$V3
FM_Frontal_cortex_hum$index<-paste0(FM_Frontal_cortex_hum$chr,"_",FM_Frontal_cortex_hum$variant_pos)
Frontal_cortex_chr14$index<-paste0("chr",Frontal_cortex_chr14$hum_chr,"_",Frontal_cortex_chr14$hum_loci)
chr14_Frontal_cortex_hum<-NULL
for(i in 1:length(same_Frontal_cortex_chr14)){
  a<-Frontal_cortex_chr14[match(same_Frontal_cortex_chr14[i],Frontal_cortex_chr14$pig_pos),]
  b<-FM_Frontal_cortex_hum[match(intersect(a$index,FM_Frontal_cortex_hum$index),FM_Frontal_cortex_hum$index),]
  chr14_Frontal_cortex_hum<-rbind(chr14_Frontal_cortex_hum,b)
}
chr14_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr14_Frontal_cortex_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr14_Frontal_cortex_hum_one2one<-chr14_Frontal_cortex_hum[match(intersect(chr14_pig2hum_one2one,chr14_Frontal_cortex_hum$gene_id),chr14_Frontal_cortex_hum$gene_id),]

chr15_Frontal_cortex_overloci_one2one<-chr15_Frontal_cortex_overloci[match(intersect(chr15_Frontal_cortex_genes,one2one_pig),chr15_Frontal_cortex_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Frontal_cortex_chr15$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Frontal_cortex_chr15$hum_chr<-tmp$V1
Frontal_cortex_chr15$hum_loci<-tmp$V3
FM_Frontal_cortex_hum$index<-paste0(FM_Frontal_cortex_hum$chr,"_",FM_Frontal_cortex_hum$variant_pos)
Frontal_cortex_chr15$index<-paste0("chr",Frontal_cortex_chr15$hum_chr,"_",Frontal_cortex_chr15$hum_loci)
chr15_Frontal_cortex_hum<-NULL
for(i in 1:length(same_Frontal_cortex_chr15)){
  a<-Frontal_cortex_chr15[match(same_Frontal_cortex_chr15[i],Frontal_cortex_chr15$pig_pos),]
  b<-FM_Frontal_cortex_hum[match(intersect(a$index,FM_Frontal_cortex_hum$index),FM_Frontal_cortex_hum$index),]
  chr15_Frontal_cortex_hum<-rbind(chr15_Frontal_cortex_hum,b)
}
chr15_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr15_Frontal_cortex_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr15_Frontal_cortex_hum_one2one<-chr15_Frontal_cortex_hum[match(intersect(chr15_pig2hum_one2one,chr15_Frontal_cortex_hum$gene_id),chr15_Frontal_cortex_hum$gene_id),]

chr16_Frontal_cortex_overloci_one2one<-chr16_Frontal_cortex_overloci[match(intersect(chr16_Frontal_cortex_genes,one2one_pig),chr16_Frontal_cortex_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Frontal_cortex_chr16$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Frontal_cortex_chr16$hum_chr<-tmp$V1
Frontal_cortex_chr16$hum_loci<-tmp$V3
FM_Frontal_cortex_hum$index<-paste0(FM_Frontal_cortex_hum$chr,"_",FM_Frontal_cortex_hum$variant_pos)
Frontal_cortex_chr16$index<-paste0("chr",Frontal_cortex_chr16$hum_chr,"_",Frontal_cortex_chr16$hum_loci)
chr16_Frontal_cortex_hum<-NULL
for(i in 1:length(same_Frontal_cortex_chr16)){
  a<-Frontal_cortex_chr16[match(same_Frontal_cortex_chr16[i],Frontal_cortex_chr16$pig_pos),]
  b<-FM_Frontal_cortex_hum[match(intersect(a$index,FM_Frontal_cortex_hum$index),FM_Frontal_cortex_hum$index),]
  chr16_Frontal_cortex_hum<-rbind(chr16_Frontal_cortex_hum,b)
}
chr16_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr16_Frontal_cortex_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr16_Frontal_cortex_hum_one2one<-chr16_Frontal_cortex_hum[match(intersect(chr16_pig2hum_one2one,chr16_Frontal_cortex_hum$gene_id),chr16_Frontal_cortex_hum$gene_id),]

chr17_Frontal_cortex_overloci_one2one<-chr17_Frontal_cortex_overloci[match(intersect(chr17_Frontal_cortex_genes,one2one_pig),chr17_Frontal_cortex_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Frontal_cortex_chr17$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Frontal_cortex_chr17$hum_chr<-tmp$V1
Frontal_cortex_chr17$hum_loci<-tmp$V3
FM_Frontal_cortex_hum$index<-paste0(FM_Frontal_cortex_hum$chr,"_",FM_Frontal_cortex_hum$variant_pos)
Frontal_cortex_chr17$index<-paste0("chr",Frontal_cortex_chr17$hum_chr,"_",Frontal_cortex_chr17$hum_loci)
chr17_Frontal_cortex_hum<-NULL
for(i in 1:length(same_Frontal_cortex_chr17)){
  a<-Frontal_cortex_chr17[match(same_Frontal_cortex_chr17[i],Frontal_cortex_chr17$pig_pos),]
  b<-FM_Frontal_cortex_hum[match(intersect(a$index,FM_Frontal_cortex_hum$index),FM_Frontal_cortex_hum$index),]
  chr17_Frontal_cortex_hum<-rbind(chr17_Frontal_cortex_hum,b)
}
chr17_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr17_Frontal_cortex_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr17_Frontal_cortex_hum_one2one<-chr17_Frontal_cortex_hum[match(intersect(chr17_pig2hum_one2one,chr17_Frontal_cortex_hum$gene_id),chr17_Frontal_cortex_hum$gene_id),]

chr18_Frontal_cortex_overloci_one2one<-chr18_Frontal_cortex_overloci[match(intersect(chr18_Frontal_cortex_genes,one2one_pig),chr18_Frontal_cortex_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Frontal_cortex_chr18$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Frontal_cortex_chr18$hum_chr<-tmp$V1
Frontal_cortex_chr18$hum_loci<-tmp$V3
FM_Frontal_cortex_hum$index<-paste0(FM_Frontal_cortex_hum$chr,"_",FM_Frontal_cortex_hum$variant_pos)
Frontal_cortex_chr18$index<-paste0("chr",Frontal_cortex_chr18$hum_chr,"_",Frontal_cortex_chr18$hum_loci)
chr18_Frontal_cortex_hum<-NULL
for(i in 1:length(same_Frontal_cortex_chr18)){
  a<-Frontal_cortex_chr18[match(same_Frontal_cortex_chr18[i],Frontal_cortex_chr18$pig_pos),]
  b<-FM_Frontal_cortex_hum[match(intersect(a$index,FM_Frontal_cortex_hum$index),FM_Frontal_cortex_hum$index),]
  chr18_Frontal_cortex_hum<-rbind(chr18_Frontal_cortex_hum,b)
}
chr18_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr18_Frontal_cortex_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr18_Frontal_cortex_hum_one2one<-chr18_Frontal_cortex_hum[match(intersect(chr18_pig2hum_one2one,chr18_Frontal_cortex_hum$gene_id),chr18_Frontal_cortex_hum$gene_id),]

Frontal_cortex_one2one_SNP_hum<-rbind(chr1_Frontal_cortex_hum_one2one,chr2_Frontal_cortex_hum_one2one,chr3_Frontal_cortex_hum_one2one,chr4_Frontal_cortex_hum_one2one,chr5_Frontal_cortex_hum_one2one,
                                      chr6_Frontal_cortex_hum_one2one,chr7_Frontal_cortex_hum_one2one,chr8_Frontal_cortex_hum_one2one,chr9_Frontal_cortex_hum_one2one,chr10_Frontal_cortex_hum_one2one,
                                      chr11_Frontal_cortex_hum_one2one,chr12_Frontal_cortex_hum_one2one,chr13_Frontal_cortex_hum_one2one,chr14_Frontal_cortex_hum_one2one,chr15_Frontal_cortex_hum_one2one,
                                      chr16_Frontal_cortex_hum_one2one,chr17_Frontal_cortex_hum_one2one,chr18_Frontal_cortex_hum_one2one)

chr1_Frontal_cortex_pig_one2one<-chr1_Frontal_cortex_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr1_Frontal_cortex_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr1_Frontal_cortex_overloci_one2one$phenotype_id),1:9]
chr2_Frontal_cortex_pig_one2one<-chr2_Frontal_cortex_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr2_Frontal_cortex_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr2_Frontal_cortex_overloci_one2one$phenotype_id),1:9]
chr3_Frontal_cortex_pig_one2one<-chr3_Frontal_cortex_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr3_Frontal_cortex_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr3_Frontal_cortex_overloci_one2one$phenotype_id),1:9]
chr4_Frontal_cortex_pig_one2one<-chr4_Frontal_cortex_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr4_Frontal_cortex_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr4_Frontal_cortex_overloci_one2one$phenotype_id),1:9]
chr5_Frontal_cortex_pig_one2one<-chr5_Frontal_cortex_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr5_Frontal_cortex_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr5_Frontal_cortex_overloci_one2one$phenotype_id),1:9]
chr6_Frontal_cortex_pig_one2one<-chr6_Frontal_cortex_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr6_Frontal_cortex_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr6_Frontal_cortex_overloci_one2one$phenotype_id),1:9]
chr7_Frontal_cortex_pig_one2one<-chr7_Frontal_cortex_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr7_Frontal_cortex_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr7_Frontal_cortex_overloci_one2one$phenotype_id),1:9]
chr8_Frontal_cortex_pig_one2one<-chr8_Frontal_cortex_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr8_Frontal_cortex_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr8_Frontal_cortex_overloci_one2one$phenotype_id),1:9]
chr9_Frontal_cortex_pig_one2one<-chr9_Frontal_cortex_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr9_Frontal_cortex_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr9_Frontal_cortex_overloci_one2one$phenotype_id),1:9]
chr10_Frontal_cortex_pig_one2one<-chr10_Frontal_cortex_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr10_Frontal_cortex_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr10_Frontal_cortex_overloci_one2one$phenotype_id),1:9]
chr11_Frontal_cortex_pig_one2one<-chr11_Frontal_cortex_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr11_Frontal_cortex_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr11_Frontal_cortex_overloci_one2one$phenotype_id),1:9]
chr12_Frontal_cortex_pig_one2one<-chr12_Frontal_cortex_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr12_Frontal_cortex_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr12_Frontal_cortex_overloci_one2one$phenotype_id),1:9]
chr13_Frontal_cortex_pig_one2one<-chr13_Frontal_cortex_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr13_Frontal_cortex_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr13_Frontal_cortex_overloci_one2one$phenotype_id),1:9]
chr14_Frontal_cortex_pig_one2one<-chr14_Frontal_cortex_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr14_Frontal_cortex_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr14_Frontal_cortex_overloci_one2one$phenotype_id),1:9]
chr15_Frontal_cortex_pig_one2one<-chr15_Frontal_cortex_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr15_Frontal_cortex_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr15_Frontal_cortex_overloci_one2one$phenotype_id),1:9]
chr16_Frontal_cortex_pig_one2one<-chr16_Frontal_cortex_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr16_Frontal_cortex_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr16_Frontal_cortex_overloci_one2one$phenotype_id),1:9]
chr17_Frontal_cortex_pig_one2one<-chr17_Frontal_cortex_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr17_Frontal_cortex_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr17_Frontal_cortex_overloci_one2one$phenotype_id),1:9]
chr18_Frontal_cortex_pig_one2one<-chr18_Frontal_cortex_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr18_Frontal_cortex_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr18_Frontal_cortex_overloci_one2one$phenotype_id),1:9]

Frontal_cortex_one2one_SNP_pig<-rbind(chr1_Frontal_cortex_pig_one2one,chr2_Frontal_cortex_pig_one2one,chr3_Frontal_cortex_pig_one2one,chr4_Frontal_cortex_pig_one2one,chr5_Frontal_cortex_pig_one2one,
                                      chr6_Frontal_cortex_pig_one2one,chr7_Frontal_cortex_pig_one2one,chr8_Frontal_cortex_pig_one2one,chr9_Frontal_cortex_pig_one2one,chr10_Frontal_cortex_pig_one2one,
                                      chr11_Frontal_cortex_pig_one2one,chr12_Frontal_cortex_pig_one2one,chr13_Frontal_cortex_pig_one2one,chr14_Frontal_cortex_pig_one2one,chr15_Frontal_cortex_pig_one2one,
                                      chr16_Frontal_cortex_pig_one2one,chr17_Frontal_cortex_pig_one2one,chr18_Frontal_cortex_pig_one2one)

Frontal_cortex_SNP_sum<-array(NA,dim=c(nrow(Frontal_cortex_one2one_SNP_hum),2))
colnames(Frontal_cortex_SNP_sum)<-c("Human","Pig")
Frontal_cortex_SNP_sum<-as.data.frame(Frontal_cortex_SNP_sum)
Frontal_cortex_SNP_sum$Human<-Frontal_cortex_one2one_SNP_hum$slope / Frontal_cortex_one2one_SNP_hum$slope_se
Frontal_cortex_SNP_sum$Pig<-Frontal_cortex_one2one_SNP_pig$slope / Frontal_cortex_one2one_SNP_pig$slope_se
cor<-cor(abs(Frontal_cortex_SNP_sum$Human),abs(Frontal_cortex_SNP_sum$Pig))
p_val<-t.test(abs(Frontal_cortex_SNP_sum$Human),abs(Frontal_cortex_SNP_sum$Pig))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/Frontal_cortex_SNP_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(Frontal_cortex_SNP_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The number of eGenes expressed in common tissues",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

save(Frontal_cortex_one2one_SNP_hum,Frontal_cortex_one2one_SNP_pig,file="/Users/baizhonghao/Downloads/human-pig GTEx/Frontal_cortex_SNP.Rdata")

#Heart_SNP_overlaploci#
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

eqtl_Heart_chr1<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Heart/Heart.cis_qtl_pairs.1.txt"))
eqtl_Heart_chr2<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Heart/Heart.cis_qtl_pairs.2.txt"))
eqtl_Heart_chr3<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Heart/Heart.cis_qtl_pairs.3.txt"))
eqtl_Heart_chr4<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Heart/Heart.cis_qtl_pairs.4.txt"))
eqtl_Heart_chr5<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Heart/Heart.cis_qtl_pairs.5.txt"))
eqtl_Heart_chr6<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Heart/Heart.cis_qtl_pairs.6.txt"))
eqtl_Heart_chr7<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Heart/Heart.cis_qtl_pairs.7.txt"))
eqtl_Heart_chr8<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Heart/Heart.cis_qtl_pairs.8.txt"))
eqtl_Heart_chr9<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Heart/Heart.cis_qtl_pairs.9.txt"))
eqtl_Heart_chr10<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Heart/Heart.cis_qtl_pairs.10.txt"))
eqtl_Heart_chr11<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Heart/Heart.cis_qtl_pairs.11.txt"))
eqtl_Heart_chr12<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Heart/Heart.cis_qtl_pairs.12.txt"))
eqtl_Heart_chr13<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Heart/Heart.cis_qtl_pairs.13.txt"))
eqtl_Heart_chr14<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Heart/Heart.cis_qtl_pairs.14.txt"))
eqtl_Heart_chr15<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Heart/Heart.cis_qtl_pairs.15.txt"))
eqtl_Heart_chr16<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Heart/Heart.cis_qtl_pairs.16.txt"))
eqtl_Heart_chr17<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Heart/Heart.cis_qtl_pairs.17.txt"))
eqtl_Heart_chr18<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Heart/Heart.cis_qtl_pairs.18.txt"))

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart_chr1$variant_id,split="_"))))
eqtl_Heart_chr1$chr<-tmp$V1
eqtl_Heart_chr1$loci<-tmp$V2

eqtl_Heart_chr1$index<-paste0(eqtl_Heart_chr1$chr,"-",eqtl_Heart_chr1$loci)
Heart_chr1$index<-paste0(Heart_chr1$chr,"-",Heart_chr1$pig_pos)
same_Heart_chr1<-intersect(Heart_chr1$pig_pos,eqtl_Heart_chr1$loci)


tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart_chr2$variant_id,split="_"))))
eqtl_Heart_chr2$chr<-tmp$V1
eqtl_Heart_chr2$loci<-tmp$V2

eqtl_Heart_chr2$index<-paste0(eqtl_Heart_chr2$chr,"-",eqtl_Heart_chr2$loci)
Heart_chr2$index<-paste0(Heart_chr2$chr,"-",Heart_chr2$pig_pos)
same_Heart_chr2<-intersect(Heart_chr2$pig_pos,eqtl_Heart_chr2$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart_chr3$variant_id,split="_"))))
eqtl_Heart_chr3$chr<-tmp$V1
eqtl_Heart_chr3$loci<-tmp$V2

eqtl_Heart_chr3$index<-paste0(eqtl_Heart_chr3$chr,"-",eqtl_Heart_chr3$loci)
Heart_chr3$index<-paste0(Heart_chr3$chr,"-",Heart_chr3$pig_pos)
same_Heart_chr3<-intersect(Heart_chr3$pig_pos,eqtl_Heart_chr3$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart_chr4$variant_id,split="_"))))
eqtl_Heart_chr4$chr<-tmp$V1
eqtl_Heart_chr4$loci<-tmp$V2

eqtl_Heart_chr4$index<-paste0(eqtl_Heart_chr4$chr,"-",eqtl_Heart_chr4$loci)
Heart_chr4$index<-paste0(Heart_chr4$chr,"-",Heart_chr4$pig_pos)
same_Heart_chr4<-intersect(Heart_chr4$pig_pos,eqtl_Heart_chr4$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart_chr5$variant_id,split="_"))))
eqtl_Heart_chr5$loci<-tmp$V2
same_Heart_chr5<-intersect(Heart_chr5$pig_pos,eqtl_Heart_chr5$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart_chr6$variant_id,split="_"))))
eqtl_Heart_chr6$loci<-tmp$V2
same_Heart_chr6<-intersect(Heart_chr6$pig_pos,eqtl_Heart_chr6$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart_chr7$variant_id,split="_"))))
eqtl_Heart_chr7$loci<-tmp$V2
same_Heart_chr7<-intersect(Heart_chr7$pig_pos,eqtl_Heart_chr7$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart_chr8$variant_id,split="_"))))
eqtl_Heart_chr8$loci<-tmp$V2
same_Heart_chr8<-intersect(Heart_chr8$pig_pos,eqtl_Heart_chr8$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart_chr9$variant_id,split="_"))))
eqtl_Heart_chr9$loci<-tmp$V2
same_Heart_chr9<-intersect(Heart_chr9$pig_pos,eqtl_Heart_chr9$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart_chr10$variant_id,split="_"))))
eqtl_Heart_chr10$loci<-tmp$V2
same_Heart_chr10<-intersect(Heart_chr10$pig_pos,eqtl_Heart_chr10$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart_chr11$variant_id,split="_"))))
eqtl_Heart_chr11$loci<-tmp$V2
same_Heart_chr11<-intersect(Heart_chr11$pig_pos,eqtl_Heart_chr11$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart_chr12$variant_id,split="_"))))
eqtl_Heart_chr12$loci<-tmp$V2
same_Heart_chr12<-intersect(Heart_chr12$pig_pos,eqtl_Heart_chr12$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart_chr13$variant_id,split="_"))))
eqtl_Heart_chr13$loci<-tmp$V2
same_Heart_chr13<-intersect(Heart_chr13$pig_pos,eqtl_Heart_chr13$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart_chr14$variant_id,split="_"))))
eqtl_Heart_chr14$loci<-tmp$V2
same_Heart_chr14<-intersect(Heart_chr14$pig_pos,eqtl_Heart_chr14$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart_chr15$variant_id,split="_"))))
eqtl_Heart_chr15$loci<-tmp$V2
same_Heart_chr15<-intersect(Heart_chr15$pig_pos,eqtl_Heart_chr15$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart_chr16$variant_id,split="_"))))
eqtl_Heart_chr16$loci<-tmp$V2
same_Heart_chr16<-intersect(Heart_chr16$pig_pos,eqtl_Heart_chr16$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart_chr17$variant_id,split="_"))))
eqtl_Heart_chr17$loci<-tmp$V2
same_Heart_chr17<-intersect(Heart_chr17$pig_pos,eqtl_Heart_chr17$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Heart_chr18$variant_id,split="_"))))
eqtl_Heart_chr18$loci<-tmp$V2
same_Heart_chr18<-intersect(Heart_chr18$pig_pos,eqtl_Heart_chr18$loci)

chr1_Heart_genes<-NULL
chr1_Heart_overloci<-NULL
if(length(same_Heart_chr1!=0)){
  for(i in 1:length(same_Heart_chr1)){
    a<-eqtl_Heart_chr1$phenotype_id[grep(same_Heart_chr1[i],eqtl_Heart_chr1$loci)]
    b<-eqtl_Heart_chr1[grep(same_Heart_chr1[i],eqtl_Heart_chr1$loci),]
    chr1_Heart_genes<-c(chr1_Heart_genes,a)
    chr1_Heart_overloci<-rbind(chr1_Heart_overloci,b)
  }
}

chr2_Heart_genes<-NULL
chr2_Heart_overloci<-NULL
if(length(same_Heart_chr2!=0)){
  for(i in 1:length(same_Heart_chr2)){
    a<-eqtl_Heart_chr2$phenotype_id[grep(same_Heart_chr2[i],eqtl_Heart_chr2$loci)]
    b<-eqtl_Heart_chr2[grep(same_Heart_chr2[i],eqtl_Heart_chr2$loci),]
    chr2_Heart_genes<-c(chr2_Heart_genes,a)
    chr2_Heart_overloci<-rbind(chr2_Heart_overloci,b)
  }
}

chr3_Heart_genes<-NULL
chr3_Heart_overloci<-NULL
if(length(same_Heart_chr3!=0)){
  for(i in 1:length(same_Heart_chr3)){
    a<-eqtl_Heart_chr3$phenotype_id[grep(same_Heart_chr3[i],eqtl_Heart_chr3$loci)]
    b<-eqtl_Heart_chr3[grep(same_Heart_chr3[i],eqtl_Heart_chr3$loci),]
    chr3_Heart_genes<-c(chr3_Heart_genes,a)
    chr3_Heart_overloci<-rbind(chr3_Heart_overloci,b)
  }
}

chr4_Heart_genes<-NULL
chr4_Heart_overloci<-NULL
if(length(same_Heart_chr4!=0)){
  for(i in 1:length(same_Heart_chr4)){
    a<-eqtl_Heart_chr4$phenotype_id[grep(same_Heart_chr4[i],eqtl_Heart_chr4$loci)]
    b<-eqtl_Heart_chr4[grep(same_Heart_chr4[i],eqtl_Heart_chr4$loci),]
    chr4_Heart_genes<-c(chr4_Heart_genes,a)
    chr4_Heart_overloci<-rbind(chr4_Heart_overloci,b)
  }
}

chr5_Heart_genes<-NULL
chr5_Heart_overloci<-NULL
if(length(same_Heart_chr5!=0)){
  for(i in 1:length(same_Heart_chr5)){
    a<-eqtl_Heart_chr5$phenotype_id[grep(same_Heart_chr5[i],eqtl_Heart_chr5$loci)]
    b<-eqtl_Heart_chr5[grep(same_Heart_chr5[i],eqtl_Heart_chr5$loci),]
    chr5_Heart_genes<-c(chr5_Heart_genes,a)
    chr5_Heart_overloci<-rbind(chr5_Heart_overloci,b)
  }
}

chr6_Heart_genes<-NULL
chr6_Heart_overloci<-NULL
if(length(same_Heart_chr6!=0)){
  for(i in 1:length(same_Heart_chr6)){
    a<-eqtl_Heart_chr6$phenotype_id[grep(same_Heart_chr6[i],eqtl_Heart_chr6$loci)]
    b<-eqtl_Heart_chr6[grep(same_Heart_chr6[i],eqtl_Heart_chr6$loci),]
    chr6_Heart_genes<-c(chr6_Heart_genes,a)
    chr6_Heart_overloci<-rbind(chr6_Heart_overloci,b)
  }
}

chr7_Heart_genes<-NULL
chr7_Heart_overloci<-NULL
if(length(same_Heart_chr7!=0)){
  for(i in 1:length(same_Heart_chr7)){
    a<-eqtl_Heart_chr7$phenotype_id[grep(same_Heart_chr7[i],eqtl_Heart_chr7$loci)]
    b<-eqtl_Heart_chr7[grep(same_Heart_chr7[i],eqtl_Heart_chr7$loci),]
    chr7_Heart_genes<-c(chr7_Heart_genes,a)
    chr7_Heart_overloci<-rbind(chr7_Heart_overloci,b)
  }
}

chr8_Heart_genes<-NULL
chr8_Heart_overloci<-NULL
if(length(same_Heart_chr8!=0)){
  for(i in 1:length(same_Heart_chr8)){
    a<-eqtl_Heart_chr8$phenotype_id[grep(same_Heart_chr8[i],eqtl_Heart_chr8$loci)]
    b<-eqtl_Heart_chr8[grep(same_Heart_chr8[i],eqtl_Heart_chr8$loci),]
    chr8_Heart_genes<-c(chr8_Heart_genes,a)
    chr8_Heart_overloci<-rbind(chr8_Heart_overloci,b)
  }
}

if(length(same_Heart_chr9!=0)){
  chr9_Heart_genes<-NULL
  chr9_Heart_overloci<-NULL
  for(i in 1:length(same_Heart_chr9)){
    a<-eqtl_Heart_chr9$phenotype_id[grep(same_Heart_chr9[i],eqtl_Heart_chr9$loci)]
    b<-eqtl_Heart_chr9[grep(same_Heart_chr9[i],eqtl_Heart_chr9$loci),]
    chr9_Heart_genes<-c(chr9_Heart_genes,a)
    chr9_Heart_overloci<-rbind(chr9_Heart_overloci,b)
  }
}

chr10_Heart_genes<-NULL
chr10_Heart_overloci<-NULL
if(length(same_Heart_chr10!=0)){
  for(i in 1:length(same_Heart_chr10)){
    a<-eqtl_Heart_chr10$phenotype_id[grep(same_Heart_chr10[i],eqtl_Heart_chr10$loci)]
    b<-eqtl_Heart_chr10[grep(same_Heart_chr10[i],eqtl_Heart_chr10$loci),]
    chr10_Heart_genes<-c(chr10_Heart_genes,a)
    chr10_Heart_overloci<-rbind(chr10_Heart_overloci,b)
  }
}

chr11_Heart_genes<-NULL
chr11_Heart_overloci<-NULL
if(length(same_Heart_chr11!=0)){
  for(i in 1:length(same_Heart_chr11)){
    a<-eqtl_Heart_chr11$phenotype_id[grep(same_Heart_chr11[i],eqtl_Heart_chr11$loci)]
    b<-eqtl_Heart_chr11[grep(same_Heart_chr11[i],eqtl_Heart_chr11$loci),]
    chr11_Heart_genes<-c(chr11_Heart_genes,a)
    chr11_Heart_overloci<-rbind(chr11_Heart_overloci,b)
  }
}

chr12_Heart_genes<-NULL
chr12_Heart_overloci<-NULL
if(length(same_Heart_chr12!=0)){
  for(i in 1:length(same_Heart_chr12)){
    a<-eqtl_Heart_chr12$phenotype_id[grep(same_Heart_chr12[i],eqtl_Heart_chr12$loci)]
    b<-eqtl_Heart_chr12[grep(same_Heart_chr12[i],eqtl_Heart_chr12$loci),]
    chr12_Heart_genes<-c(chr12_Heart_genes,a)
    chr12_Heart_overloci<-rbind(chr12_Heart_overloci,b)
  }
}

chr13_Heart_genes<-NULL
chr13_Heart_overloci<-NULL
if(length(same_Heart_chr13!=0)){
  for(i in 1:length(same_Heart_chr13)){
    a<-eqtl_Heart_chr13$phenotype_id[grep(same_Heart_chr13[i],eqtl_Heart_chr13$loci)]
    b<-eqtl_Heart_chr13[grep(same_Heart_chr13[i],eqtl_Heart_chr13$loci),]
    chr13_Heart_genes<-c(chr13_Heart_genes,a)
    chr13_Heart_overloci<-rbind(chr13_Heart_overloci,b)
  }
}

chr14_Heart_genes<-NULL
chr14_Heart_overloci<-NULL
if(length(same_Heart_chr14!=0)){
  for(i in 1:length(same_Heart_chr14)){
    a<-eqtl_Heart_chr14$phenotype_id[grep(same_Heart_chr14[i],eqtl_Heart_chr14$loci)]
    b<-eqtl_Heart_chr14[grep(same_Heart_chr14[i],eqtl_Heart_chr14$loci),]
    chr14_Heart_genes<-c(chr14_Heart_genes,a)
    chr14_Heart_overloci<-rbind(chr14_Heart_overloci,b)
  }
}

chr15_Heart_genes<-NULL
chr15_Heart_overloci<-NULL
if(length(same_Heart_chr15!=0)){
  for(i in 1:length(same_Heart_chr15)){
    a<-eqtl_Heart_chr15$phenotype_id[grep(same_Heart_chr15[i],eqtl_Heart_chr15$loci)]
    b<-eqtl_Heart_chr15[grep(same_Heart_chr15[i],eqtl_Heart_chr15$loci),]
    chr15_Heart_genes<-c(chr15_Heart_genes,a)
    chr15_Heart_overloci<-rbind(chr15_Heart_overloci,b)
  }
}

chr16_Heart_genes<-NULL
chr16_Heart_overloci<-NULL
if(length(same_Heart_chr16!=0)){
  for(i in 1:length(same_Heart_chr16)){
    a<-eqtl_Heart_chr16$phenotype_id[grep(same_Heart_chr16[i],eqtl_Heart_chr16$loci)]
    b<-eqtl_Heart_chr16[grep(same_Heart_chr16[i],eqtl_Heart_chr16$loci),]
    chr16_Heart_genes<-c(chr16_Heart_genes,a)
    chr16_Heart_overloci<-rbind(chr16_Heart_overloci,b)
  }
}

chr17_Heart_genes<-NULL
chr17_Heart_overloci<-NULL
if(length(same_Heart_chr17!=0)){
  for(i in 1:length(same_Heart_chr17)){
    a<-eqtl_Heart_chr17$phenotype_id[grep(same_Heart_chr17[i],eqtl_Heart_chr17$loci)]
    b<-eqtl_Heart_chr17[grep(same_Heart_chr17[i],eqtl_Heart_chr17$loci),]
    chr17_Heart_genes<-c(chr17_Heart_genes,a)
    chr17_Heart_overloci<-rbind(chr17_Heart_overloci,b)
  }
}

chr18_Heart_genes<-NULL
chr18_Heart_overloci<-NULL
if(length(same_Heart_chr18!=0)){
  for(i in 1:length(same_Heart_chr18)){
    a<-eqtl_Heart_chr18$phenotype_id[grep(same_Heart_chr18[i],eqtl_Heart_chr18$loci)]
    b<-eqtl_Heart_chr18[grep(same_Heart_chr18[i],eqtl_Heart_chr18$loci),]
    chr18_Heart_genes<-c(chr18_Heart_genes,a)
    chr18_Heart_overloci<-rbind(chr18_Heart_overloci,b)
  }
}
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
one2one_pig<-annotation$`Pig gene stable ID`[which(annotation$`Pig homology type`=="ortholog_one2one")]

chr1_Heart_overloci_one2one<-chr1_Heart_overloci[match(intersect(chr1_Heart_genes,one2one_pig),chr1_Heart_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Heart_chr1$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Heart_chr1$hum_chr<-tmp$V1
Heart_chr1$hum_loci<-tmp$V3
FM_Heart_hum$index<-paste0(FM_Heart_hum$chr,"_",FM_Heart_hum$variant_pos)
Heart_chr1$index<-paste0("chr",Heart_chr1$hum_chr,"_",Heart_chr1$hum_loci)
chr1_Heart_hum<-NULL
for(i in 1:length(same_Heart_chr1)){
  a<-Heart_chr1[match(same_Heart_chr1[i],Heart_chr1$pig_pos),]
  b<-FM_Heart_hum[match(intersect(a$index,FM_Heart_hum$index),FM_Heart_hum$index),]
  chr1_Heart_hum<-rbind(chr1_Heart_hum,b)
}
chr1_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr1_Heart_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr1_Heart_hum_one2one<-chr1_Heart_hum[match(intersect(chr1_pig2hum_one2one,chr1_Heart_hum$gene_id),chr1_Heart_hum$gene_id),]

chr2_Heart_overloci_one2one<-chr2_Heart_overloci[match(intersect(chr2_Heart_genes,one2one_pig),chr2_Heart_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Heart_chr2$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Heart_chr2$hum_chr<-tmp$V1
Heart_chr2$hum_loci<-tmp$V3
FM_Heart_hum$index<-paste0(FM_Heart_hum$chr,"_",FM_Heart_hum$variant_pos)
Heart_chr2$index<-paste0("chr",Heart_chr2$hum_chr,"_",Heart_chr2$hum_loci)
chr2_Heart_hum<-NULL
for(i in 1:length(same_Heart_chr2)){
  a<-Heart_chr2[match(same_Heart_chr2[i],Heart_chr2$pig_pos),]
  b<-FM_Heart_hum[match(intersect(a$index,FM_Heart_hum$index),FM_Heart_hum$index),]
  chr2_Heart_hum<-rbind(chr2_Heart_hum,b)
}
chr2_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr2_Heart_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr2_Heart_hum_one2one<-chr2_Heart_hum[match(intersect(chr2_pig2hum_one2one,chr2_Heart_hum$gene_id),chr2_Heart_hum$gene_id),]

chr3_Heart_overloci_one2one<-chr3_Heart_overloci[match(intersect(chr3_Heart_genes,one2one_pig),chr3_Heart_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Heart_chr3$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Heart_chr3$hum_chr<-tmp$V1
Heart_chr3$hum_loci<-tmp$V3
FM_Heart_hum$index<-paste0(FM_Heart_hum$chr,"_",FM_Heart_hum$variant_pos)
Heart_chr3$index<-paste0("chr",Heart_chr3$hum_chr,"_",Heart_chr3$hum_loci)
chr3_Heart_hum<-NULL
for(i in 1:length(same_Heart_chr3)){
  a<-Heart_chr3[match(same_Heart_chr3[i],Heart_chr3$pig_pos),]
  b<-FM_Heart_hum[match(intersect(a$index,FM_Heart_hum$index),FM_Heart_hum$index),]
  chr3_Heart_hum<-rbind(chr3_Heart_hum,b)
}
chr3_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr3_Heart_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr3_Heart_hum_one2one<-chr3_Heart_hum[match(intersect(chr3_pig2hum_one2one,chr3_Heart_hum$gene_id),chr3_Heart_hum$gene_id),]

chr4_Heart_overloci_one2one<-chr4_Heart_overloci[match(intersect(chr4_Heart_genes,one2one_pig),chr4_Heart_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Heart_chr4$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Heart_chr4$hum_chr<-tmp$V1
Heart_chr4$hum_loci<-tmp$V3
FM_Heart_hum$index<-paste0(FM_Heart_hum$chr,"_",FM_Heart_hum$variant_pos)
Heart_chr4$index<-paste0("chr",Heart_chr4$hum_chr,"_",Heart_chr4$hum_loci)
chr4_Heart_hum<-NULL
for(i in 1:length(same_Heart_chr4)){
  a<-Heart_chr4[match(same_Heart_chr4[i],Heart_chr4$pig_pos),]
  b<-FM_Heart_hum[match(intersect(a$index,FM_Heart_hum$index),FM_Heart_hum$index),]
  chr4_Heart_hum<-rbind(chr4_Heart_hum,b)
}
chr4_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr4_Heart_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr4_Heart_hum_one2one<-chr4_Heart_hum[match(intersect(chr4_pig2hum_one2one,chr4_Heart_hum$gene_id),chr4_Heart_hum$gene_id),]

chr5_Heart_overloci_one2one<-chr5_Heart_overloci[match(intersect(chr5_Heart_genes,one2one_pig),chr5_Heart_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Heart_chr5$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Heart_chr5$hum_chr<-tmp$V1
Heart_chr5$hum_loci<-tmp$V3
FM_Heart_hum$index<-paste0(FM_Heart_hum$chr,"_",FM_Heart_hum$variant_pos)
Heart_chr5$index<-paste0("chr",Heart_chr5$hum_chr,"_",Heart_chr5$hum_loci)
chr5_Heart_hum<-NULL
for(i in 1:length(same_Heart_chr5)){
  a<-Heart_chr5[match(same_Heart_chr5[i],Heart_chr5$pig_pos),]
  b<-FM_Heart_hum[match(intersect(a$index,FM_Heart_hum$index),FM_Heart_hum$index),]
  chr5_Heart_hum<-rbind(chr5_Heart_hum,b)
}
chr5_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr5_Heart_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr5_Heart_hum_one2one<-chr5_Heart_hum[match(intersect(chr5_pig2hum_one2one,chr5_Heart_hum$gene_id),chr5_Heart_hum$gene_id),]

chr6_Heart_overloci_one2one<-chr6_Heart_overloci[match(intersect(chr6_Heart_genes,one2one_pig),chr6_Heart_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Heart_chr6$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Heart_chr6$hum_chr<-tmp$V1
Heart_chr6$hum_loci<-tmp$V3
FM_Heart_hum$index<-paste0(FM_Heart_hum$chr,"_",FM_Heart_hum$variant_pos)
Heart_chr6$index<-paste0("chr",Heart_chr6$hum_chr,"_",Heart_chr6$hum_loci)
chr6_Heart_hum<-NULL
for(i in 1:length(same_Heart_chr6)){
  a<-Heart_chr6[match(same_Heart_chr6[i],Heart_chr6$pig_pos),]
  b<-FM_Heart_hum[match(intersect(a$index,FM_Heart_hum$index),FM_Heart_hum$index),]
  chr6_Heart_hum<-rbind(chr6_Heart_hum,b)
}
chr6_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr6_Heart_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr6_Heart_hum_one2one<-chr6_Heart_hum[match(intersect(chr6_pig2hum_one2one,chr6_Heart_hum$gene_id),chr6_Heart_hum$gene_id),]

chr7_Heart_overloci_one2one<-chr7_Heart_overloci[match(intersect(chr7_Heart_genes,one2one_pig),chr7_Heart_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Heart_chr7$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Heart_chr7$hum_chr<-tmp$V1
Heart_chr7$hum_loci<-tmp$V3
FM_Heart_hum$index<-paste0(FM_Heart_hum$chr,"_",FM_Heart_hum$variant_pos)
Heart_chr7$index<-paste0("chr",Heart_chr7$hum_chr,"_",Heart_chr7$hum_loci)
chr7_Heart_hum<-NULL
for(i in 1:length(same_Heart_chr7)){
  a<-Heart_chr7[match(same_Heart_chr7[i],Heart_chr7$pig_pos),]
  b<-FM_Heart_hum[match(intersect(a$index,FM_Heart_hum$index),FM_Heart_hum$index),]
  chr7_Heart_hum<-rbind(chr7_Heart_hum,b)
}
chr7_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr7_Heart_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr7_Heart_hum_one2one<-chr7_Heart_hum[match(intersect(chr7_pig2hum_one2one,chr7_Heart_hum$gene_id),chr7_Heart_hum$gene_id),]

chr8_Heart_overloci_one2one<-chr8_Heart_overloci[match(intersect(chr8_Heart_genes,one2one_pig),chr8_Heart_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Heart_chr8$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Heart_chr8$hum_chr<-tmp$V1
Heart_chr8$hum_loci<-tmp$V3
FM_Heart_hum$index<-paste0(FM_Heart_hum$chr,"_",FM_Heart_hum$variant_pos)
Heart_chr8$index<-paste0("chr",Heart_chr8$hum_chr,"_",Heart_chr8$hum_loci)
chr8_Heart_hum<-NULL
for(i in 1:length(same_Heart_chr8)){
  a<-Heart_chr8[match(same_Heart_chr8[i],Heart_chr8$pig_pos),]
  b<-FM_Heart_hum[match(intersect(a$index,FM_Heart_hum$index),FM_Heart_hum$index),]
  chr8_Heart_hum<-rbind(chr8_Heart_hum,b)
}
chr8_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr8_Heart_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr8_Heart_hum_one2one<-chr8_Heart_hum[match(intersect(chr8_pig2hum_one2one,chr8_Heart_hum$gene_id),chr8_Heart_hum$gene_id),]

chr9_Heart_overloci_one2one<-chr9_Heart_overloci[match(intersect(chr9_Heart_genes,one2one_pig),chr9_Heart_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Heart_chr9$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Heart_chr9$hum_chr<-tmp$V1
Heart_chr9$hum_loci<-tmp$V3
FM_Heart_hum$index<-paste0(FM_Heart_hum$chr,"_",FM_Heart_hum$variant_pos)
Heart_chr9$index<-paste0("chr",Heart_chr9$hum_chr,"_",Heart_chr9$hum_loci)
chr9_Heart_hum<-NULL
for(i in 1:length(same_Heart_chr9)){
  a<-Heart_chr9[match(same_Heart_chr9[i],Heart_chr9$pig_pos),]
  b<-FM_Heart_hum[match(intersect(a$index,FM_Heart_hum$index),FM_Heart_hum$index),]
  chr9_Heart_hum<-rbind(chr9_Heart_hum,b)
}
chr9_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr9_Heart_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr9_Heart_hum_one2one<-chr9_Heart_hum[match(intersect(chr9_pig2hum_one2one,chr9_Heart_hum$gene_id),chr9_Heart_hum$gene_id),]

chr10_Heart_overloci_one2one<-chr10_Heart_overloci[match(intersect(chr10_Heart_genes,one2one_pig),chr10_Heart_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Heart_chr10$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Heart_chr10$hum_chr<-tmp$V1
Heart_chr10$hum_loci<-tmp$V3
FM_Heart_hum$index<-paste0(FM_Heart_hum$chr,"_",FM_Heart_hum$variant_pos)
Heart_chr10$index<-paste0("chr",Heart_chr10$hum_chr,"_",Heart_chr10$hum_loci)
chr10_Heart_hum<-NULL
for(i in 1:length(same_Heart_chr10)){
  a<-Heart_chr10[match(same_Heart_chr10[i],Heart_chr10$pig_pos),]
  b<-FM_Heart_hum[match(intersect(a$index,FM_Heart_hum$index),FM_Heart_hum$index),]
  chr10_Heart_hum<-rbind(chr10_Heart_hum,b)
}
chr10_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr10_Heart_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr10_Heart_hum_one2one<-chr10_Heart_hum[match(intersect(chr10_pig2hum_one2one,chr10_Heart_hum$gene_id),chr10_Heart_hum$gene_id),]

chr11_Heart_overloci_one2one<-chr11_Heart_overloci[match(intersect(chr11_Heart_genes,one2one_pig),chr11_Heart_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Heart_chr11$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Heart_chr11$hum_chr<-tmp$V1
Heart_chr11$hum_loci<-tmp$V3
FM_Heart_hum$index<-paste0(FM_Heart_hum$chr,"_",FM_Heart_hum$variant_pos)
Heart_chr11$index<-paste0("chr",Heart_chr11$hum_chr,"_",Heart_chr11$hum_loci)
chr11_Heart_hum<-NULL
for(i in 1:length(same_Heart_chr11)){
  a<-Heart_chr11[match(same_Heart_chr11[i],Heart_chr11$pig_pos),]
  b<-FM_Heart_hum[match(intersect(a$index,FM_Heart_hum$index),FM_Heart_hum$index),]
  chr11_Heart_hum<-rbind(chr11_Heart_hum,b)
}
chr11_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr11_Heart_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr11_Heart_hum_one2one<-chr11_Heart_hum[match(intersect(chr11_pig2hum_one2one,chr11_Heart_hum$gene_id),chr11_Heart_hum$gene_id),]

chr12_Heart_overloci_one2one<-chr12_Heart_overloci[match(intersect(chr12_Heart_genes,one2one_pig),chr12_Heart_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Heart_chr12$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Heart_chr12$hum_chr<-tmp$V1
Heart_chr12$hum_loci<-tmp$V3
FM_Heart_hum$index<-paste0(FM_Heart_hum$chr,"_",FM_Heart_hum$variant_pos)
Heart_chr12$index<-paste0("chr",Heart_chr12$hum_chr,"_",Heart_chr12$hum_loci)
chr12_Heart_hum<-NULL
for(i in 1:length(same_Heart_chr12)){
  a<-Heart_chr12[match(same_Heart_chr12[i],Heart_chr12$pig_pos),]
  b<-FM_Heart_hum[match(intersect(a$index,FM_Heart_hum$index),FM_Heart_hum$index),]
  chr12_Heart_hum<-rbind(chr12_Heart_hum,b)
}
chr12_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr12_Heart_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr12_Heart_hum_one2one<-chr12_Heart_hum[match(intersect(chr12_pig2hum_one2one,chr12_Heart_hum$gene_id),chr12_Heart_hum$gene_id),]

chr13_Heart_overloci_one2one<-chr13_Heart_overloci[match(intersect(chr13_Heart_genes,one2one_pig),chr13_Heart_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Heart_chr13$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Heart_chr13$hum_chr<-tmp$V1
Heart_chr13$hum_loci<-tmp$V3
FM_Heart_hum$index<-paste0(FM_Heart_hum$chr,"_",FM_Heart_hum$variant_pos)
Heart_chr13$index<-paste0("chr",Heart_chr13$hum_chr,"_",Heart_chr13$hum_loci)
chr13_Heart_hum<-NULL
for(i in 1:length(same_Heart_chr13)){
  a<-Heart_chr13[match(same_Heart_chr13[i],Heart_chr13$pig_pos),]
  b<-FM_Heart_hum[match(intersect(a$index,FM_Heart_hum$index),FM_Heart_hum$index),]
  chr13_Heart_hum<-rbind(chr13_Heart_hum,b)
}
chr13_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr13_Heart_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr13_Heart_hum_one2one<-chr13_Heart_hum[match(intersect(chr13_pig2hum_one2one,chr13_Heart_hum$gene_id),chr13_Heart_hum$gene_id),]

chr14_Heart_overloci_one2one<-chr14_Heart_overloci[match(intersect(chr14_Heart_genes,one2one_pig),chr14_Heart_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Heart_chr14$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Heart_chr14$hum_chr<-tmp$V1
Heart_chr14$hum_loci<-tmp$V3
FM_Heart_hum$index<-paste0(FM_Heart_hum$chr,"_",FM_Heart_hum$variant_pos)
Heart_chr14$index<-paste0("chr",Heart_chr14$hum_chr,"_",Heart_chr14$hum_loci)
chr14_Heart_hum<-NULL
for(i in 1:length(same_Heart_chr14)){
  a<-Heart_chr14[match(same_Heart_chr14[i],Heart_chr14$pig_pos),]
  b<-FM_Heart_hum[match(intersect(a$index,FM_Heart_hum$index),FM_Heart_hum$index),]
  chr14_Heart_hum<-rbind(chr14_Heart_hum,b)
}
chr14_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr14_Heart_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr14_Heart_hum_one2one<-chr14_Heart_hum[match(intersect(chr14_pig2hum_one2one,chr14_Heart_hum$gene_id),chr14_Heart_hum$gene_id),]

chr15_Heart_overloci_one2one<-chr15_Heart_overloci[match(intersect(chr15_Heart_genes,one2one_pig),chr15_Heart_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Heart_chr15$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Heart_chr15$hum_chr<-tmp$V1
Heart_chr15$hum_loci<-tmp$V3
FM_Heart_hum$index<-paste0(FM_Heart_hum$chr,"_",FM_Heart_hum$variant_pos)
Heart_chr15$index<-paste0("chr",Heart_chr15$hum_chr,"_",Heart_chr15$hum_loci)
chr15_Heart_hum<-NULL
for(i in 1:length(same_Heart_chr15)){
  a<-Heart_chr15[match(same_Heart_chr15[i],Heart_chr15$pig_pos),]
  b<-FM_Heart_hum[match(intersect(a$index,FM_Heart_hum$index),FM_Heart_hum$index),]
  chr15_Heart_hum<-rbind(chr15_Heart_hum,b)
}
chr15_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr15_Heart_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr15_Heart_hum_one2one<-chr15_Heart_hum[match(intersect(chr15_pig2hum_one2one,chr15_Heart_hum$gene_id),chr15_Heart_hum$gene_id),]

chr16_Heart_overloci_one2one<-chr16_Heart_overloci[match(intersect(chr16_Heart_genes,one2one_pig),chr16_Heart_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Heart_chr16$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Heart_chr16$hum_chr<-tmp$V1
Heart_chr16$hum_loci<-tmp$V3
FM_Heart_hum$index<-paste0(FM_Heart_hum$chr,"_",FM_Heart_hum$variant_pos)
Heart_chr16$index<-paste0("chr",Heart_chr16$hum_chr,"_",Heart_chr16$hum_loci)
chr16_Heart_hum<-NULL
for(i in 1:length(same_Heart_chr16)){
  a<-Heart_chr16[match(same_Heart_chr16[i],Heart_chr16$pig_pos),]
  b<-FM_Heart_hum[match(intersect(a$index,FM_Heart_hum$index),FM_Heart_hum$index),]
  chr16_Heart_hum<-rbind(chr16_Heart_hum,b)
}
chr16_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr16_Heart_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr16_Heart_hum_one2one<-chr16_Heart_hum[match(intersect(chr16_pig2hum_one2one,chr16_Heart_hum$gene_id),chr16_Heart_hum$gene_id),]

chr17_Heart_overloci_one2one<-chr17_Heart_overloci[match(intersect(chr17_Heart_genes,one2one_pig),chr17_Heart_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Heart_chr17$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Heart_chr17$hum_chr<-tmp$V1
Heart_chr17$hum_loci<-tmp$V3
FM_Heart_hum$index<-paste0(FM_Heart_hum$chr,"_",FM_Heart_hum$variant_pos)
Heart_chr17$index<-paste0("chr",Heart_chr17$hum_chr,"_",Heart_chr17$hum_loci)
chr17_Heart_hum<-NULL
for(i in 1:length(same_Heart_chr17)){
  a<-Heart_chr17[match(same_Heart_chr17[i],Heart_chr17$pig_pos),]
  b<-FM_Heart_hum[match(intersect(a$index,FM_Heart_hum$index),FM_Heart_hum$index),]
  chr17_Heart_hum<-rbind(chr17_Heart_hum,b)
}
chr17_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr17_Heart_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr17_Heart_hum_one2one<-chr17_Heart_hum[match(intersect(chr17_pig2hum_one2one,chr17_Heart_hum$gene_id),chr17_Heart_hum$gene_id),]

chr18_Heart_overloci_one2one<-chr18_Heart_overloci[match(intersect(chr18_Heart_genes,one2one_pig),chr18_Heart_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Heart_chr18$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Heart_chr18$hum_chr<-tmp$V1
Heart_chr18$hum_loci<-tmp$V3
FM_Heart_hum$index<-paste0(FM_Heart_hum$chr,"_",FM_Heart_hum$variant_pos)
Heart_chr18$index<-paste0("chr",Heart_chr18$hum_chr,"_",Heart_chr18$hum_loci)
chr18_Heart_hum<-NULL
for(i in 1:length(same_Heart_chr18)){
  a<-Heart_chr18[match(same_Heart_chr18[i],Heart_chr18$pig_pos),]
  b<-FM_Heart_hum[match(intersect(a$index,FM_Heart_hum$index),FM_Heart_hum$index),]
  chr18_Heart_hum<-rbind(chr18_Heart_hum,b)
}
chr18_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr18_Heart_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr18_Heart_hum_one2one<-chr18_Heart_hum[match(intersect(chr18_pig2hum_one2one,chr18_Heart_hum$gene_id),chr18_Heart_hum$gene_id),]

Heart_one2one_SNP_hum<-rbind(chr1_Heart_hum_one2one,chr2_Heart_hum_one2one,chr3_Heart_hum_one2one,chr4_Heart_hum_one2one,chr5_Heart_hum_one2one,
                             chr6_Heart_hum_one2one,chr7_Heart_hum_one2one,chr8_Heart_hum_one2one,chr9_Heart_hum_one2one,chr10_Heart_hum_one2one,
                             chr11_Heart_hum_one2one,chr12_Heart_hum_one2one,chr13_Heart_hum_one2one,chr14_Heart_hum_one2one,chr15_Heart_hum_one2one,
                             chr16_Heart_hum_one2one,chr17_Heart_hum_one2one,chr18_Heart_hum_one2one)

chr1_Heart_pig_one2one<-chr1_Heart_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr1_Heart_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr1_Heart_overloci_one2one$phenotype_id),1:9]
chr2_Heart_pig_one2one<-chr2_Heart_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr2_Heart_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr2_Heart_overloci_one2one$phenotype_id),1:9]
chr3_Heart_pig_one2one<-chr3_Heart_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr3_Heart_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr3_Heart_overloci_one2one$phenotype_id),1:9]
chr4_Heart_pig_one2one<-chr4_Heart_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr4_Heart_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr4_Heart_overloci_one2one$phenotype_id),1:9]
chr5_Heart_pig_one2one<-chr5_Heart_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr5_Heart_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr5_Heart_overloci_one2one$phenotype_id),1:9]
chr6_Heart_pig_one2one<-chr6_Heart_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr6_Heart_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr6_Heart_overloci_one2one$phenotype_id),1:9]
chr7_Heart_pig_one2one<-chr7_Heart_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr7_Heart_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr7_Heart_overloci_one2one$phenotype_id),1:9]
chr8_Heart_pig_one2one<-chr8_Heart_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr8_Heart_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr8_Heart_overloci_one2one$phenotype_id),1:9]
chr9_Heart_pig_one2one<-chr9_Heart_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr9_Heart_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr9_Heart_overloci_one2one$phenotype_id),1:9]
chr10_Heart_pig_one2one<-chr10_Heart_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr10_Heart_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr10_Heart_overloci_one2one$phenotype_id),1:9]
chr11_Heart_pig_one2one<-chr11_Heart_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr11_Heart_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr11_Heart_overloci_one2one$phenotype_id),1:9]
chr12_Heart_pig_one2one<-chr12_Heart_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr12_Heart_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr12_Heart_overloci_one2one$phenotype_id),1:9]
chr13_Heart_pig_one2one<-chr13_Heart_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr13_Heart_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr13_Heart_overloci_one2one$phenotype_id),1:9]
chr14_Heart_pig_one2one<-chr14_Heart_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr14_Heart_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr14_Heart_overloci_one2one$phenotype_id),1:9]
chr15_Heart_pig_one2one<-chr15_Heart_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr15_Heart_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr15_Heart_overloci_one2one$phenotype_id),1:9]
chr16_Heart_pig_one2one<-chr16_Heart_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr16_Heart_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr16_Heart_overloci_one2one$phenotype_id),1:9]
chr17_Heart_pig_one2one<-chr17_Heart_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr17_Heart_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr17_Heart_overloci_one2one$phenotype_id),1:9]
chr18_Heart_pig_one2one<-chr18_Heart_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr18_Heart_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr18_Heart_overloci_one2one$phenotype_id),1:9]

Heart_one2one_SNP_pig<-rbind(chr1_Heart_pig_one2one,chr2_Heart_pig_one2one,chr3_Heart_pig_one2one,chr4_Heart_pig_one2one,chr5_Heart_pig_one2one,
                             chr6_Heart_pig_one2one,chr7_Heart_pig_one2one,chr8_Heart_pig_one2one,chr9_Heart_pig_one2one,chr10_Heart_pig_one2one,
                             chr11_Heart_pig_one2one,chr12_Heart_pig_one2one,chr13_Heart_pig_one2one,chr14_Heart_pig_one2one,chr15_Heart_pig_one2one,
                             chr16_Heart_pig_one2one,chr17_Heart_pig_one2one,chr18_Heart_pig_one2one)

Heart_SNP_sum<-array(NA,dim=c(nrow(Heart_one2one_SNP_hum),2))
colnames(Heart_SNP_sum)<-c("Human","Pig")
Heart_SNP_sum<-as.data.frame(Heart_SNP_sum)
Heart_SNP_sum$Human<-Heart_one2one_SNP_hum$slope / Heart_one2one_SNP_hum$slope_se
Heart_SNP_sum$Pig<-Heart_one2one_SNP_pig$slope / Heart_one2one_SNP_pig$slope_se
cor<-cor(abs(Heart_SNP_sum$Human),abs(Heart_SNP_sum$Pig))
p_val<-t.test(abs(Heart_SNP_sum$Human),abs(Heart_SNP_sum$Pig))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/Heart_SNP_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(Heart_SNP_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The number of eGenes expressed in common tissues",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

save(Heart_one2one_SNP_hum,Heart_one2one_SNP_pig,file="/Users/baizhonghao/Downloads/human-pig GTEx/Heart_SNP.Rdata")

#Hypothalamus_SNP_overlaploci#
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

eqtl_Hypothalamus_chr1<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Hypothalamus/Hypothalamus.cis_qtl_pairs.1.txt"))
eqtl_Hypothalamus_chr2<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Hypothalamus/Hypothalamus.cis_qtl_pairs.2.txt"))
eqtl_Hypothalamus_chr3<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Hypothalamus/Hypothalamus.cis_qtl_pairs.3.txt"))
eqtl_Hypothalamus_chr4<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Hypothalamus/Hypothalamus.cis_qtl_pairs.4.txt"))
eqtl_Hypothalamus_chr5<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Hypothalamus/Hypothalamus.cis_qtl_pairs.5.txt"))
eqtl_Hypothalamus_chr6<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Hypothalamus/Hypothalamus.cis_qtl_pairs.6.txt"))
eqtl_Hypothalamus_chr7<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Hypothalamus/Hypothalamus.cis_qtl_pairs.7.txt"))
eqtl_Hypothalamus_chr8<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Hypothalamus/Hypothalamus.cis_qtl_pairs.8.txt"))
eqtl_Hypothalamus_chr9<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Hypothalamus/Hypothalamus.cis_qtl_pairs.9.txt"))
eqtl_Hypothalamus_chr10<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Hypothalamus/Hypothalamus.cis_qtl_pairs.10.txt"))
eqtl_Hypothalamus_chr11<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Hypothalamus/Hypothalamus.cis_qtl_pairs.11.txt"))
eqtl_Hypothalamus_chr12<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Hypothalamus/Hypothalamus.cis_qtl_pairs.12.txt"))
eqtl_Hypothalamus_chr13<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Hypothalamus/Hypothalamus.cis_qtl_pairs.13.txt"))
eqtl_Hypothalamus_chr14<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Hypothalamus/Hypothalamus.cis_qtl_pairs.14.txt"))
eqtl_Hypothalamus_chr15<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Hypothalamus/Hypothalamus.cis_qtl_pairs.15.txt"))
eqtl_Hypothalamus_chr16<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Hypothalamus/Hypothalamus.cis_qtl_pairs.16.txt"))
eqtl_Hypothalamus_chr17<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Hypothalamus/Hypothalamus.cis_qtl_pairs.17.txt"))
eqtl_Hypothalamus_chr18<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Hypothalamus/Hypothalamus.cis_qtl_pairs.18.txt"))

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus_chr1$variant_id,split="_"))))
eqtl_Hypothalamus_chr1$chr<-tmp$V1
eqtl_Hypothalamus_chr1$loci<-tmp$V2

eqtl_Hypothalamus_chr1$index<-paste0(eqtl_Hypothalamus_chr1$chr,"-",eqtl_Hypothalamus_chr1$loci)
Hypothalamus_chr1$index<-paste0(Hypothalamus_chr1$chr,"-",Hypothalamus_chr1$pig_pos)
same_Hypothalamus_chr1<-intersect(Hypothalamus_chr1$pig_pos,eqtl_Hypothalamus_chr1$loci)


tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus_chr2$variant_id,split="_"))))
eqtl_Hypothalamus_chr2$chr<-tmp$V1
eqtl_Hypothalamus_chr2$loci<-tmp$V2

eqtl_Hypothalamus_chr2$index<-paste0(eqtl_Hypothalamus_chr2$chr,"-",eqtl_Hypothalamus_chr2$loci)
Hypothalamus_chr2$index<-paste0(Hypothalamus_chr2$chr,"-",Hypothalamus_chr2$pig_pos)
same_Hypothalamus_chr2<-intersect(Hypothalamus_chr2$pig_pos,eqtl_Hypothalamus_chr2$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus_chr3$variant_id,split="_"))))
eqtl_Hypothalamus_chr3$chr<-tmp$V1
eqtl_Hypothalamus_chr3$loci<-tmp$V2

eqtl_Hypothalamus_chr3$index<-paste0(eqtl_Hypothalamus_chr3$chr,"-",eqtl_Hypothalamus_chr3$loci)
Hypothalamus_chr3$index<-paste0(Hypothalamus_chr3$chr,"-",Hypothalamus_chr3$pig_pos)
same_Hypothalamus_chr3<-intersect(Hypothalamus_chr3$pig_pos,eqtl_Hypothalamus_chr3$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus_chr4$variant_id,split="_"))))
eqtl_Hypothalamus_chr4$chr<-tmp$V1
eqtl_Hypothalamus_chr4$loci<-tmp$V2

eqtl_Hypothalamus_chr4$index<-paste0(eqtl_Hypothalamus_chr4$chr,"-",eqtl_Hypothalamus_chr4$loci)
Hypothalamus_chr4$index<-paste0(Hypothalamus_chr4$chr,"-",Hypothalamus_chr4$pig_pos)
same_Hypothalamus_chr4<-intersect(Hypothalamus_chr4$pig_pos,eqtl_Hypothalamus_chr4$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus_chr5$variant_id,split="_"))))
eqtl_Hypothalamus_chr5$loci<-tmp$V2
same_Hypothalamus_chr5<-intersect(Hypothalamus_chr5$pig_pos,eqtl_Hypothalamus_chr5$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus_chr6$variant_id,split="_"))))
eqtl_Hypothalamus_chr6$loci<-tmp$V2
same_Hypothalamus_chr6<-intersect(Hypothalamus_chr6$pig_pos,eqtl_Hypothalamus_chr6$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus_chr7$variant_id,split="_"))))
eqtl_Hypothalamus_chr7$loci<-tmp$V2
same_Hypothalamus_chr7<-intersect(Hypothalamus_chr7$pig_pos,eqtl_Hypothalamus_chr7$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus_chr8$variant_id,split="_"))))
eqtl_Hypothalamus_chr8$loci<-tmp$V2
same_Hypothalamus_chr8<-intersect(Hypothalamus_chr8$pig_pos,eqtl_Hypothalamus_chr8$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus_chr9$variant_id,split="_"))))
eqtl_Hypothalamus_chr9$loci<-tmp$V2
same_Hypothalamus_chr9<-intersect(Hypothalamus_chr9$pig_pos,eqtl_Hypothalamus_chr9$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus_chr10$variant_id,split="_"))))
eqtl_Hypothalamus_chr10$loci<-tmp$V2
same_Hypothalamus_chr10<-intersect(Hypothalamus_chr10$pig_pos,eqtl_Hypothalamus_chr10$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus_chr11$variant_id,split="_"))))
eqtl_Hypothalamus_chr11$loci<-tmp$V2
same_Hypothalamus_chr11<-intersect(Hypothalamus_chr11$pig_pos,eqtl_Hypothalamus_chr11$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus_chr12$variant_id,split="_"))))
eqtl_Hypothalamus_chr12$loci<-tmp$V2
same_Hypothalamus_chr12<-intersect(Hypothalamus_chr12$pig_pos,eqtl_Hypothalamus_chr12$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus_chr13$variant_id,split="_"))))
eqtl_Hypothalamus_chr13$loci<-tmp$V2
same_Hypothalamus_chr13<-intersect(Hypothalamus_chr13$pig_pos,eqtl_Hypothalamus_chr13$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus_chr14$variant_id,split="_"))))
eqtl_Hypothalamus_chr14$loci<-tmp$V2
same_Hypothalamus_chr14<-intersect(Hypothalamus_chr14$pig_pos,eqtl_Hypothalamus_chr14$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus_chr15$variant_id,split="_"))))
eqtl_Hypothalamus_chr15$loci<-tmp$V2
same_Hypothalamus_chr15<-intersect(Hypothalamus_chr15$pig_pos,eqtl_Hypothalamus_chr15$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus_chr16$variant_id,split="_"))))
eqtl_Hypothalamus_chr16$loci<-tmp$V2
same_Hypothalamus_chr16<-intersect(Hypothalamus_chr16$pig_pos,eqtl_Hypothalamus_chr16$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus_chr17$variant_id,split="_"))))
eqtl_Hypothalamus_chr17$loci<-tmp$V2
same_Hypothalamus_chr17<-intersect(Hypothalamus_chr17$pig_pos,eqtl_Hypothalamus_chr17$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Hypothalamus_chr18$variant_id,split="_"))))
eqtl_Hypothalamus_chr18$loci<-tmp$V2
same_Hypothalamus_chr18<-intersect(Hypothalamus_chr18$pig_pos,eqtl_Hypothalamus_chr18$loci)

chr1_Hypothalamus_genes<-NULL
chr1_Hypothalamus_overloci<-NULL
if(length(same_Hypothalamus_chr1!=0)){
  for(i in 1:length(same_Hypothalamus_chr1)){
    a<-eqtl_Hypothalamus_chr1$phenotype_id[grep(same_Hypothalamus_chr1[i],eqtl_Hypothalamus_chr1$loci)]
    b<-eqtl_Hypothalamus_chr1[grep(same_Hypothalamus_chr1[i],eqtl_Hypothalamus_chr1$loci),]
    chr1_Hypothalamus_genes<-c(chr1_Hypothalamus_genes,a)
    chr1_Hypothalamus_overloci<-rbind(chr1_Hypothalamus_overloci,b)
  }
}

chr2_Hypothalamus_genes<-NULL
chr2_Hypothalamus_overloci<-NULL
if(length(same_Hypothalamus_chr2!=0)){
  for(i in 1:length(same_Hypothalamus_chr2)){
    a<-eqtl_Hypothalamus_chr2$phenotype_id[grep(same_Hypothalamus_chr2[i],eqtl_Hypothalamus_chr2$loci)]
    b<-eqtl_Hypothalamus_chr2[grep(same_Hypothalamus_chr2[i],eqtl_Hypothalamus_chr2$loci),]
    chr2_Hypothalamus_genes<-c(chr2_Hypothalamus_genes,a)
    chr2_Hypothalamus_overloci<-rbind(chr2_Hypothalamus_overloci,b)
  }
}

chr3_Hypothalamus_genes<-NULL
chr3_Hypothalamus_overloci<-NULL
if(length(same_Hypothalamus_chr3!=0)){
  for(i in 1:length(same_Hypothalamus_chr3)){
    a<-eqtl_Hypothalamus_chr3$phenotype_id[grep(same_Hypothalamus_chr3[i],eqtl_Hypothalamus_chr3$loci)]
    b<-eqtl_Hypothalamus_chr3[grep(same_Hypothalamus_chr3[i],eqtl_Hypothalamus_chr3$loci),]
    chr3_Hypothalamus_genes<-c(chr3_Hypothalamus_genes,a)
    chr3_Hypothalamus_overloci<-rbind(chr3_Hypothalamus_overloci,b)
  }
}

chr4_Hypothalamus_genes<-NULL
chr4_Hypothalamus_overloci<-NULL
if(length(same_Hypothalamus_chr4!=0)){
  for(i in 1:length(same_Hypothalamus_chr4)){
    a<-eqtl_Hypothalamus_chr4$phenotype_id[grep(same_Hypothalamus_chr4[i],eqtl_Hypothalamus_chr4$loci)]
    b<-eqtl_Hypothalamus_chr4[grep(same_Hypothalamus_chr4[i],eqtl_Hypothalamus_chr4$loci),]
    chr4_Hypothalamus_genes<-c(chr4_Hypothalamus_genes,a)
    chr4_Hypothalamus_overloci<-rbind(chr4_Hypothalamus_overloci,b)
  }
}

chr5_Hypothalamus_genes<-NULL
chr5_Hypothalamus_overloci<-NULL
if(length(same_Hypothalamus_chr5!=0)){
  for(i in 1:length(same_Hypothalamus_chr5)){
    a<-eqtl_Hypothalamus_chr5$phenotype_id[grep(same_Hypothalamus_chr5[i],eqtl_Hypothalamus_chr5$loci)]
    b<-eqtl_Hypothalamus_chr5[grep(same_Hypothalamus_chr5[i],eqtl_Hypothalamus_chr5$loci),]
    chr5_Hypothalamus_genes<-c(chr5_Hypothalamus_genes,a)
    chr5_Hypothalamus_overloci<-rbind(chr5_Hypothalamus_overloci,b)
  }
}

chr6_Hypothalamus_genes<-NULL
chr6_Hypothalamus_overloci<-NULL
if(length(same_Hypothalamus_chr6!=0)){
  for(i in 1:length(same_Hypothalamus_chr6)){
    a<-eqtl_Hypothalamus_chr6$phenotype_id[grep(same_Hypothalamus_chr6[i],eqtl_Hypothalamus_chr6$loci)]
    b<-eqtl_Hypothalamus_chr6[grep(same_Hypothalamus_chr6[i],eqtl_Hypothalamus_chr6$loci),]
    chr6_Hypothalamus_genes<-c(chr6_Hypothalamus_genes,a)
    chr6_Hypothalamus_overloci<-rbind(chr6_Hypothalamus_overloci,b)
  }
}

chr7_Hypothalamus_genes<-NULL
chr7_Hypothalamus_overloci<-NULL
if(length(same_Hypothalamus_chr7!=0)){
  for(i in 1:length(same_Hypothalamus_chr7)){
    a<-eqtl_Hypothalamus_chr7$phenotype_id[grep(same_Hypothalamus_chr7[i],eqtl_Hypothalamus_chr7$loci)]
    b<-eqtl_Hypothalamus_chr7[grep(same_Hypothalamus_chr7[i],eqtl_Hypothalamus_chr7$loci),]
    chr7_Hypothalamus_genes<-c(chr7_Hypothalamus_genes,a)
    chr7_Hypothalamus_overloci<-rbind(chr7_Hypothalamus_overloci,b)
  }
}

chr8_Hypothalamus_genes<-NULL
chr8_Hypothalamus_overloci<-NULL
if(length(same_Hypothalamus_chr8!=0)){
  for(i in 1:length(same_Hypothalamus_chr8)){
    a<-eqtl_Hypothalamus_chr8$phenotype_id[grep(same_Hypothalamus_chr8[i],eqtl_Hypothalamus_chr8$loci)]
    b<-eqtl_Hypothalamus_chr8[grep(same_Hypothalamus_chr8[i],eqtl_Hypothalamus_chr8$loci),]
    chr8_Hypothalamus_genes<-c(chr8_Hypothalamus_genes,a)
    chr8_Hypothalamus_overloci<-rbind(chr8_Hypothalamus_overloci,b)
  }
}

if(length(same_Hypothalamus_chr9!=0)){
  chr9_Hypothalamus_genes<-NULL
  chr9_Hypothalamus_overloci<-NULL
  for(i in 1:length(same_Hypothalamus_chr9)){
    a<-eqtl_Hypothalamus_chr9$phenotype_id[grep(same_Hypothalamus_chr9[i],eqtl_Hypothalamus_chr9$loci)]
    b<-eqtl_Hypothalamus_chr9[grep(same_Hypothalamus_chr9[i],eqtl_Hypothalamus_chr9$loci),]
    chr9_Hypothalamus_genes<-c(chr9_Hypothalamus_genes,a)
    chr9_Hypothalamus_overloci<-rbind(chr9_Hypothalamus_overloci,b)
  }
}

chr10_Hypothalamus_genes<-NULL
chr10_Hypothalamus_overloci<-NULL
if(length(same_Hypothalamus_chr10!=0)){
  for(i in 1:length(same_Hypothalamus_chr10)){
    a<-eqtl_Hypothalamus_chr10$phenotype_id[grep(same_Hypothalamus_chr10[i],eqtl_Hypothalamus_chr10$loci)]
    b<-eqtl_Hypothalamus_chr10[grep(same_Hypothalamus_chr10[i],eqtl_Hypothalamus_chr10$loci),]
    chr10_Hypothalamus_genes<-c(chr10_Hypothalamus_genes,a)
    chr10_Hypothalamus_overloci<-rbind(chr10_Hypothalamus_overloci,b)
  }
}

chr11_Hypothalamus_genes<-NULL
chr11_Hypothalamus_overloci<-NULL
if(length(same_Hypothalamus_chr11!=0)){
  for(i in 1:length(same_Hypothalamus_chr11)){
    a<-eqtl_Hypothalamus_chr11$phenotype_id[grep(same_Hypothalamus_chr11[i],eqtl_Hypothalamus_chr11$loci)]
    b<-eqtl_Hypothalamus_chr11[grep(same_Hypothalamus_chr11[i],eqtl_Hypothalamus_chr11$loci),]
    chr11_Hypothalamus_genes<-c(chr11_Hypothalamus_genes,a)
    chr11_Hypothalamus_overloci<-rbind(chr11_Hypothalamus_overloci,b)
  }
}

chr12_Hypothalamus_genes<-NULL
chr12_Hypothalamus_overloci<-NULL
if(length(same_Hypothalamus_chr12!=0)){
  for(i in 1:length(same_Hypothalamus_chr12)){
    a<-eqtl_Hypothalamus_chr12$phenotype_id[grep(same_Hypothalamus_chr12[i],eqtl_Hypothalamus_chr12$loci)]
    b<-eqtl_Hypothalamus_chr12[grep(same_Hypothalamus_chr12[i],eqtl_Hypothalamus_chr12$loci),]
    chr12_Hypothalamus_genes<-c(chr12_Hypothalamus_genes,a)
    chr12_Hypothalamus_overloci<-rbind(chr12_Hypothalamus_overloci,b)
  }
}

chr13_Hypothalamus_genes<-NULL
chr13_Hypothalamus_overloci<-NULL
if(length(same_Hypothalamus_chr13!=0)){
  for(i in 1:length(same_Hypothalamus_chr13)){
    a<-eqtl_Hypothalamus_chr13$phenotype_id[grep(same_Hypothalamus_chr13[i],eqtl_Hypothalamus_chr13$loci)]
    b<-eqtl_Hypothalamus_chr13[grep(same_Hypothalamus_chr13[i],eqtl_Hypothalamus_chr13$loci),]
    chr13_Hypothalamus_genes<-c(chr13_Hypothalamus_genes,a)
    chr13_Hypothalamus_overloci<-rbind(chr13_Hypothalamus_overloci,b)
  }
}

chr14_Hypothalamus_genes<-NULL
chr14_Hypothalamus_overloci<-NULL
if(length(same_Hypothalamus_chr14!=0)){
  for(i in 1:length(same_Hypothalamus_chr14)){
    a<-eqtl_Hypothalamus_chr14$phenotype_id[grep(same_Hypothalamus_chr14[i],eqtl_Hypothalamus_chr14$loci)]
    b<-eqtl_Hypothalamus_chr14[grep(same_Hypothalamus_chr14[i],eqtl_Hypothalamus_chr14$loci),]
    chr14_Hypothalamus_genes<-c(chr14_Hypothalamus_genes,a)
    chr14_Hypothalamus_overloci<-rbind(chr14_Hypothalamus_overloci,b)
  }
}

chr15_Hypothalamus_genes<-NULL
chr15_Hypothalamus_overloci<-NULL
if(length(same_Hypothalamus_chr15!=0)){
  for(i in 1:length(same_Hypothalamus_chr15)){
    a<-eqtl_Hypothalamus_chr15$phenotype_id[grep(same_Hypothalamus_chr15[i],eqtl_Hypothalamus_chr15$loci)]
    b<-eqtl_Hypothalamus_chr15[grep(same_Hypothalamus_chr15[i],eqtl_Hypothalamus_chr15$loci),]
    chr15_Hypothalamus_genes<-c(chr15_Hypothalamus_genes,a)
    chr15_Hypothalamus_overloci<-rbind(chr15_Hypothalamus_overloci,b)
  }
}

chr16_Hypothalamus_genes<-NULL
chr16_Hypothalamus_overloci<-NULL
if(length(same_Hypothalamus_chr16!=0)){
  for(i in 1:length(same_Hypothalamus_chr16)){
    a<-eqtl_Hypothalamus_chr16$phenotype_id[grep(same_Hypothalamus_chr16[i],eqtl_Hypothalamus_chr16$loci)]
    b<-eqtl_Hypothalamus_chr16[grep(same_Hypothalamus_chr16[i],eqtl_Hypothalamus_chr16$loci),]
    chr16_Hypothalamus_genes<-c(chr16_Hypothalamus_genes,a)
    chr16_Hypothalamus_overloci<-rbind(chr16_Hypothalamus_overloci,b)
  }
}

chr17_Hypothalamus_genes<-NULL
chr17_Hypothalamus_overloci<-NULL
if(length(same_Hypothalamus_chr17!=0)){
  for(i in 1:length(same_Hypothalamus_chr17)){
    a<-eqtl_Hypothalamus_chr17$phenotype_id[grep(same_Hypothalamus_chr17[i],eqtl_Hypothalamus_chr17$loci)]
    b<-eqtl_Hypothalamus_chr17[grep(same_Hypothalamus_chr17[i],eqtl_Hypothalamus_chr17$loci),]
    chr17_Hypothalamus_genes<-c(chr17_Hypothalamus_genes,a)
    chr17_Hypothalamus_overloci<-rbind(chr17_Hypothalamus_overloci,b)
  }
}

chr18_Hypothalamus_genes<-NULL
chr18_Hypothalamus_overloci<-NULL
if(length(same_Hypothalamus_chr18!=0)){
  for(i in 1:length(same_Hypothalamus_chr18)){
    a<-eqtl_Hypothalamus_chr18$phenotype_id[grep(same_Hypothalamus_chr18[i],eqtl_Hypothalamus_chr18$loci)]
    b<-eqtl_Hypothalamus_chr18[grep(same_Hypothalamus_chr18[i],eqtl_Hypothalamus_chr18$loci),]
    chr18_Hypothalamus_genes<-c(chr18_Hypothalamus_genes,a)
    chr18_Hypothalamus_overloci<-rbind(chr18_Hypothalamus_overloci,b)
  }
}
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
one2one_pig<-annotation$`Pig gene stable ID`[which(annotation$`Pig homology type`=="ortholog_one2one")]

chr1_Hypothalamus_overloci_one2one<-chr1_Hypothalamus_overloci[match(intersect(chr1_Hypothalamus_genes,one2one_pig),chr1_Hypothalamus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Hypothalamus_chr1$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Hypothalamus_chr1$hum_chr<-tmp$V1
Hypothalamus_chr1$hum_loci<-tmp$V3
FM_Hypothalamus_hum$index<-paste0(FM_Hypothalamus_hum$chr,"_",FM_Hypothalamus_hum$variant_pos)
Hypothalamus_chr1$index<-paste0("chr",Hypothalamus_chr1$hum_chr,"_",Hypothalamus_chr1$hum_loci)
chr1_Hypothalamus_hum<-NULL
for(i in 1:length(same_Hypothalamus_chr1)){
  a<-Hypothalamus_chr1[match(same_Hypothalamus_chr1[i],Hypothalamus_chr1$pig_pos),]
  b<-FM_Hypothalamus_hum[match(intersect(a$index,FM_Hypothalamus_hum$index),FM_Hypothalamus_hum$index),]
  chr1_Hypothalamus_hum<-rbind(chr1_Hypothalamus_hum,b)
}
chr1_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr1_Hypothalamus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr1_Hypothalamus_hum_one2one<-chr1_Hypothalamus_hum[match(intersect(chr1_pig2hum_one2one,chr1_Hypothalamus_hum$gene_id),chr1_Hypothalamus_hum$gene_id),]

chr2_Hypothalamus_overloci_one2one<-chr2_Hypothalamus_overloci[match(intersect(chr2_Hypothalamus_genes,one2one_pig),chr2_Hypothalamus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Hypothalamus_chr2$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Hypothalamus_chr2$hum_chr<-tmp$V1
Hypothalamus_chr2$hum_loci<-tmp$V3
FM_Hypothalamus_hum$index<-paste0(FM_Hypothalamus_hum$chr,"_",FM_Hypothalamus_hum$variant_pos)
Hypothalamus_chr2$index<-paste0("chr",Hypothalamus_chr2$hum_chr,"_",Hypothalamus_chr2$hum_loci)
chr2_Hypothalamus_hum<-NULL
for(i in 1:length(same_Hypothalamus_chr2)){
  a<-Hypothalamus_chr2[match(same_Hypothalamus_chr2[i],Hypothalamus_chr2$pig_pos),]
  b<-FM_Hypothalamus_hum[match(intersect(a$index,FM_Hypothalamus_hum$index),FM_Hypothalamus_hum$index),]
  chr2_Hypothalamus_hum<-rbind(chr2_Hypothalamus_hum,b)
}
chr2_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr2_Hypothalamus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr2_Hypothalamus_hum_one2one<-chr2_Hypothalamus_hum[match(intersect(chr2_pig2hum_one2one,chr2_Hypothalamus_hum$gene_id),chr2_Hypothalamus_hum$gene_id),]

chr3_Hypothalamus_overloci_one2one<-chr3_Hypothalamus_overloci[match(intersect(chr3_Hypothalamus_genes,one2one_pig),chr3_Hypothalamus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Hypothalamus_chr3$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Hypothalamus_chr3$hum_chr<-tmp$V1
Hypothalamus_chr3$hum_loci<-tmp$V3
FM_Hypothalamus_hum$index<-paste0(FM_Hypothalamus_hum$chr,"_",FM_Hypothalamus_hum$variant_pos)
Hypothalamus_chr3$index<-paste0("chr",Hypothalamus_chr3$hum_chr,"_",Hypothalamus_chr3$hum_loci)
chr3_Hypothalamus_hum<-NULL
for(i in 1:length(same_Hypothalamus_chr3)){
  a<-Hypothalamus_chr3[match(same_Hypothalamus_chr3[i],Hypothalamus_chr3$pig_pos),]
  b<-FM_Hypothalamus_hum[match(intersect(a$index,FM_Hypothalamus_hum$index),FM_Hypothalamus_hum$index),]
  chr3_Hypothalamus_hum<-rbind(chr3_Hypothalamus_hum,b)
}
chr3_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr3_Hypothalamus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr3_Hypothalamus_hum_one2one<-chr3_Hypothalamus_hum[match(intersect(chr3_pig2hum_one2one,chr3_Hypothalamus_hum$gene_id),chr3_Hypothalamus_hum$gene_id),]

chr4_Hypothalamus_overloci_one2one<-chr4_Hypothalamus_overloci[match(intersect(chr4_Hypothalamus_genes,one2one_pig),chr4_Hypothalamus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Hypothalamus_chr4$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Hypothalamus_chr4$hum_chr<-tmp$V1
Hypothalamus_chr4$hum_loci<-tmp$V3
FM_Hypothalamus_hum$index<-paste0(FM_Hypothalamus_hum$chr,"_",FM_Hypothalamus_hum$variant_pos)
Hypothalamus_chr4$index<-paste0("chr",Hypothalamus_chr4$hum_chr,"_",Hypothalamus_chr4$hum_loci)
chr4_Hypothalamus_hum<-NULL
for(i in 1:length(same_Hypothalamus_chr4)){
  a<-Hypothalamus_chr4[match(same_Hypothalamus_chr4[i],Hypothalamus_chr4$pig_pos),]
  b<-FM_Hypothalamus_hum[match(intersect(a$index,FM_Hypothalamus_hum$index),FM_Hypothalamus_hum$index),]
  chr4_Hypothalamus_hum<-rbind(chr4_Hypothalamus_hum,b)
}
chr4_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr4_Hypothalamus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr4_Hypothalamus_hum_one2one<-chr4_Hypothalamus_hum[match(intersect(chr4_pig2hum_one2one,chr4_Hypothalamus_hum$gene_id),chr4_Hypothalamus_hum$gene_id),]

chr5_Hypothalamus_overloci_one2one<-chr5_Hypothalamus_overloci[match(intersect(chr5_Hypothalamus_genes,one2one_pig),chr5_Hypothalamus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Hypothalamus_chr5$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Hypothalamus_chr5$hum_chr<-tmp$V1
Hypothalamus_chr5$hum_loci<-tmp$V3
FM_Hypothalamus_hum$index<-paste0(FM_Hypothalamus_hum$chr,"_",FM_Hypothalamus_hum$variant_pos)
Hypothalamus_chr5$index<-paste0("chr",Hypothalamus_chr5$hum_chr,"_",Hypothalamus_chr5$hum_loci)
chr5_Hypothalamus_hum<-NULL
for(i in 1:length(same_Hypothalamus_chr5)){
  a<-Hypothalamus_chr5[match(same_Hypothalamus_chr5[i],Hypothalamus_chr5$pig_pos),]
  b<-FM_Hypothalamus_hum[match(intersect(a$index,FM_Hypothalamus_hum$index),FM_Hypothalamus_hum$index),]
  chr5_Hypothalamus_hum<-rbind(chr5_Hypothalamus_hum,b)
}
chr5_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr5_Hypothalamus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr5_Hypothalamus_hum_one2one<-chr5_Hypothalamus_hum[match(intersect(chr5_pig2hum_one2one,chr5_Hypothalamus_hum$gene_id),chr5_Hypothalamus_hum$gene_id),]

chr6_Hypothalamus_overloci_one2one<-chr6_Hypothalamus_overloci[match(intersect(chr6_Hypothalamus_genes,one2one_pig),chr6_Hypothalamus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Hypothalamus_chr6$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Hypothalamus_chr6$hum_chr<-tmp$V1
Hypothalamus_chr6$hum_loci<-tmp$V3
FM_Hypothalamus_hum$index<-paste0(FM_Hypothalamus_hum$chr,"_",FM_Hypothalamus_hum$variant_pos)
Hypothalamus_chr6$index<-paste0("chr",Hypothalamus_chr6$hum_chr,"_",Hypothalamus_chr6$hum_loci)
chr6_Hypothalamus_hum<-NULL
for(i in 1:length(same_Hypothalamus_chr6)){
  a<-Hypothalamus_chr6[match(same_Hypothalamus_chr6[i],Hypothalamus_chr6$pig_pos),]
  b<-FM_Hypothalamus_hum[match(intersect(a$index,FM_Hypothalamus_hum$index),FM_Hypothalamus_hum$index),]
  chr6_Hypothalamus_hum<-rbind(chr6_Hypothalamus_hum,b)
}
chr6_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr6_Hypothalamus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr6_Hypothalamus_hum_one2one<-chr6_Hypothalamus_hum[match(intersect(chr6_pig2hum_one2one,chr6_Hypothalamus_hum$gene_id),chr6_Hypothalamus_hum$gene_id),]

chr7_Hypothalamus_overloci_one2one<-chr7_Hypothalamus_overloci[match(intersect(chr7_Hypothalamus_genes,one2one_pig),chr7_Hypothalamus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Hypothalamus_chr7$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Hypothalamus_chr7$hum_chr<-tmp$V1
Hypothalamus_chr7$hum_loci<-tmp$V3
FM_Hypothalamus_hum$index<-paste0(FM_Hypothalamus_hum$chr,"_",FM_Hypothalamus_hum$variant_pos)
Hypothalamus_chr7$index<-paste0("chr",Hypothalamus_chr7$hum_chr,"_",Hypothalamus_chr7$hum_loci)
chr7_Hypothalamus_hum<-NULL
for(i in 1:length(same_Hypothalamus_chr7)){
  a<-Hypothalamus_chr7[match(same_Hypothalamus_chr7[i],Hypothalamus_chr7$pig_pos),]
  b<-FM_Hypothalamus_hum[match(intersect(a$index,FM_Hypothalamus_hum$index),FM_Hypothalamus_hum$index),]
  chr7_Hypothalamus_hum<-rbind(chr7_Hypothalamus_hum,b)
}
chr7_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr7_Hypothalamus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr7_Hypothalamus_hum_one2one<-chr7_Hypothalamus_hum[match(intersect(chr7_pig2hum_one2one,chr7_Hypothalamus_hum$gene_id),chr7_Hypothalamus_hum$gene_id),]

chr8_Hypothalamus_overloci_one2one<-chr8_Hypothalamus_overloci[match(intersect(chr8_Hypothalamus_genes,one2one_pig),chr8_Hypothalamus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Hypothalamus_chr8$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Hypothalamus_chr8$hum_chr<-tmp$V1
Hypothalamus_chr8$hum_loci<-tmp$V3
FM_Hypothalamus_hum$index<-paste0(FM_Hypothalamus_hum$chr,"_",FM_Hypothalamus_hum$variant_pos)
Hypothalamus_chr8$index<-paste0("chr",Hypothalamus_chr8$hum_chr,"_",Hypothalamus_chr8$hum_loci)
chr8_Hypothalamus_hum<-NULL
for(i in 1:length(same_Hypothalamus_chr8)){
  a<-Hypothalamus_chr8[match(same_Hypothalamus_chr8[i],Hypothalamus_chr8$pig_pos),]
  b<-FM_Hypothalamus_hum[match(intersect(a$index,FM_Hypothalamus_hum$index),FM_Hypothalamus_hum$index),]
  chr8_Hypothalamus_hum<-rbind(chr8_Hypothalamus_hum,b)
}
chr8_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr8_Hypothalamus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr8_Hypothalamus_hum_one2one<-chr8_Hypothalamus_hum[match(intersect(chr8_pig2hum_one2one,chr8_Hypothalamus_hum$gene_id),chr8_Hypothalamus_hum$gene_id),]

chr9_Hypothalamus_overloci_one2one<-chr9_Hypothalamus_overloci[match(intersect(chr9_Hypothalamus_genes,one2one_pig),chr9_Hypothalamus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Hypothalamus_chr9$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Hypothalamus_chr9$hum_chr<-tmp$V1
Hypothalamus_chr9$hum_loci<-tmp$V3
FM_Hypothalamus_hum$index<-paste0(FM_Hypothalamus_hum$chr,"_",FM_Hypothalamus_hum$variant_pos)
Hypothalamus_chr9$index<-paste0("chr",Hypothalamus_chr9$hum_chr,"_",Hypothalamus_chr9$hum_loci)
chr9_Hypothalamus_hum<-NULL
for(i in 1:length(same_Hypothalamus_chr9)){
  a<-Hypothalamus_chr9[match(same_Hypothalamus_chr9[i],Hypothalamus_chr9$pig_pos),]
  b<-FM_Hypothalamus_hum[match(intersect(a$index,FM_Hypothalamus_hum$index),FM_Hypothalamus_hum$index),]
  chr9_Hypothalamus_hum<-rbind(chr9_Hypothalamus_hum,b)
}
chr9_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr9_Hypothalamus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr9_Hypothalamus_hum_one2one<-chr9_Hypothalamus_hum[match(intersect(chr9_pig2hum_one2one,chr9_Hypothalamus_hum$gene_id),chr9_Hypothalamus_hum$gene_id),]

chr10_Hypothalamus_overloci_one2one<-chr10_Hypothalamus_overloci[match(intersect(chr10_Hypothalamus_genes,one2one_pig),chr10_Hypothalamus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Hypothalamus_chr10$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Hypothalamus_chr10$hum_chr<-tmp$V1
Hypothalamus_chr10$hum_loci<-tmp$V3
FM_Hypothalamus_hum$index<-paste0(FM_Hypothalamus_hum$chr,"_",FM_Hypothalamus_hum$variant_pos)
Hypothalamus_chr10$index<-paste0("chr",Hypothalamus_chr10$hum_chr,"_",Hypothalamus_chr10$hum_loci)
chr10_Hypothalamus_hum<-NULL
for(i in 1:length(same_Hypothalamus_chr10)){
  a<-Hypothalamus_chr10[match(same_Hypothalamus_chr10[i],Hypothalamus_chr10$pig_pos),]
  b<-FM_Hypothalamus_hum[match(intersect(a$index,FM_Hypothalamus_hum$index),FM_Hypothalamus_hum$index),]
  chr10_Hypothalamus_hum<-rbind(chr10_Hypothalamus_hum,b)
}
chr10_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr10_Hypothalamus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr10_Hypothalamus_hum_one2one<-chr10_Hypothalamus_hum[match(intersect(chr10_pig2hum_one2one,chr10_Hypothalamus_hum$gene_id),chr10_Hypothalamus_hum$gene_id),]

chr11_Hypothalamus_overloci_one2one<-chr11_Hypothalamus_overloci[match(intersect(chr11_Hypothalamus_genes,one2one_pig),chr11_Hypothalamus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Hypothalamus_chr11$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Hypothalamus_chr11$hum_chr<-tmp$V1
Hypothalamus_chr11$hum_loci<-tmp$V3
FM_Hypothalamus_hum$index<-paste0(FM_Hypothalamus_hum$chr,"_",FM_Hypothalamus_hum$variant_pos)
Hypothalamus_chr11$index<-paste0("chr",Hypothalamus_chr11$hum_chr,"_",Hypothalamus_chr11$hum_loci)
chr11_Hypothalamus_hum<-NULL
for(i in 1:length(same_Hypothalamus_chr11)){
  a<-Hypothalamus_chr11[match(same_Hypothalamus_chr11[i],Hypothalamus_chr11$pig_pos),]
  b<-FM_Hypothalamus_hum[match(intersect(a$index,FM_Hypothalamus_hum$index),FM_Hypothalamus_hum$index),]
  chr11_Hypothalamus_hum<-rbind(chr11_Hypothalamus_hum,b)
}
chr11_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr11_Hypothalamus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr11_Hypothalamus_hum_one2one<-chr11_Hypothalamus_hum[match(intersect(chr11_pig2hum_one2one,chr11_Hypothalamus_hum$gene_id),chr11_Hypothalamus_hum$gene_id),]

chr12_Hypothalamus_overloci_one2one<-chr12_Hypothalamus_overloci[match(intersect(chr12_Hypothalamus_genes,one2one_pig),chr12_Hypothalamus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Hypothalamus_chr12$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Hypothalamus_chr12$hum_chr<-tmp$V1
Hypothalamus_chr12$hum_loci<-tmp$V3
FM_Hypothalamus_hum$index<-paste0(FM_Hypothalamus_hum$chr,"_",FM_Hypothalamus_hum$variant_pos)
Hypothalamus_chr12$index<-paste0("chr",Hypothalamus_chr12$hum_chr,"_",Hypothalamus_chr12$hum_loci)
chr12_Hypothalamus_hum<-NULL
for(i in 1:length(same_Hypothalamus_chr12)){
  a<-Hypothalamus_chr12[match(same_Hypothalamus_chr12[i],Hypothalamus_chr12$pig_pos),]
  b<-FM_Hypothalamus_hum[match(intersect(a$index,FM_Hypothalamus_hum$index),FM_Hypothalamus_hum$index),]
  chr12_Hypothalamus_hum<-rbind(chr12_Hypothalamus_hum,b)
}
chr12_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr12_Hypothalamus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr12_Hypothalamus_hum_one2one<-chr12_Hypothalamus_hum[match(intersect(chr12_pig2hum_one2one,chr12_Hypothalamus_hum$gene_id),chr12_Hypothalamus_hum$gene_id),]

chr13_Hypothalamus_overloci_one2one<-chr13_Hypothalamus_overloci[match(intersect(chr13_Hypothalamus_genes,one2one_pig),chr13_Hypothalamus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Hypothalamus_chr13$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Hypothalamus_chr13$hum_chr<-tmp$V1
Hypothalamus_chr13$hum_loci<-tmp$V3
FM_Hypothalamus_hum$index<-paste0(FM_Hypothalamus_hum$chr,"_",FM_Hypothalamus_hum$variant_pos)
Hypothalamus_chr13$index<-paste0("chr",Hypothalamus_chr13$hum_chr,"_",Hypothalamus_chr13$hum_loci)
chr13_Hypothalamus_hum<-NULL
for(i in 1:length(same_Hypothalamus_chr13)){
  a<-Hypothalamus_chr13[match(same_Hypothalamus_chr13[i],Hypothalamus_chr13$pig_pos),]
  b<-FM_Hypothalamus_hum[match(intersect(a$index,FM_Hypothalamus_hum$index),FM_Hypothalamus_hum$index),]
  chr13_Hypothalamus_hum<-rbind(chr13_Hypothalamus_hum,b)
}
chr13_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr13_Hypothalamus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr13_Hypothalamus_hum_one2one<-chr13_Hypothalamus_hum[match(intersect(chr13_pig2hum_one2one,chr13_Hypothalamus_hum$gene_id),chr13_Hypothalamus_hum$gene_id),]

chr14_Hypothalamus_overloci_one2one<-chr14_Hypothalamus_overloci[match(intersect(chr14_Hypothalamus_genes,one2one_pig),chr14_Hypothalamus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Hypothalamus_chr14$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Hypothalamus_chr14$hum_chr<-tmp$V1
Hypothalamus_chr14$hum_loci<-tmp$V3
FM_Hypothalamus_hum$index<-paste0(FM_Hypothalamus_hum$chr,"_",FM_Hypothalamus_hum$variant_pos)
Hypothalamus_chr14$index<-paste0("chr",Hypothalamus_chr14$hum_chr,"_",Hypothalamus_chr14$hum_loci)
chr14_Hypothalamus_hum<-NULL
for(i in 1:length(same_Hypothalamus_chr14)){
  a<-Hypothalamus_chr14[match(same_Hypothalamus_chr14[i],Hypothalamus_chr14$pig_pos),]
  b<-FM_Hypothalamus_hum[match(intersect(a$index,FM_Hypothalamus_hum$index),FM_Hypothalamus_hum$index),]
  chr14_Hypothalamus_hum<-rbind(chr14_Hypothalamus_hum,b)
}
chr14_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr14_Hypothalamus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr14_Hypothalamus_hum_one2one<-chr14_Hypothalamus_hum[match(intersect(chr14_pig2hum_one2one,chr14_Hypothalamus_hum$gene_id),chr14_Hypothalamus_hum$gene_id),]

chr15_Hypothalamus_overloci_one2one<-chr15_Hypothalamus_overloci[match(intersect(chr15_Hypothalamus_genes,one2one_pig),chr15_Hypothalamus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Hypothalamus_chr15$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Hypothalamus_chr15$hum_chr<-tmp$V1
Hypothalamus_chr15$hum_loci<-tmp$V3
FM_Hypothalamus_hum$index<-paste0(FM_Hypothalamus_hum$chr,"_",FM_Hypothalamus_hum$variant_pos)
Hypothalamus_chr15$index<-paste0("chr",Hypothalamus_chr15$hum_chr,"_",Hypothalamus_chr15$hum_loci)
chr15_Hypothalamus_hum<-NULL
for(i in 1:length(same_Hypothalamus_chr15)){
  a<-Hypothalamus_chr15[match(same_Hypothalamus_chr15[i],Hypothalamus_chr15$pig_pos),]
  b<-FM_Hypothalamus_hum[match(intersect(a$index,FM_Hypothalamus_hum$index),FM_Hypothalamus_hum$index),]
  chr15_Hypothalamus_hum<-rbind(chr15_Hypothalamus_hum,b)
}
chr15_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr15_Hypothalamus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr15_Hypothalamus_hum_one2one<-chr15_Hypothalamus_hum[match(intersect(chr15_pig2hum_one2one,chr15_Hypothalamus_hum$gene_id),chr15_Hypothalamus_hum$gene_id),]

chr16_Hypothalamus_overloci_one2one<-chr16_Hypothalamus_overloci[match(intersect(chr16_Hypothalamus_genes,one2one_pig),chr16_Hypothalamus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Hypothalamus_chr16$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Hypothalamus_chr16$hum_chr<-tmp$V1
Hypothalamus_chr16$hum_loci<-tmp$V3
FM_Hypothalamus_hum$index<-paste0(FM_Hypothalamus_hum$chr,"_",FM_Hypothalamus_hum$variant_pos)
Hypothalamus_chr16$index<-paste0("chr",Hypothalamus_chr16$hum_chr,"_",Hypothalamus_chr16$hum_loci)
chr16_Hypothalamus_hum<-NULL
for(i in 1:length(same_Hypothalamus_chr16)){
  a<-Hypothalamus_chr16[match(same_Hypothalamus_chr16[i],Hypothalamus_chr16$pig_pos),]
  b<-FM_Hypothalamus_hum[match(intersect(a$index,FM_Hypothalamus_hum$index),FM_Hypothalamus_hum$index),]
  chr16_Hypothalamus_hum<-rbind(chr16_Hypothalamus_hum,b)
}
chr16_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr16_Hypothalamus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr16_Hypothalamus_hum_one2one<-chr16_Hypothalamus_hum[match(intersect(chr16_pig2hum_one2one,chr16_Hypothalamus_hum$gene_id),chr16_Hypothalamus_hum$gene_id),]

chr17_Hypothalamus_overloci_one2one<-chr17_Hypothalamus_overloci[match(intersect(chr17_Hypothalamus_genes,one2one_pig),chr17_Hypothalamus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Hypothalamus_chr17$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Hypothalamus_chr17$hum_chr<-tmp$V1
Hypothalamus_chr17$hum_loci<-tmp$V3
FM_Hypothalamus_hum$index<-paste0(FM_Hypothalamus_hum$chr,"_",FM_Hypothalamus_hum$variant_pos)
Hypothalamus_chr17$index<-paste0("chr",Hypothalamus_chr17$hum_chr,"_",Hypothalamus_chr17$hum_loci)
chr17_Hypothalamus_hum<-NULL
for(i in 1:length(same_Hypothalamus_chr17)){
  a<-Hypothalamus_chr17[match(same_Hypothalamus_chr17[i],Hypothalamus_chr17$pig_pos),]
  b<-FM_Hypothalamus_hum[match(intersect(a$index,FM_Hypothalamus_hum$index),FM_Hypothalamus_hum$index),]
  chr17_Hypothalamus_hum<-rbind(chr17_Hypothalamus_hum,b)
}
chr17_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr17_Hypothalamus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr17_Hypothalamus_hum_one2one<-chr17_Hypothalamus_hum[match(intersect(chr17_pig2hum_one2one,chr17_Hypothalamus_hum$gene_id),chr17_Hypothalamus_hum$gene_id),]

chr18_Hypothalamus_overloci_one2one<-chr18_Hypothalamus_overloci[match(intersect(chr18_Hypothalamus_genes,one2one_pig),chr18_Hypothalamus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Hypothalamus_chr18$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Hypothalamus_chr18$hum_chr<-tmp$V1
Hypothalamus_chr18$hum_loci<-tmp$V3
FM_Hypothalamus_hum$index<-paste0(FM_Hypothalamus_hum$chr,"_",FM_Hypothalamus_hum$variant_pos)
Hypothalamus_chr18$index<-paste0("chr",Hypothalamus_chr18$hum_chr,"_",Hypothalamus_chr18$hum_loci)
chr18_Hypothalamus_hum<-NULL
for(i in 1:length(same_Hypothalamus_chr18)){
  a<-Hypothalamus_chr18[match(same_Hypothalamus_chr18[i],Hypothalamus_chr18$pig_pos),]
  b<-FM_Hypothalamus_hum[match(intersect(a$index,FM_Hypothalamus_hum$index),FM_Hypothalamus_hum$index),]
  chr18_Hypothalamus_hum<-rbind(chr18_Hypothalamus_hum,b)
}
chr18_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr18_Hypothalamus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr18_Hypothalamus_hum_one2one<-chr18_Hypothalamus_hum[match(intersect(chr18_pig2hum_one2one,chr18_Hypothalamus_hum$gene_id),chr18_Hypothalamus_hum$gene_id),]

Hypothalamus_one2one_SNP_hum<-rbind(chr1_Hypothalamus_hum_one2one,chr2_Hypothalamus_hum_one2one,chr3_Hypothalamus_hum_one2one,chr4_Hypothalamus_hum_one2one,chr5_Hypothalamus_hum_one2one,
                                    chr6_Hypothalamus_hum_one2one,chr7_Hypothalamus_hum_one2one,chr8_Hypothalamus_hum_one2one,chr9_Hypothalamus_hum_one2one,chr10_Hypothalamus_hum_one2one,
                                    chr11_Hypothalamus_hum_one2one,chr12_Hypothalamus_hum_one2one,chr13_Hypothalamus_hum_one2one,chr14_Hypothalamus_hum_one2one,chr15_Hypothalamus_hum_one2one,
                                    chr16_Hypothalamus_hum_one2one,chr17_Hypothalamus_hum_one2one,chr18_Hypothalamus_hum_one2one)

chr1_Hypothalamus_pig_one2one<-chr1_Hypothalamus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr1_Hypothalamus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr1_Hypothalamus_overloci_one2one$phenotype_id),1:9]
chr2_Hypothalamus_pig_one2one<-chr2_Hypothalamus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr2_Hypothalamus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr2_Hypothalamus_overloci_one2one$phenotype_id),1:9]
chr3_Hypothalamus_pig_one2one<-chr3_Hypothalamus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr3_Hypothalamus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr3_Hypothalamus_overloci_one2one$phenotype_id),1:9]
chr4_Hypothalamus_pig_one2one<-chr4_Hypothalamus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr4_Hypothalamus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr4_Hypothalamus_overloci_one2one$phenotype_id),1:9]
chr5_Hypothalamus_pig_one2one<-chr5_Hypothalamus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr5_Hypothalamus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr5_Hypothalamus_overloci_one2one$phenotype_id),1:9]
chr6_Hypothalamus_pig_one2one<-chr6_Hypothalamus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr6_Hypothalamus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr6_Hypothalamus_overloci_one2one$phenotype_id),1:9]
chr7_Hypothalamus_pig_one2one<-chr7_Hypothalamus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr7_Hypothalamus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr7_Hypothalamus_overloci_one2one$phenotype_id),1:9]
chr8_Hypothalamus_pig_one2one<-chr8_Hypothalamus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr8_Hypothalamus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr8_Hypothalamus_overloci_one2one$phenotype_id),1:9]
chr9_Hypothalamus_pig_one2one<-chr9_Hypothalamus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr9_Hypothalamus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr9_Hypothalamus_overloci_one2one$phenotype_id),1:9]
chr10_Hypothalamus_pig_one2one<-chr10_Hypothalamus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr10_Hypothalamus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr10_Hypothalamus_overloci_one2one$phenotype_id),1:9]
chr11_Hypothalamus_pig_one2one<-chr11_Hypothalamus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr11_Hypothalamus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr11_Hypothalamus_overloci_one2one$phenotype_id),1:9]
chr12_Hypothalamus_pig_one2one<-chr12_Hypothalamus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr12_Hypothalamus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr12_Hypothalamus_overloci_one2one$phenotype_id),1:9]
chr13_Hypothalamus_pig_one2one<-chr13_Hypothalamus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr13_Hypothalamus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr13_Hypothalamus_overloci_one2one$phenotype_id),1:9]
chr14_Hypothalamus_pig_one2one<-chr14_Hypothalamus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr14_Hypothalamus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr14_Hypothalamus_overloci_one2one$phenotype_id),1:9]
chr15_Hypothalamus_pig_one2one<-chr15_Hypothalamus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr15_Hypothalamus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr15_Hypothalamus_overloci_one2one$phenotype_id),1:9]
chr16_Hypothalamus_pig_one2one<-chr16_Hypothalamus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr16_Hypothalamus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr16_Hypothalamus_overloci_one2one$phenotype_id),1:9]
chr17_Hypothalamus_pig_one2one<-chr17_Hypothalamus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr17_Hypothalamus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr17_Hypothalamus_overloci_one2one$phenotype_id),1:9]
chr18_Hypothalamus_pig_one2one<-chr18_Hypothalamus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr18_Hypothalamus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr18_Hypothalamus_overloci_one2one$phenotype_id),1:9]

Hypothalamus_one2one_SNP_pig<-rbind(chr1_Hypothalamus_pig_one2one,chr2_Hypothalamus_pig_one2one,chr3_Hypothalamus_pig_one2one,chr4_Hypothalamus_pig_one2one,chr5_Hypothalamus_pig_one2one,
                                    chr6_Hypothalamus_pig_one2one,chr7_Hypothalamus_pig_one2one,chr8_Hypothalamus_pig_one2one,chr9_Hypothalamus_pig_one2one,chr10_Hypothalamus_pig_one2one,
                                    chr11_Hypothalamus_pig_one2one,chr12_Hypothalamus_pig_one2one,chr13_Hypothalamus_pig_one2one,chr14_Hypothalamus_pig_one2one,chr15_Hypothalamus_pig_one2one,
                                    chr16_Hypothalamus_pig_one2one,chr17_Hypothalamus_pig_one2one,chr18_Hypothalamus_pig_one2one)

Hypothalamus_SNP_sum<-array(NA,dim=c(nrow(Hypothalamus_one2one_SNP_hum),2))
colnames(Hypothalamus_SNP_sum)<-c("Human","Pig")
Hypothalamus_SNP_sum<-as.data.frame(Hypothalamus_SNP_sum)
Hypothalamus_SNP_sum$Human<-Hypothalamus_one2one_SNP_hum$slope / Hypothalamus_one2one_SNP_hum$slope_se
Hypothalamus_SNP_sum$Pig<-Hypothalamus_one2one_SNP_pig$slope / Hypothalamus_one2one_SNP_pig$slope_se
cor<-cor(abs(Hypothalamus_SNP_sum$Human),abs(Hypothalamus_SNP_sum$Pig))
p_val<-t.test(abs(Hypothalamus_SNP_sum$Human),abs(Hypothalamus_SNP_sum$Pig))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/Hypothalamus_SNP_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(Hypothalamus_SNP_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The number of eGenes expressed in common tissues",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

save(Hypothalamus_one2one_SNP_hum,Hypothalamus_one2one_SNP_pig,file="/Users/baizhonghao/Downloads/human-pig GTEx/Hypothalamus_SNP.Rdata")

#Ileum_SNP_overlaploci#
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

eqtl_Ileum_chr1<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ileum/Ileum.cis_qtl_pairs.1.txt"))
eqtl_Ileum_chr2<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ileum/Ileum.cis_qtl_pairs.2.txt"))
eqtl_Ileum_chr3<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ileum/Ileum.cis_qtl_pairs.3.txt"))
eqtl_Ileum_chr4<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ileum/Ileum.cis_qtl_pairs.4.txt"))
eqtl_Ileum_chr5<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ileum/Ileum.cis_qtl_pairs.5.txt"))
eqtl_Ileum_chr6<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ileum/Ileum.cis_qtl_pairs.6.txt"))
eqtl_Ileum_chr7<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ileum/Ileum.cis_qtl_pairs.7.txt"))
eqtl_Ileum_chr8<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ileum/Ileum.cis_qtl_pairs.8.txt"))
eqtl_Ileum_chr9<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ileum/Ileum.cis_qtl_pairs.9.txt"))
eqtl_Ileum_chr10<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ileum/Ileum.cis_qtl_pairs.10.txt"))
eqtl_Ileum_chr11<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ileum/Ileum.cis_qtl_pairs.11.txt"))
eqtl_Ileum_chr12<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ileum/Ileum.cis_qtl_pairs.12.txt"))
eqtl_Ileum_chr13<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ileum/Ileum.cis_qtl_pairs.13.txt"))
eqtl_Ileum_chr14<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ileum/Ileum.cis_qtl_pairs.14.txt"))
eqtl_Ileum_chr15<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ileum/Ileum.cis_qtl_pairs.15.txt"))
eqtl_Ileum_chr16<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ileum/Ileum.cis_qtl_pairs.16.txt"))
eqtl_Ileum_chr17<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ileum/Ileum.cis_qtl_pairs.17.txt"))
eqtl_Ileum_chr18<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ileum/Ileum.cis_qtl_pairs.18.txt"))

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum_chr1$variant_id,split="_"))))
eqtl_Ileum_chr1$chr<-tmp$V1
eqtl_Ileum_chr1$loci<-tmp$V2

eqtl_Ileum_chr1$index<-paste0(eqtl_Ileum_chr1$chr,"-",eqtl_Ileum_chr1$loci)
Ileum_chr1$index<-paste0(Ileum_chr1$chr,"-",Ileum_chr1$pig_pos)
same_Ileum_chr1<-intersect(Ileum_chr1$pig_pos,eqtl_Ileum_chr1$loci)


tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum_chr2$variant_id,split="_"))))
eqtl_Ileum_chr2$chr<-tmp$V1
eqtl_Ileum_chr2$loci<-tmp$V2

eqtl_Ileum_chr2$index<-paste0(eqtl_Ileum_chr2$chr,"-",eqtl_Ileum_chr2$loci)
Ileum_chr2$index<-paste0(Ileum_chr2$chr,"-",Ileum_chr2$pig_pos)
same_Ileum_chr2<-intersect(Ileum_chr2$pig_pos,eqtl_Ileum_chr2$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum_chr3$variant_id,split="_"))))
eqtl_Ileum_chr3$chr<-tmp$V1
eqtl_Ileum_chr3$loci<-tmp$V2

eqtl_Ileum_chr3$index<-paste0(eqtl_Ileum_chr3$chr,"-",eqtl_Ileum_chr3$loci)
Ileum_chr3$index<-paste0(Ileum_chr3$chr,"-",Ileum_chr3$pig_pos)
same_Ileum_chr3<-intersect(Ileum_chr3$pig_pos,eqtl_Ileum_chr3$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum_chr4$variant_id,split="_"))))
eqtl_Ileum_chr4$chr<-tmp$V1
eqtl_Ileum_chr4$loci<-tmp$V2

eqtl_Ileum_chr4$index<-paste0(eqtl_Ileum_chr4$chr,"-",eqtl_Ileum_chr4$loci)
Ileum_chr4$index<-paste0(Ileum_chr4$chr,"-",Ileum_chr4$pig_pos)
same_Ileum_chr4<-intersect(Ileum_chr4$pig_pos,eqtl_Ileum_chr4$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum_chr5$variant_id,split="_"))))
eqtl_Ileum_chr5$loci<-tmp$V2
same_Ileum_chr5<-intersect(Ileum_chr5$pig_pos,eqtl_Ileum_chr5$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum_chr6$variant_id,split="_"))))
eqtl_Ileum_chr6$loci<-tmp$V2
same_Ileum_chr6<-intersect(Ileum_chr6$pig_pos,eqtl_Ileum_chr6$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum_chr7$variant_id,split="_"))))
eqtl_Ileum_chr7$loci<-tmp$V2
same_Ileum_chr7<-intersect(Ileum_chr7$pig_pos,eqtl_Ileum_chr7$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum_chr8$variant_id,split="_"))))
eqtl_Ileum_chr8$loci<-tmp$V2
same_Ileum_chr8<-intersect(Ileum_chr8$pig_pos,eqtl_Ileum_chr8$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum_chr9$variant_id,split="_"))))
eqtl_Ileum_chr9$loci<-tmp$V2
same_Ileum_chr9<-intersect(Ileum_chr9$pig_pos,eqtl_Ileum_chr9$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum_chr10$variant_id,split="_"))))
eqtl_Ileum_chr10$loci<-tmp$V2
same_Ileum_chr10<-intersect(Ileum_chr10$pig_pos,eqtl_Ileum_chr10$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum_chr11$variant_id,split="_"))))
eqtl_Ileum_chr11$loci<-tmp$V2
same_Ileum_chr11<-intersect(Ileum_chr11$pig_pos,eqtl_Ileum_chr11$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum_chr12$variant_id,split="_"))))
eqtl_Ileum_chr12$loci<-tmp$V2
same_Ileum_chr12<-intersect(Ileum_chr12$pig_pos,eqtl_Ileum_chr12$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum_chr13$variant_id,split="_"))))
eqtl_Ileum_chr13$loci<-tmp$V2
same_Ileum_chr13<-intersect(Ileum_chr13$pig_pos,eqtl_Ileum_chr13$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum_chr14$variant_id,split="_"))))
eqtl_Ileum_chr14$loci<-tmp$V2
same_Ileum_chr14<-intersect(Ileum_chr14$pig_pos,eqtl_Ileum_chr14$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum_chr15$variant_id,split="_"))))
eqtl_Ileum_chr15$loci<-tmp$V2
same_Ileum_chr15<-intersect(Ileum_chr15$pig_pos,eqtl_Ileum_chr15$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum_chr16$variant_id,split="_"))))
eqtl_Ileum_chr16$loci<-tmp$V2
same_Ileum_chr16<-intersect(Ileum_chr16$pig_pos,eqtl_Ileum_chr16$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum_chr17$variant_id,split="_"))))
eqtl_Ileum_chr17$loci<-tmp$V2
same_Ileum_chr17<-intersect(Ileum_chr17$pig_pos,eqtl_Ileum_chr17$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ileum_chr18$variant_id,split="_"))))
eqtl_Ileum_chr18$loci<-tmp$V2
same_Ileum_chr18<-intersect(Ileum_chr18$pig_pos,eqtl_Ileum_chr18$loci)

chr1_Ileum_genes<-NULL
chr1_Ileum_overloci<-NULL
if(length(same_Ileum_chr1!=0)){
  for(i in 1:length(same_Ileum_chr1)){
    a<-eqtl_Ileum_chr1$phenotype_id[grep(same_Ileum_chr1[i],eqtl_Ileum_chr1$loci)]
    b<-eqtl_Ileum_chr1[grep(same_Ileum_chr1[i],eqtl_Ileum_chr1$loci),]
    chr1_Ileum_genes<-c(chr1_Ileum_genes,a)
    chr1_Ileum_overloci<-rbind(chr1_Ileum_overloci,b)
  }
}

chr2_Ileum_genes<-NULL
chr2_Ileum_overloci<-NULL
if(length(same_Ileum_chr2!=0)){
  for(i in 1:length(same_Ileum_chr2)){
    a<-eqtl_Ileum_chr2$phenotype_id[grep(same_Ileum_chr2[i],eqtl_Ileum_chr2$loci)]
    b<-eqtl_Ileum_chr2[grep(same_Ileum_chr2[i],eqtl_Ileum_chr2$loci),]
    chr2_Ileum_genes<-c(chr2_Ileum_genes,a)
    chr2_Ileum_overloci<-rbind(chr2_Ileum_overloci,b)
  }
}

chr3_Ileum_genes<-NULL
chr3_Ileum_overloci<-NULL
if(length(same_Ileum_chr3!=0)){
  for(i in 1:length(same_Ileum_chr3)){
    a<-eqtl_Ileum_chr3$phenotype_id[grep(same_Ileum_chr3[i],eqtl_Ileum_chr3$loci)]
    b<-eqtl_Ileum_chr3[grep(same_Ileum_chr3[i],eqtl_Ileum_chr3$loci),]
    chr3_Ileum_genes<-c(chr3_Ileum_genes,a)
    chr3_Ileum_overloci<-rbind(chr3_Ileum_overloci,b)
  }
}

chr4_Ileum_genes<-NULL
chr4_Ileum_overloci<-NULL
if(length(same_Ileum_chr4!=0)){
  for(i in 1:length(same_Ileum_chr4)){
    a<-eqtl_Ileum_chr4$phenotype_id[grep(same_Ileum_chr4[i],eqtl_Ileum_chr4$loci)]
    b<-eqtl_Ileum_chr4[grep(same_Ileum_chr4[i],eqtl_Ileum_chr4$loci),]
    chr4_Ileum_genes<-c(chr4_Ileum_genes,a)
    chr4_Ileum_overloci<-rbind(chr4_Ileum_overloci,b)
  }
}

chr5_Ileum_genes<-NULL
chr5_Ileum_overloci<-NULL
if(length(same_Ileum_chr5!=0)){
  for(i in 1:length(same_Ileum_chr5)){
    a<-eqtl_Ileum_chr5$phenotype_id[grep(same_Ileum_chr5[i],eqtl_Ileum_chr5$loci)]
    b<-eqtl_Ileum_chr5[grep(same_Ileum_chr5[i],eqtl_Ileum_chr5$loci),]
    chr5_Ileum_genes<-c(chr5_Ileum_genes,a)
    chr5_Ileum_overloci<-rbind(chr5_Ileum_overloci,b)
  }
}

chr6_Ileum_genes<-NULL
chr6_Ileum_overloci<-NULL
if(length(same_Ileum_chr6!=0)){
  for(i in 1:length(same_Ileum_chr6)){
    a<-eqtl_Ileum_chr6$phenotype_id[grep(same_Ileum_chr6[i],eqtl_Ileum_chr6$loci)]
    b<-eqtl_Ileum_chr6[grep(same_Ileum_chr6[i],eqtl_Ileum_chr6$loci),]
    chr6_Ileum_genes<-c(chr6_Ileum_genes,a)
    chr6_Ileum_overloci<-rbind(chr6_Ileum_overloci,b)
  }
}

chr7_Ileum_genes<-NULL
chr7_Ileum_overloci<-NULL
if(length(same_Ileum_chr7!=0)){
  for(i in 1:length(same_Ileum_chr7)){
    a<-eqtl_Ileum_chr7$phenotype_id[grep(same_Ileum_chr7[i],eqtl_Ileum_chr7$loci)]
    b<-eqtl_Ileum_chr7[grep(same_Ileum_chr7[i],eqtl_Ileum_chr7$loci),]
    chr7_Ileum_genes<-c(chr7_Ileum_genes,a)
    chr7_Ileum_overloci<-rbind(chr7_Ileum_overloci,b)
  }
}

chr8_Ileum_genes<-NULL
chr8_Ileum_overloci<-NULL
if(length(same_Ileum_chr8!=0)){
  for(i in 1:length(same_Ileum_chr8)){
    a<-eqtl_Ileum_chr8$phenotype_id[grep(same_Ileum_chr8[i],eqtl_Ileum_chr8$loci)]
    b<-eqtl_Ileum_chr8[grep(same_Ileum_chr8[i],eqtl_Ileum_chr8$loci),]
    chr8_Ileum_genes<-c(chr8_Ileum_genes,a)
    chr8_Ileum_overloci<-rbind(chr8_Ileum_overloci,b)
  }
}

if(length(same_Ileum_chr9!=0)){
  chr9_Ileum_genes<-NULL
  chr9_Ileum_overloci<-NULL
  for(i in 1:length(same_Ileum_chr9)){
    a<-eqtl_Ileum_chr9$phenotype_id[grep(same_Ileum_chr9[i],eqtl_Ileum_chr9$loci)]
    b<-eqtl_Ileum_chr9[grep(same_Ileum_chr9[i],eqtl_Ileum_chr9$loci),]
    chr9_Ileum_genes<-c(chr9_Ileum_genes,a)
    chr9_Ileum_overloci<-rbind(chr9_Ileum_overloci,b)
  }
}

chr10_Ileum_genes<-NULL
chr10_Ileum_overloci<-NULL
if(length(same_Ileum_chr10!=0)){
  for(i in 1:length(same_Ileum_chr10)){
    a<-eqtl_Ileum_chr10$phenotype_id[grep(same_Ileum_chr10[i],eqtl_Ileum_chr10$loci)]
    b<-eqtl_Ileum_chr10[grep(same_Ileum_chr10[i],eqtl_Ileum_chr10$loci),]
    chr10_Ileum_genes<-c(chr10_Ileum_genes,a)
    chr10_Ileum_overloci<-rbind(chr10_Ileum_overloci,b)
  }
}

chr11_Ileum_genes<-NULL
chr11_Ileum_overloci<-NULL
if(length(same_Ileum_chr11!=0)){
  for(i in 1:length(same_Ileum_chr11)){
    a<-eqtl_Ileum_chr11$phenotype_id[grep(same_Ileum_chr11[i],eqtl_Ileum_chr11$loci)]
    b<-eqtl_Ileum_chr11[grep(same_Ileum_chr11[i],eqtl_Ileum_chr11$loci),]
    chr11_Ileum_genes<-c(chr11_Ileum_genes,a)
    chr11_Ileum_overloci<-rbind(chr11_Ileum_overloci,b)
  }
}

chr12_Ileum_genes<-NULL
chr12_Ileum_overloci<-NULL
if(length(same_Ileum_chr12!=0)){
  for(i in 1:length(same_Ileum_chr12)){
    a<-eqtl_Ileum_chr12$phenotype_id[grep(same_Ileum_chr12[i],eqtl_Ileum_chr12$loci)]
    b<-eqtl_Ileum_chr12[grep(same_Ileum_chr12[i],eqtl_Ileum_chr12$loci),]
    chr12_Ileum_genes<-c(chr12_Ileum_genes,a)
    chr12_Ileum_overloci<-rbind(chr12_Ileum_overloci,b)
  }
}

chr13_Ileum_genes<-NULL
chr13_Ileum_overloci<-NULL
if(length(same_Ileum_chr13!=0)){
  for(i in 1:length(same_Ileum_chr13)){
    a<-eqtl_Ileum_chr13$phenotype_id[grep(same_Ileum_chr13[i],eqtl_Ileum_chr13$loci)]
    b<-eqtl_Ileum_chr13[grep(same_Ileum_chr13[i],eqtl_Ileum_chr13$loci),]
    chr13_Ileum_genes<-c(chr13_Ileum_genes,a)
    chr13_Ileum_overloci<-rbind(chr13_Ileum_overloci,b)
  }
}

chr14_Ileum_genes<-NULL
chr14_Ileum_overloci<-NULL
if(length(same_Ileum_chr14!=0)){
  for(i in 1:length(same_Ileum_chr14)){
    a<-eqtl_Ileum_chr14$phenotype_id[grep(same_Ileum_chr14[i],eqtl_Ileum_chr14$loci)]
    b<-eqtl_Ileum_chr14[grep(same_Ileum_chr14[i],eqtl_Ileum_chr14$loci),]
    chr14_Ileum_genes<-c(chr14_Ileum_genes,a)
    chr14_Ileum_overloci<-rbind(chr14_Ileum_overloci,b)
  }
}

chr15_Ileum_genes<-NULL
chr15_Ileum_overloci<-NULL
if(length(same_Ileum_chr15!=0)){
  for(i in 1:length(same_Ileum_chr15)){
    a<-eqtl_Ileum_chr15$phenotype_id[grep(same_Ileum_chr15[i],eqtl_Ileum_chr15$loci)]
    b<-eqtl_Ileum_chr15[grep(same_Ileum_chr15[i],eqtl_Ileum_chr15$loci),]
    chr15_Ileum_genes<-c(chr15_Ileum_genes,a)
    chr15_Ileum_overloci<-rbind(chr15_Ileum_overloci,b)
  }
}

chr16_Ileum_genes<-NULL
chr16_Ileum_overloci<-NULL
if(length(same_Ileum_chr16!=0)){
  for(i in 1:length(same_Ileum_chr16)){
    a<-eqtl_Ileum_chr16$phenotype_id[grep(same_Ileum_chr16[i],eqtl_Ileum_chr16$loci)]
    b<-eqtl_Ileum_chr16[grep(same_Ileum_chr16[i],eqtl_Ileum_chr16$loci),]
    chr16_Ileum_genes<-c(chr16_Ileum_genes,a)
    chr16_Ileum_overloci<-rbind(chr16_Ileum_overloci,b)
  }
}

chr17_Ileum_genes<-NULL
chr17_Ileum_overloci<-NULL
if(length(same_Ileum_chr17!=0)){
  for(i in 1:length(same_Ileum_chr17)){
    a<-eqtl_Ileum_chr17$phenotype_id[grep(same_Ileum_chr17[i],eqtl_Ileum_chr17$loci)]
    b<-eqtl_Ileum_chr17[grep(same_Ileum_chr17[i],eqtl_Ileum_chr17$loci),]
    chr17_Ileum_genes<-c(chr17_Ileum_genes,a)
    chr17_Ileum_overloci<-rbind(chr17_Ileum_overloci,b)
  }
}

chr18_Ileum_genes<-NULL
chr18_Ileum_overloci<-NULL
if(length(same_Ileum_chr18!=0)){
  for(i in 1:length(same_Ileum_chr18)){
    a<-eqtl_Ileum_chr18$phenotype_id[grep(same_Ileum_chr18[i],eqtl_Ileum_chr18$loci)]
    b<-eqtl_Ileum_chr18[grep(same_Ileum_chr18[i],eqtl_Ileum_chr18$loci),]
    chr18_Ileum_genes<-c(chr18_Ileum_genes,a)
    chr18_Ileum_overloci<-rbind(chr18_Ileum_overloci,b)
  }
}
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
one2one_pig<-annotation$`Pig gene stable ID`[which(annotation$`Pig homology type`=="ortholog_one2one")]

chr1_Ileum_overloci_one2one<-chr1_Ileum_overloci[match(intersect(chr1_Ileum_genes,one2one_pig),chr1_Ileum_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ileum_chr1$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ileum_chr1$hum_chr<-tmp$V1
Ileum_chr1$hum_loci<-tmp$V3
FM_Ileum_hum$index<-paste0(FM_Ileum_hum$chr,"_",FM_Ileum_hum$variant_pos)
Ileum_chr1$index<-paste0("chr",Ileum_chr1$hum_chr,"_",Ileum_chr1$hum_loci)
chr1_Ileum_hum<-NULL
for(i in 1:length(same_Ileum_chr1)){
  a<-Ileum_chr1[match(same_Ileum_chr1[i],Ileum_chr1$pig_pos),]
  b<-FM_Ileum_hum[match(intersect(a$index,FM_Ileum_hum$index),FM_Ileum_hum$index),]
  chr1_Ileum_hum<-rbind(chr1_Ileum_hum,b)
}
chr1_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr1_Ileum_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr1_Ileum_hum_one2one<-chr1_Ileum_hum[match(intersect(chr1_pig2hum_one2one,chr1_Ileum_hum$gene_id),chr1_Ileum_hum$gene_id),]

chr2_Ileum_overloci_one2one<-chr2_Ileum_overloci[match(intersect(chr2_Ileum_genes,one2one_pig),chr2_Ileum_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ileum_chr2$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ileum_chr2$hum_chr<-tmp$V1
Ileum_chr2$hum_loci<-tmp$V3
FM_Ileum_hum$index<-paste0(FM_Ileum_hum$chr,"_",FM_Ileum_hum$variant_pos)
Ileum_chr2$index<-paste0("chr",Ileum_chr2$hum_chr,"_",Ileum_chr2$hum_loci)
chr2_Ileum_hum<-NULL
for(i in 1:length(same_Ileum_chr2)){
  a<-Ileum_chr2[match(same_Ileum_chr2[i],Ileum_chr2$pig_pos),]
  b<-FM_Ileum_hum[match(intersect(a$index,FM_Ileum_hum$index),FM_Ileum_hum$index),]
  chr2_Ileum_hum<-rbind(chr2_Ileum_hum,b)
}
chr2_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr2_Ileum_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr2_Ileum_hum_one2one<-chr2_Ileum_hum[match(intersect(chr2_pig2hum_one2one,chr2_Ileum_hum$gene_id),chr2_Ileum_hum$gene_id),]

chr3_Ileum_overloci_one2one<-chr3_Ileum_overloci[match(intersect(chr3_Ileum_genes,one2one_pig),chr3_Ileum_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ileum_chr3$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ileum_chr3$hum_chr<-tmp$V1
Ileum_chr3$hum_loci<-tmp$V3
FM_Ileum_hum$index<-paste0(FM_Ileum_hum$chr,"_",FM_Ileum_hum$variant_pos)
Ileum_chr3$index<-paste0("chr",Ileum_chr3$hum_chr,"_",Ileum_chr3$hum_loci)
chr3_Ileum_hum<-NULL
for(i in 1:length(same_Ileum_chr3)){
  a<-Ileum_chr3[match(same_Ileum_chr3[i],Ileum_chr3$pig_pos),]
  b<-FM_Ileum_hum[match(intersect(a$index,FM_Ileum_hum$index),FM_Ileum_hum$index),]
  chr3_Ileum_hum<-rbind(chr3_Ileum_hum,b)
}
chr3_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr3_Ileum_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr3_Ileum_hum_one2one<-chr3_Ileum_hum[match(intersect(chr3_pig2hum_one2one,chr3_Ileum_hum$gene_id),chr3_Ileum_hum$gene_id),]

chr4_Ileum_overloci_one2one<-chr4_Ileum_overloci[match(intersect(chr4_Ileum_genes,one2one_pig),chr4_Ileum_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ileum_chr4$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ileum_chr4$hum_chr<-tmp$V1
Ileum_chr4$hum_loci<-tmp$V3
FM_Ileum_hum$index<-paste0(FM_Ileum_hum$chr,"_",FM_Ileum_hum$variant_pos)
Ileum_chr4$index<-paste0("chr",Ileum_chr4$hum_chr,"_",Ileum_chr4$hum_loci)
chr4_Ileum_hum<-NULL
for(i in 1:length(same_Ileum_chr4)){
  a<-Ileum_chr4[match(same_Ileum_chr4[i],Ileum_chr4$pig_pos),]
  b<-FM_Ileum_hum[match(intersect(a$index,FM_Ileum_hum$index),FM_Ileum_hum$index),]
  chr4_Ileum_hum<-rbind(chr4_Ileum_hum,b)
}
chr4_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr4_Ileum_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr4_Ileum_hum_one2one<-chr4_Ileum_hum[match(intersect(chr4_pig2hum_one2one,chr4_Ileum_hum$gene_id),chr4_Ileum_hum$gene_id),]

chr5_Ileum_overloci_one2one<-chr5_Ileum_overloci[match(intersect(chr5_Ileum_genes,one2one_pig),chr5_Ileum_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ileum_chr5$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ileum_chr5$hum_chr<-tmp$V1
Ileum_chr5$hum_loci<-tmp$V3
FM_Ileum_hum$index<-paste0(FM_Ileum_hum$chr,"_",FM_Ileum_hum$variant_pos)
Ileum_chr5$index<-paste0("chr",Ileum_chr5$hum_chr,"_",Ileum_chr5$hum_loci)
chr5_Ileum_hum<-NULL
for(i in 1:length(same_Ileum_chr5)){
  a<-Ileum_chr5[match(same_Ileum_chr5[i],Ileum_chr5$pig_pos),]
  b<-FM_Ileum_hum[match(intersect(a$index,FM_Ileum_hum$index),FM_Ileum_hum$index),]
  chr5_Ileum_hum<-rbind(chr5_Ileum_hum,b)
}
chr5_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr5_Ileum_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr5_Ileum_hum_one2one<-chr5_Ileum_hum[match(intersect(chr5_pig2hum_one2one,chr5_Ileum_hum$gene_id),chr5_Ileum_hum$gene_id),]

chr6_Ileum_overloci_one2one<-chr6_Ileum_overloci[match(intersect(chr6_Ileum_genes,one2one_pig),chr6_Ileum_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ileum_chr6$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ileum_chr6$hum_chr<-tmp$V1
Ileum_chr6$hum_loci<-tmp$V3
FM_Ileum_hum$index<-paste0(FM_Ileum_hum$chr,"_",FM_Ileum_hum$variant_pos)
Ileum_chr6$index<-paste0("chr",Ileum_chr6$hum_chr,"_",Ileum_chr6$hum_loci)
chr6_Ileum_hum<-NULL
for(i in 1:length(same_Ileum_chr6)){
  a<-Ileum_chr6[match(same_Ileum_chr6[i],Ileum_chr6$pig_pos),]
  b<-FM_Ileum_hum[match(intersect(a$index,FM_Ileum_hum$index),FM_Ileum_hum$index),]
  chr6_Ileum_hum<-rbind(chr6_Ileum_hum,b)
}
chr6_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr6_Ileum_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr6_Ileum_hum_one2one<-chr6_Ileum_hum[match(intersect(chr6_pig2hum_one2one,chr6_Ileum_hum$gene_id),chr6_Ileum_hum$gene_id),]

chr7_Ileum_overloci_one2one<-chr7_Ileum_overloci[match(intersect(chr7_Ileum_genes,one2one_pig),chr7_Ileum_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ileum_chr7$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ileum_chr7$hum_chr<-tmp$V1
Ileum_chr7$hum_loci<-tmp$V3
FM_Ileum_hum$index<-paste0(FM_Ileum_hum$chr,"_",FM_Ileum_hum$variant_pos)
Ileum_chr7$index<-paste0("chr",Ileum_chr7$hum_chr,"_",Ileum_chr7$hum_loci)
chr7_Ileum_hum<-NULL
for(i in 1:length(same_Ileum_chr7)){
  a<-Ileum_chr7[match(same_Ileum_chr7[i],Ileum_chr7$pig_pos),]
  b<-FM_Ileum_hum[match(intersect(a$index,FM_Ileum_hum$index),FM_Ileum_hum$index),]
  chr7_Ileum_hum<-rbind(chr7_Ileum_hum,b)
}
chr7_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr7_Ileum_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr7_Ileum_hum_one2one<-chr7_Ileum_hum[match(intersect(chr7_pig2hum_one2one,chr7_Ileum_hum$gene_id),chr7_Ileum_hum$gene_id),]

chr8_Ileum_overloci_one2one<-chr8_Ileum_overloci[match(intersect(chr8_Ileum_genes,one2one_pig),chr8_Ileum_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ileum_chr8$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ileum_chr8$hum_chr<-tmp$V1
Ileum_chr8$hum_loci<-tmp$V3
FM_Ileum_hum$index<-paste0(FM_Ileum_hum$chr,"_",FM_Ileum_hum$variant_pos)
Ileum_chr8$index<-paste0("chr",Ileum_chr8$hum_chr,"_",Ileum_chr8$hum_loci)
chr8_Ileum_hum<-NULL
for(i in 1:length(same_Ileum_chr8)){
  a<-Ileum_chr8[match(same_Ileum_chr8[i],Ileum_chr8$pig_pos),]
  b<-FM_Ileum_hum[match(intersect(a$index,FM_Ileum_hum$index),FM_Ileum_hum$index),]
  chr8_Ileum_hum<-rbind(chr8_Ileum_hum,b)
}
chr8_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr8_Ileum_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr8_Ileum_hum_one2one<-chr8_Ileum_hum[match(intersect(chr8_pig2hum_one2one,chr8_Ileum_hum$gene_id),chr8_Ileum_hum$gene_id),]

chr9_Ileum_overloci_one2one<-chr9_Ileum_overloci[match(intersect(chr9_Ileum_genes,one2one_pig),chr9_Ileum_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ileum_chr9$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ileum_chr9$hum_chr<-tmp$V1
Ileum_chr9$hum_loci<-tmp$V3
FM_Ileum_hum$index<-paste0(FM_Ileum_hum$chr,"_",FM_Ileum_hum$variant_pos)
Ileum_chr9$index<-paste0("chr",Ileum_chr9$hum_chr,"_",Ileum_chr9$hum_loci)
chr9_Ileum_hum<-NULL
for(i in 1:length(same_Ileum_chr9)){
  a<-Ileum_chr9[match(same_Ileum_chr9[i],Ileum_chr9$pig_pos),]
  b<-FM_Ileum_hum[match(intersect(a$index,FM_Ileum_hum$index),FM_Ileum_hum$index),]
  chr9_Ileum_hum<-rbind(chr9_Ileum_hum,b)
}
chr9_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr9_Ileum_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr9_Ileum_hum_one2one<-chr9_Ileum_hum[match(intersect(chr9_pig2hum_one2one,chr9_Ileum_hum$gene_id),chr9_Ileum_hum$gene_id),]

chr10_Ileum_overloci_one2one<-chr10_Ileum_overloci[match(intersect(chr10_Ileum_genes,one2one_pig),chr10_Ileum_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ileum_chr10$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ileum_chr10$hum_chr<-tmp$V1
Ileum_chr10$hum_loci<-tmp$V3
FM_Ileum_hum$index<-paste0(FM_Ileum_hum$chr,"_",FM_Ileum_hum$variant_pos)
Ileum_chr10$index<-paste0("chr",Ileum_chr10$hum_chr,"_",Ileum_chr10$hum_loci)
chr10_Ileum_hum<-NULL
for(i in 1:length(same_Ileum_chr10)){
  a<-Ileum_chr10[match(same_Ileum_chr10[i],Ileum_chr10$pig_pos),]
  b<-FM_Ileum_hum[match(intersect(a$index,FM_Ileum_hum$index),FM_Ileum_hum$index),]
  chr10_Ileum_hum<-rbind(chr10_Ileum_hum,b)
}
chr10_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr10_Ileum_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr10_Ileum_hum_one2one<-chr10_Ileum_hum[match(intersect(chr10_pig2hum_one2one,chr10_Ileum_hum$gene_id),chr10_Ileum_hum$gene_id),]

chr11_Ileum_overloci_one2one<-chr11_Ileum_overloci[match(intersect(chr11_Ileum_genes,one2one_pig),chr11_Ileum_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ileum_chr11$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ileum_chr11$hum_chr<-tmp$V1
Ileum_chr11$hum_loci<-tmp$V3
FM_Ileum_hum$index<-paste0(FM_Ileum_hum$chr,"_",FM_Ileum_hum$variant_pos)
Ileum_chr11$index<-paste0("chr",Ileum_chr11$hum_chr,"_",Ileum_chr11$hum_loci)
chr11_Ileum_hum<-NULL
for(i in 1:length(same_Ileum_chr11)){
  a<-Ileum_chr11[match(same_Ileum_chr11[i],Ileum_chr11$pig_pos),]
  b<-FM_Ileum_hum[match(intersect(a$index,FM_Ileum_hum$index),FM_Ileum_hum$index),]
  chr11_Ileum_hum<-rbind(chr11_Ileum_hum,b)
}
chr11_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr11_Ileum_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr11_Ileum_hum_one2one<-chr11_Ileum_hum[match(intersect(chr11_pig2hum_one2one,chr11_Ileum_hum$gene_id),chr11_Ileum_hum$gene_id),]

chr12_Ileum_overloci_one2one<-chr12_Ileum_overloci[match(intersect(chr12_Ileum_genes,one2one_pig),chr12_Ileum_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ileum_chr12$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ileum_chr12$hum_chr<-tmp$V1
Ileum_chr12$hum_loci<-tmp$V3
FM_Ileum_hum$index<-paste0(FM_Ileum_hum$chr,"_",FM_Ileum_hum$variant_pos)
Ileum_chr12$index<-paste0("chr",Ileum_chr12$hum_chr,"_",Ileum_chr12$hum_loci)
chr12_Ileum_hum<-NULL
for(i in 1:length(same_Ileum_chr12)){
  a<-Ileum_chr12[match(same_Ileum_chr12[i],Ileum_chr12$pig_pos),]
  b<-FM_Ileum_hum[match(intersect(a$index,FM_Ileum_hum$index),FM_Ileum_hum$index),]
  chr12_Ileum_hum<-rbind(chr12_Ileum_hum,b)
}
chr12_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr12_Ileum_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr12_Ileum_hum_one2one<-chr12_Ileum_hum[match(intersect(chr12_pig2hum_one2one,chr12_Ileum_hum$gene_id),chr12_Ileum_hum$gene_id),]

chr13_Ileum_overloci_one2one<-chr13_Ileum_overloci[match(intersect(chr13_Ileum_genes,one2one_pig),chr13_Ileum_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ileum_chr13$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ileum_chr13$hum_chr<-tmp$V1
Ileum_chr13$hum_loci<-tmp$V3
FM_Ileum_hum$index<-paste0(FM_Ileum_hum$chr,"_",FM_Ileum_hum$variant_pos)
Ileum_chr13$index<-paste0("chr",Ileum_chr13$hum_chr,"_",Ileum_chr13$hum_loci)
chr13_Ileum_hum<-NULL
for(i in 1:length(same_Ileum_chr13)){
  a<-Ileum_chr13[match(same_Ileum_chr13[i],Ileum_chr13$pig_pos),]
  b<-FM_Ileum_hum[match(intersect(a$index,FM_Ileum_hum$index),FM_Ileum_hum$index),]
  chr13_Ileum_hum<-rbind(chr13_Ileum_hum,b)
}
chr13_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr13_Ileum_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr13_Ileum_hum_one2one<-chr13_Ileum_hum[match(intersect(chr13_pig2hum_one2one,chr13_Ileum_hum$gene_id),chr13_Ileum_hum$gene_id),]

chr14_Ileum_overloci_one2one<-chr14_Ileum_overloci[match(intersect(chr14_Ileum_genes,one2one_pig),chr14_Ileum_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ileum_chr14$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ileum_chr14$hum_chr<-tmp$V1
Ileum_chr14$hum_loci<-tmp$V3
FM_Ileum_hum$index<-paste0(FM_Ileum_hum$chr,"_",FM_Ileum_hum$variant_pos)
Ileum_chr14$index<-paste0("chr",Ileum_chr14$hum_chr,"_",Ileum_chr14$hum_loci)
chr14_Ileum_hum<-NULL
for(i in 1:length(same_Ileum_chr14)){
  a<-Ileum_chr14[match(same_Ileum_chr14[i],Ileum_chr14$pig_pos),]
  b<-FM_Ileum_hum[match(intersect(a$index,FM_Ileum_hum$index),FM_Ileum_hum$index),]
  chr14_Ileum_hum<-rbind(chr14_Ileum_hum,b)
}
chr14_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr14_Ileum_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr14_Ileum_hum_one2one<-chr14_Ileum_hum[match(intersect(chr14_pig2hum_one2one,chr14_Ileum_hum$gene_id),chr14_Ileum_hum$gene_id),]

chr15_Ileum_overloci_one2one<-chr15_Ileum_overloci[match(intersect(chr15_Ileum_genes,one2one_pig),chr15_Ileum_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ileum_chr15$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ileum_chr15$hum_chr<-tmp$V1
Ileum_chr15$hum_loci<-tmp$V3
FM_Ileum_hum$index<-paste0(FM_Ileum_hum$chr,"_",FM_Ileum_hum$variant_pos)
Ileum_chr15$index<-paste0("chr",Ileum_chr15$hum_chr,"_",Ileum_chr15$hum_loci)
chr15_Ileum_hum<-NULL
for(i in 1:length(same_Ileum_chr15)){
  a<-Ileum_chr15[match(same_Ileum_chr15[i],Ileum_chr15$pig_pos),]
  b<-FM_Ileum_hum[match(intersect(a$index,FM_Ileum_hum$index),FM_Ileum_hum$index),]
  chr15_Ileum_hum<-rbind(chr15_Ileum_hum,b)
}
chr15_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr15_Ileum_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr15_Ileum_hum_one2one<-chr15_Ileum_hum[match(intersect(chr15_pig2hum_one2one,chr15_Ileum_hum$gene_id),chr15_Ileum_hum$gene_id),]

chr16_Ileum_overloci_one2one<-chr16_Ileum_overloci[match(intersect(chr16_Ileum_genes,one2one_pig),chr16_Ileum_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ileum_chr16$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ileum_chr16$hum_chr<-tmp$V1
Ileum_chr16$hum_loci<-tmp$V3
FM_Ileum_hum$index<-paste0(FM_Ileum_hum$chr,"_",FM_Ileum_hum$variant_pos)
Ileum_chr16$index<-paste0("chr",Ileum_chr16$hum_chr,"_",Ileum_chr16$hum_loci)
chr16_Ileum_hum<-NULL
for(i in 1:length(same_Ileum_chr16)){
  a<-Ileum_chr16[match(same_Ileum_chr16[i],Ileum_chr16$pig_pos),]
  b<-FM_Ileum_hum[match(intersect(a$index,FM_Ileum_hum$index),FM_Ileum_hum$index),]
  chr16_Ileum_hum<-rbind(chr16_Ileum_hum,b)
}
chr16_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr16_Ileum_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr16_Ileum_hum_one2one<-chr16_Ileum_hum[match(intersect(chr16_pig2hum_one2one,chr16_Ileum_hum$gene_id),chr16_Ileum_hum$gene_id),]

chr17_Ileum_overloci_one2one<-chr17_Ileum_overloci[match(intersect(chr17_Ileum_genes,one2one_pig),chr17_Ileum_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ileum_chr17$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ileum_chr17$hum_chr<-tmp$V1
Ileum_chr17$hum_loci<-tmp$V3
FM_Ileum_hum$index<-paste0(FM_Ileum_hum$chr,"_",FM_Ileum_hum$variant_pos)
Ileum_chr17$index<-paste0("chr",Ileum_chr17$hum_chr,"_",Ileum_chr17$hum_loci)
chr17_Ileum_hum<-NULL
for(i in 1:length(same_Ileum_chr17)){
  a<-Ileum_chr17[match(same_Ileum_chr17[i],Ileum_chr17$pig_pos),]
  b<-FM_Ileum_hum[match(intersect(a$index,FM_Ileum_hum$index),FM_Ileum_hum$index),]
  chr17_Ileum_hum<-rbind(chr17_Ileum_hum,b)
}
chr17_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr17_Ileum_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr17_Ileum_hum_one2one<-chr17_Ileum_hum[match(intersect(chr17_pig2hum_one2one,chr17_Ileum_hum$gene_id),chr17_Ileum_hum$gene_id),]

chr18_Ileum_overloci_one2one<-chr18_Ileum_overloci[match(intersect(chr18_Ileum_genes,one2one_pig),chr18_Ileum_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ileum_chr18$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ileum_chr18$hum_chr<-tmp$V1
Ileum_chr18$hum_loci<-tmp$V3
FM_Ileum_hum$index<-paste0(FM_Ileum_hum$chr,"_",FM_Ileum_hum$variant_pos)
Ileum_chr18$index<-paste0("chr",Ileum_chr18$hum_chr,"_",Ileum_chr18$hum_loci)
chr18_Ileum_hum<-NULL
for(i in 1:length(same_Ileum_chr18)){
  a<-Ileum_chr18[match(same_Ileum_chr18[i],Ileum_chr18$pig_pos),]
  b<-FM_Ileum_hum[match(intersect(a$index,FM_Ileum_hum$index),FM_Ileum_hum$index),]
  chr18_Ileum_hum<-rbind(chr18_Ileum_hum,b)
}
chr18_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr18_Ileum_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr18_Ileum_hum_one2one<-chr18_Ileum_hum[match(intersect(chr18_pig2hum_one2one,chr18_Ileum_hum$gene_id),chr18_Ileum_hum$gene_id),]

Ileum_one2one_SNP_hum<-rbind(chr1_Ileum_hum_one2one,chr2_Ileum_hum_one2one,chr3_Ileum_hum_one2one,chr4_Ileum_hum_one2one,chr5_Ileum_hum_one2one,
                             chr6_Ileum_hum_one2one,chr7_Ileum_hum_one2one,chr8_Ileum_hum_one2one,chr9_Ileum_hum_one2one,chr10_Ileum_hum_one2one,
                             chr11_Ileum_hum_one2one,chr12_Ileum_hum_one2one,chr13_Ileum_hum_one2one,chr14_Ileum_hum_one2one,chr15_Ileum_hum_one2one,
                             chr16_Ileum_hum_one2one,chr17_Ileum_hum_one2one,chr18_Ileum_hum_one2one)

chr1_Ileum_pig_one2one<-chr1_Ileum_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr1_Ileum_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr1_Ileum_overloci_one2one$phenotype_id),1:9]
chr2_Ileum_pig_one2one<-chr2_Ileum_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr2_Ileum_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr2_Ileum_overloci_one2one$phenotype_id),1:9]
chr3_Ileum_pig_one2one<-chr3_Ileum_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr3_Ileum_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr3_Ileum_overloci_one2one$phenotype_id),1:9]
chr4_Ileum_pig_one2one<-chr4_Ileum_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr4_Ileum_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr4_Ileum_overloci_one2one$phenotype_id),1:9]
chr5_Ileum_pig_one2one<-chr5_Ileum_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr5_Ileum_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr5_Ileum_overloci_one2one$phenotype_id),1:9]
chr6_Ileum_pig_one2one<-chr6_Ileum_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr6_Ileum_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr6_Ileum_overloci_one2one$phenotype_id),1:9]
chr7_Ileum_pig_one2one<-chr7_Ileum_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr7_Ileum_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr7_Ileum_overloci_one2one$phenotype_id),1:9]
chr8_Ileum_pig_one2one<-chr8_Ileum_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr8_Ileum_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr8_Ileum_overloci_one2one$phenotype_id),1:9]
chr9_Ileum_pig_one2one<-chr9_Ileum_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr9_Ileum_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr9_Ileum_overloci_one2one$phenotype_id),1:9]
chr10_Ileum_pig_one2one<-chr10_Ileum_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr10_Ileum_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr10_Ileum_overloci_one2one$phenotype_id),1:9]
chr11_Ileum_pig_one2one<-chr11_Ileum_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr11_Ileum_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr11_Ileum_overloci_one2one$phenotype_id),1:9]
chr12_Ileum_pig_one2one<-chr12_Ileum_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr12_Ileum_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr12_Ileum_overloci_one2one$phenotype_id),1:9]
chr13_Ileum_pig_one2one<-chr13_Ileum_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr13_Ileum_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr13_Ileum_overloci_one2one$phenotype_id),1:9]
chr14_Ileum_pig_one2one<-chr14_Ileum_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr14_Ileum_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr14_Ileum_overloci_one2one$phenotype_id),1:9]
chr15_Ileum_pig_one2one<-chr15_Ileum_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr15_Ileum_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr15_Ileum_overloci_one2one$phenotype_id),1:9]
chr16_Ileum_pig_one2one<-chr16_Ileum_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr16_Ileum_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr16_Ileum_overloci_one2one$phenotype_id),1:9]
chr17_Ileum_pig_one2one<-chr17_Ileum_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr17_Ileum_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr17_Ileum_overloci_one2one$phenotype_id),1:9]
chr18_Ileum_pig_one2one<-chr18_Ileum_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr18_Ileum_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr18_Ileum_overloci_one2one$phenotype_id),1:9]

Ileum_one2one_SNP_pig<-rbind(chr1_Ileum_pig_one2one,chr2_Ileum_pig_one2one,chr3_Ileum_pig_one2one,chr4_Ileum_pig_one2one,chr5_Ileum_pig_one2one,
                             chr6_Ileum_pig_one2one,chr7_Ileum_pig_one2one,chr8_Ileum_pig_one2one,chr9_Ileum_pig_one2one,chr10_Ileum_pig_one2one,
                             chr11_Ileum_pig_one2one,chr12_Ileum_pig_one2one,chr13_Ileum_pig_one2one,chr14_Ileum_pig_one2one,chr15_Ileum_pig_one2one,
                             chr16_Ileum_pig_one2one,chr17_Ileum_pig_one2one,chr18_Ileum_pig_one2one)

Ileum_SNP_sum<-array(NA,dim=c(nrow(Ileum_one2one_SNP_hum),2))
colnames(Ileum_SNP_sum)<-c("Human","Pig")
Ileum_SNP_sum<-as.data.frame(Ileum_SNP_sum)
Ileum_SNP_sum$Human<-Ileum_one2one_SNP_hum$slope / Ileum_one2one_SNP_hum$slope_se
Ileum_SNP_sum$Pig<-Ileum_one2one_SNP_pig$slope / Ileum_one2one_SNP_pig$slope_se
cor<-cor(abs(Ileum_SNP_sum$Human),abs(Ileum_SNP_sum$Pig))
p_val<-t.test(abs(Ileum_SNP_sum$Human),abs(Ileum_SNP_sum$Pig))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/Ileum_SNP_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(Ileum_SNP_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The number of eGenes expressed in common tissues",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

save(Ileum_one2one_SNP_hum,Ileum_one2one_SNP_pig,file="/Users/baizhonghao/Downloads/human-pig GTEx/Ileum_SNP.Rdata")

#Kidney_SNP_overlaploci#
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

eqtl_Kidney_chr1<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Kidney/Kidney.cis_qtl_pairs.1.txt"))
eqtl_Kidney_chr2<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Kidney/Kidney.cis_qtl_pairs.2.txt"))
eqtl_Kidney_chr3<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Kidney/Kidney.cis_qtl_pairs.3.txt"))
eqtl_Kidney_chr4<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Kidney/Kidney.cis_qtl_pairs.4.txt"))
eqtl_Kidney_chr5<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Kidney/Kidney.cis_qtl_pairs.5.txt"))
eqtl_Kidney_chr6<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Kidney/Kidney.cis_qtl_pairs.6.txt"))
eqtl_Kidney_chr7<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Kidney/Kidney.cis_qtl_pairs.7.txt"))
eqtl_Kidney_chr8<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Kidney/Kidney.cis_qtl_pairs.8.txt"))
eqtl_Kidney_chr9<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Kidney/Kidney.cis_qtl_pairs.9.txt"))
eqtl_Kidney_chr10<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Kidney/Kidney.cis_qtl_pairs.10.txt"))
eqtl_Kidney_chr11<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Kidney/Kidney.cis_qtl_pairs.11.txt"))
eqtl_Kidney_chr12<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Kidney/Kidney.cis_qtl_pairs.12.txt"))
eqtl_Kidney_chr13<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Kidney/Kidney.cis_qtl_pairs.13.txt"))
eqtl_Kidney_chr14<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Kidney/Kidney.cis_qtl_pairs.14.txt"))
eqtl_Kidney_chr15<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Kidney/Kidney.cis_qtl_pairs.15.txt"))
eqtl_Kidney_chr16<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Kidney/Kidney.cis_qtl_pairs.16.txt"))
eqtl_Kidney_chr17<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Kidney/Kidney.cis_qtl_pairs.17.txt"))
eqtl_Kidney_chr18<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Kidney/Kidney.cis_qtl_pairs.18.txt"))

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney_chr1$variant_id,split="_"))))
eqtl_Kidney_chr1$chr<-tmp$V1
eqtl_Kidney_chr1$loci<-tmp$V2

eqtl_Kidney_chr1$index<-paste0(eqtl_Kidney_chr1$chr,"-",eqtl_Kidney_chr1$loci)
Kidney_chr1$index<-paste0(Kidney_chr1$chr,"-",Kidney_chr1$pig_pos)
same_Kidney_chr1<-intersect(Kidney_chr1$pig_pos,eqtl_Kidney_chr1$loci)


tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney_chr2$variant_id,split="_"))))
eqtl_Kidney_chr2$chr<-tmp$V1
eqtl_Kidney_chr2$loci<-tmp$V2

eqtl_Kidney_chr2$index<-paste0(eqtl_Kidney_chr2$chr,"-",eqtl_Kidney_chr2$loci)
Kidney_chr2$index<-paste0(Kidney_chr2$chr,"-",Kidney_chr2$pig_pos)
same_Kidney_chr2<-intersect(Kidney_chr2$pig_pos,eqtl_Kidney_chr2$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney_chr3$variant_id,split="_"))))
eqtl_Kidney_chr3$chr<-tmp$V1
eqtl_Kidney_chr3$loci<-tmp$V2

eqtl_Kidney_chr3$index<-paste0(eqtl_Kidney_chr3$chr,"-",eqtl_Kidney_chr3$loci)
Kidney_chr3$index<-paste0(Kidney_chr3$chr,"-",Kidney_chr3$pig_pos)
same_Kidney_chr3<-intersect(Kidney_chr3$pig_pos,eqtl_Kidney_chr3$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney_chr4$variant_id,split="_"))))
eqtl_Kidney_chr4$chr<-tmp$V1
eqtl_Kidney_chr4$loci<-tmp$V2

eqtl_Kidney_chr4$index<-paste0(eqtl_Kidney_chr4$chr,"-",eqtl_Kidney_chr4$loci)
Kidney_chr4$index<-paste0(Kidney_chr4$chr,"-",Kidney_chr4$pig_pos)
same_Kidney_chr4<-intersect(Kidney_chr4$pig_pos,eqtl_Kidney_chr4$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney_chr5$variant_id,split="_"))))
eqtl_Kidney_chr5$loci<-tmp$V2
same_Kidney_chr5<-intersect(Kidney_chr5$pig_pos,eqtl_Kidney_chr5$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney_chr6$variant_id,split="_"))))
eqtl_Kidney_chr6$loci<-tmp$V2
same_Kidney_chr6<-intersect(Kidney_chr6$pig_pos,eqtl_Kidney_chr6$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney_chr7$variant_id,split="_"))))
eqtl_Kidney_chr7$loci<-tmp$V2
same_Kidney_chr7<-intersect(Kidney_chr7$pig_pos,eqtl_Kidney_chr7$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney_chr8$variant_id,split="_"))))
eqtl_Kidney_chr8$loci<-tmp$V2
same_Kidney_chr8<-intersect(Kidney_chr8$pig_pos,eqtl_Kidney_chr8$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney_chr9$variant_id,split="_"))))
eqtl_Kidney_chr9$loci<-tmp$V2
same_Kidney_chr9<-intersect(Kidney_chr9$pig_pos,eqtl_Kidney_chr9$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney_chr10$variant_id,split="_"))))
eqtl_Kidney_chr10$loci<-tmp$V2
same_Kidney_chr10<-intersect(Kidney_chr10$pig_pos,eqtl_Kidney_chr10$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney_chr11$variant_id,split="_"))))
eqtl_Kidney_chr11$loci<-tmp$V2
same_Kidney_chr11<-intersect(Kidney_chr11$pig_pos,eqtl_Kidney_chr11$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney_chr12$variant_id,split="_"))))
eqtl_Kidney_chr12$loci<-tmp$V2
same_Kidney_chr12<-intersect(Kidney_chr12$pig_pos,eqtl_Kidney_chr12$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney_chr13$variant_id,split="_"))))
eqtl_Kidney_chr13$loci<-tmp$V2
same_Kidney_chr13<-intersect(Kidney_chr13$pig_pos,eqtl_Kidney_chr13$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney_chr14$variant_id,split="_"))))
eqtl_Kidney_chr14$loci<-tmp$V2
same_Kidney_chr14<-intersect(Kidney_chr14$pig_pos,eqtl_Kidney_chr14$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney_chr15$variant_id,split="_"))))
eqtl_Kidney_chr15$loci<-tmp$V2
same_Kidney_chr15<-intersect(Kidney_chr15$pig_pos,eqtl_Kidney_chr15$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney_chr16$variant_id,split="_"))))
eqtl_Kidney_chr16$loci<-tmp$V2
same_Kidney_chr16<-intersect(Kidney_chr16$pig_pos,eqtl_Kidney_chr16$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney_chr17$variant_id,split="_"))))
eqtl_Kidney_chr17$loci<-tmp$V2
same_Kidney_chr17<-intersect(Kidney_chr17$pig_pos,eqtl_Kidney_chr17$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Kidney_chr18$variant_id,split="_"))))
eqtl_Kidney_chr18$loci<-tmp$V2
same_Kidney_chr18<-intersect(Kidney_chr18$pig_pos,eqtl_Kidney_chr18$loci)

chr1_Kidney_genes<-NULL
chr1_Kidney_overloci<-NULL
if(length(same_Kidney_chr1!=0)){
  for(i in 1:length(same_Kidney_chr1)){
    a<-eqtl_Kidney_chr1$phenotype_id[grep(same_Kidney_chr1[i],eqtl_Kidney_chr1$loci)]
    b<-eqtl_Kidney_chr1[grep(same_Kidney_chr1[i],eqtl_Kidney_chr1$loci),]
    chr1_Kidney_genes<-c(chr1_Kidney_genes,a)
    chr1_Kidney_overloci<-rbind(chr1_Kidney_overloci,b)
  }
}

chr2_Kidney_genes<-NULL
chr2_Kidney_overloci<-NULL
if(length(same_Kidney_chr2!=0)){
  for(i in 1:length(same_Kidney_chr2)){
    a<-eqtl_Kidney_chr2$phenotype_id[grep(same_Kidney_chr2[i],eqtl_Kidney_chr2$loci)]
    b<-eqtl_Kidney_chr2[grep(same_Kidney_chr2[i],eqtl_Kidney_chr2$loci),]
    chr2_Kidney_genes<-c(chr2_Kidney_genes,a)
    chr2_Kidney_overloci<-rbind(chr2_Kidney_overloci,b)
  }
}

chr3_Kidney_genes<-NULL
chr3_Kidney_overloci<-NULL
if(length(same_Kidney_chr3!=0)){
  for(i in 1:length(same_Kidney_chr3)){
    a<-eqtl_Kidney_chr3$phenotype_id[grep(same_Kidney_chr3[i],eqtl_Kidney_chr3$loci)]
    b<-eqtl_Kidney_chr3[grep(same_Kidney_chr3[i],eqtl_Kidney_chr3$loci),]
    chr3_Kidney_genes<-c(chr3_Kidney_genes,a)
    chr3_Kidney_overloci<-rbind(chr3_Kidney_overloci,b)
  }
}

chr4_Kidney_genes<-NULL
chr4_Kidney_overloci<-NULL
if(length(same_Kidney_chr4!=0)){
  for(i in 1:length(same_Kidney_chr4)){
    a<-eqtl_Kidney_chr4$phenotype_id[grep(same_Kidney_chr4[i],eqtl_Kidney_chr4$loci)]
    b<-eqtl_Kidney_chr4[grep(same_Kidney_chr4[i],eqtl_Kidney_chr4$loci),]
    chr4_Kidney_genes<-c(chr4_Kidney_genes,a)
    chr4_Kidney_overloci<-rbind(chr4_Kidney_overloci,b)
  }
}

chr5_Kidney_genes<-NULL
chr5_Kidney_overloci<-NULL
if(length(same_Kidney_chr5!=0)){
  for(i in 1:length(same_Kidney_chr5)){
    a<-eqtl_Kidney_chr5$phenotype_id[grep(same_Kidney_chr5[i],eqtl_Kidney_chr5$loci)]
    b<-eqtl_Kidney_chr5[grep(same_Kidney_chr5[i],eqtl_Kidney_chr5$loci),]
    chr5_Kidney_genes<-c(chr5_Kidney_genes,a)
    chr5_Kidney_overloci<-rbind(chr5_Kidney_overloci,b)
  }
}

chr6_Kidney_genes<-NULL
chr6_Kidney_overloci<-NULL
if(length(same_Kidney_chr6!=0)){
  for(i in 1:length(same_Kidney_chr6)){
    a<-eqtl_Kidney_chr6$phenotype_id[grep(same_Kidney_chr6[i],eqtl_Kidney_chr6$loci)]
    b<-eqtl_Kidney_chr6[grep(same_Kidney_chr6[i],eqtl_Kidney_chr6$loci),]
    chr6_Kidney_genes<-c(chr6_Kidney_genes,a)
    chr6_Kidney_overloci<-rbind(chr6_Kidney_overloci,b)
  }
}

chr7_Kidney_genes<-NULL
chr7_Kidney_overloci<-NULL
if(length(same_Kidney_chr7!=0)){
  for(i in 1:length(same_Kidney_chr7)){
    a<-eqtl_Kidney_chr7$phenotype_id[grep(same_Kidney_chr7[i],eqtl_Kidney_chr7$loci)]
    b<-eqtl_Kidney_chr7[grep(same_Kidney_chr7[i],eqtl_Kidney_chr7$loci),]
    chr7_Kidney_genes<-c(chr7_Kidney_genes,a)
    chr7_Kidney_overloci<-rbind(chr7_Kidney_overloci,b)
  }
}

chr8_Kidney_genes<-NULL
chr8_Kidney_overloci<-NULL
if(length(same_Kidney_chr8!=0)){
  for(i in 1:length(same_Kidney_chr8)){
    a<-eqtl_Kidney_chr8$phenotype_id[grep(same_Kidney_chr8[i],eqtl_Kidney_chr8$loci)]
    b<-eqtl_Kidney_chr8[grep(same_Kidney_chr8[i],eqtl_Kidney_chr8$loci),]
    chr8_Kidney_genes<-c(chr8_Kidney_genes,a)
    chr8_Kidney_overloci<-rbind(chr8_Kidney_overloci,b)
  }
}

if(length(same_Kidney_chr9!=0)){
  chr9_Kidney_genes<-NULL
  chr9_Kidney_overloci<-NULL
  for(i in 1:length(same_Kidney_chr9)){
    a<-eqtl_Kidney_chr9$phenotype_id[grep(same_Kidney_chr9[i],eqtl_Kidney_chr9$loci)]
    b<-eqtl_Kidney_chr9[grep(same_Kidney_chr9[i],eqtl_Kidney_chr9$loci),]
    chr9_Kidney_genes<-c(chr9_Kidney_genes,a)
    chr9_Kidney_overloci<-rbind(chr9_Kidney_overloci,b)
  }
}

chr10_Kidney_genes<-NULL
chr10_Kidney_overloci<-NULL
if(length(same_Kidney_chr10!=0)){
  for(i in 1:length(same_Kidney_chr10)){
    a<-eqtl_Kidney_chr10$phenotype_id[grep(same_Kidney_chr10[i],eqtl_Kidney_chr10$loci)]
    b<-eqtl_Kidney_chr10[grep(same_Kidney_chr10[i],eqtl_Kidney_chr10$loci),]
    chr10_Kidney_genes<-c(chr10_Kidney_genes,a)
    chr10_Kidney_overloci<-rbind(chr10_Kidney_overloci,b)
  }
}

chr11_Kidney_genes<-NULL
chr11_Kidney_overloci<-NULL
if(length(same_Kidney_chr11!=0)){
  for(i in 1:length(same_Kidney_chr11)){
    a<-eqtl_Kidney_chr11$phenotype_id[grep(same_Kidney_chr11[i],eqtl_Kidney_chr11$loci)]
    b<-eqtl_Kidney_chr11[grep(same_Kidney_chr11[i],eqtl_Kidney_chr11$loci),]
    chr11_Kidney_genes<-c(chr11_Kidney_genes,a)
    chr11_Kidney_overloci<-rbind(chr11_Kidney_overloci,b)
  }
}

chr12_Kidney_genes<-NULL
chr12_Kidney_overloci<-NULL
if(length(same_Kidney_chr12!=0)){
  for(i in 1:length(same_Kidney_chr12)){
    a<-eqtl_Kidney_chr12$phenotype_id[grep(same_Kidney_chr12[i],eqtl_Kidney_chr12$loci)]
    b<-eqtl_Kidney_chr12[grep(same_Kidney_chr12[i],eqtl_Kidney_chr12$loci),]
    chr12_Kidney_genes<-c(chr12_Kidney_genes,a)
    chr12_Kidney_overloci<-rbind(chr12_Kidney_overloci,b)
  }
}

chr13_Kidney_genes<-NULL
chr13_Kidney_overloci<-NULL
if(length(same_Kidney_chr13!=0)){
  for(i in 1:length(same_Kidney_chr13)){
    a<-eqtl_Kidney_chr13$phenotype_id[grep(same_Kidney_chr13[i],eqtl_Kidney_chr13$loci)]
    b<-eqtl_Kidney_chr13[grep(same_Kidney_chr13[i],eqtl_Kidney_chr13$loci),]
    chr13_Kidney_genes<-c(chr13_Kidney_genes,a)
    chr13_Kidney_overloci<-rbind(chr13_Kidney_overloci,b)
  }
}

chr14_Kidney_genes<-NULL
chr14_Kidney_overloci<-NULL
if(length(same_Kidney_chr14!=0)){
  for(i in 1:length(same_Kidney_chr14)){
    a<-eqtl_Kidney_chr14$phenotype_id[grep(same_Kidney_chr14[i],eqtl_Kidney_chr14$loci)]
    b<-eqtl_Kidney_chr14[grep(same_Kidney_chr14[i],eqtl_Kidney_chr14$loci),]
    chr14_Kidney_genes<-c(chr14_Kidney_genes,a)
    chr14_Kidney_overloci<-rbind(chr14_Kidney_overloci,b)
  }
}

chr15_Kidney_genes<-NULL
chr15_Kidney_overloci<-NULL
if(length(same_Kidney_chr15!=0)){
  for(i in 1:length(same_Kidney_chr15)){
    a<-eqtl_Kidney_chr15$phenotype_id[grep(same_Kidney_chr15[i],eqtl_Kidney_chr15$loci)]
    b<-eqtl_Kidney_chr15[grep(same_Kidney_chr15[i],eqtl_Kidney_chr15$loci),]
    chr15_Kidney_genes<-c(chr15_Kidney_genes,a)
    chr15_Kidney_overloci<-rbind(chr15_Kidney_overloci,b)
  }
}

chr16_Kidney_genes<-NULL
chr16_Kidney_overloci<-NULL
if(length(same_Kidney_chr16!=0)){
  for(i in 1:length(same_Kidney_chr16)){
    a<-eqtl_Kidney_chr16$phenotype_id[grep(same_Kidney_chr16[i],eqtl_Kidney_chr16$loci)]
    b<-eqtl_Kidney_chr16[grep(same_Kidney_chr16[i],eqtl_Kidney_chr16$loci),]
    chr16_Kidney_genes<-c(chr16_Kidney_genes,a)
    chr16_Kidney_overloci<-rbind(chr16_Kidney_overloci,b)
  }
}

chr17_Kidney_genes<-NULL
chr17_Kidney_overloci<-NULL
if(length(same_Kidney_chr17!=0)){
  for(i in 1:length(same_Kidney_chr17)){
    a<-eqtl_Kidney_chr17$phenotype_id[grep(same_Kidney_chr17[i],eqtl_Kidney_chr17$loci)]
    b<-eqtl_Kidney_chr17[grep(same_Kidney_chr17[i],eqtl_Kidney_chr17$loci),]
    chr17_Kidney_genes<-c(chr17_Kidney_genes,a)
    chr17_Kidney_overloci<-rbind(chr17_Kidney_overloci,b)
  }
}

chr18_Kidney_genes<-NULL
chr18_Kidney_overloci<-NULL
if(length(same_Kidney_chr18!=0)){
  for(i in 1:length(same_Kidney_chr18)){
    a<-eqtl_Kidney_chr18$phenotype_id[grep(same_Kidney_chr18[i],eqtl_Kidney_chr18$loci)]
    b<-eqtl_Kidney_chr18[grep(same_Kidney_chr18[i],eqtl_Kidney_chr18$loci),]
    chr18_Kidney_genes<-c(chr18_Kidney_genes,a)
    chr18_Kidney_overloci<-rbind(chr18_Kidney_overloci,b)
  }
}
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
one2one_pig<-annotation$`Pig gene stable ID`[which(annotation$`Pig homology type`=="ortholog_one2one")]

chr1_Kidney_overloci_one2one<-chr1_Kidney_overloci[match(intersect(chr1_Kidney_genes,one2one_pig),chr1_Kidney_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Kidney_chr1$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Kidney_chr1$hum_chr<-tmp$V1
Kidney_chr1$hum_loci<-tmp$V3
FM_Kidney_hum$index<-paste0(FM_Kidney_hum$chr,"_",FM_Kidney_hum$variant_pos)
Kidney_chr1$index<-paste0("chr",Kidney_chr1$hum_chr,"_",Kidney_chr1$hum_loci)
chr1_Kidney_hum<-NULL
for(i in 1:length(same_Kidney_chr1)){
  a<-Kidney_chr1[match(same_Kidney_chr1[i],Kidney_chr1$pig_pos),]
  b<-FM_Kidney_hum[match(intersect(a$index,FM_Kidney_hum$index),FM_Kidney_hum$index),]
  chr1_Kidney_hum<-rbind(chr1_Kidney_hum,b)
}
chr1_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr1_Kidney_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr1_Kidney_hum_one2one<-chr1_Kidney_hum[match(intersect(chr1_pig2hum_one2one,chr1_Kidney_hum$gene_id),chr1_Kidney_hum$gene_id),]

chr2_Kidney_overloci_one2one<-chr2_Kidney_overloci[match(intersect(chr2_Kidney_genes,one2one_pig),chr2_Kidney_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Kidney_chr2$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Kidney_chr2$hum_chr<-tmp$V1
Kidney_chr2$hum_loci<-tmp$V3
FM_Kidney_hum$index<-paste0(FM_Kidney_hum$chr,"_",FM_Kidney_hum$variant_pos)
Kidney_chr2$index<-paste0("chr",Kidney_chr2$hum_chr,"_",Kidney_chr2$hum_loci)
chr2_Kidney_hum<-NULL
for(i in 1:length(same_Kidney_chr2)){
  a<-Kidney_chr2[match(same_Kidney_chr2[i],Kidney_chr2$pig_pos),]
  b<-FM_Kidney_hum[match(intersect(a$index,FM_Kidney_hum$index),FM_Kidney_hum$index),]
  chr2_Kidney_hum<-rbind(chr2_Kidney_hum,b)
}
chr2_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr2_Kidney_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr2_Kidney_hum_one2one<-chr2_Kidney_hum[match(intersect(chr2_pig2hum_one2one,chr2_Kidney_hum$gene_id),chr2_Kidney_hum$gene_id),]

chr3_Kidney_overloci_one2one<-chr3_Kidney_overloci[match(intersect(chr3_Kidney_genes,one2one_pig),chr3_Kidney_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Kidney_chr3$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Kidney_chr3$hum_chr<-tmp$V1
Kidney_chr3$hum_loci<-tmp$V3
FM_Kidney_hum$index<-paste0(FM_Kidney_hum$chr,"_",FM_Kidney_hum$variant_pos)
Kidney_chr3$index<-paste0("chr",Kidney_chr3$hum_chr,"_",Kidney_chr3$hum_loci)
chr3_Kidney_hum<-NULL
for(i in 1:length(same_Kidney_chr3)){
  a<-Kidney_chr3[match(same_Kidney_chr3[i],Kidney_chr3$pig_pos),]
  b<-FM_Kidney_hum[match(intersect(a$index,FM_Kidney_hum$index),FM_Kidney_hum$index),]
  chr3_Kidney_hum<-rbind(chr3_Kidney_hum,b)
}
chr3_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr3_Kidney_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr3_Kidney_hum_one2one<-chr3_Kidney_hum[match(intersect(chr3_pig2hum_one2one,chr3_Kidney_hum$gene_id),chr3_Kidney_hum$gene_id),]

chr4_Kidney_overloci_one2one<-chr4_Kidney_overloci[match(intersect(chr4_Kidney_genes,one2one_pig),chr4_Kidney_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Kidney_chr4$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Kidney_chr4$hum_chr<-tmp$V1
Kidney_chr4$hum_loci<-tmp$V3
FM_Kidney_hum$index<-paste0(FM_Kidney_hum$chr,"_",FM_Kidney_hum$variant_pos)
Kidney_chr4$index<-paste0("chr",Kidney_chr4$hum_chr,"_",Kidney_chr4$hum_loci)
chr4_Kidney_hum<-NULL
for(i in 1:length(same_Kidney_chr4)){
  a<-Kidney_chr4[match(same_Kidney_chr4[i],Kidney_chr4$pig_pos),]
  b<-FM_Kidney_hum[match(intersect(a$index,FM_Kidney_hum$index),FM_Kidney_hum$index),]
  chr4_Kidney_hum<-rbind(chr4_Kidney_hum,b)
}
chr4_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr4_Kidney_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr4_Kidney_hum_one2one<-chr4_Kidney_hum[match(intersect(chr4_pig2hum_one2one,chr4_Kidney_hum$gene_id),chr4_Kidney_hum$gene_id),]

chr5_Kidney_overloci_one2one<-chr5_Kidney_overloci[match(intersect(chr5_Kidney_genes,one2one_pig),chr5_Kidney_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Kidney_chr5$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Kidney_chr5$hum_chr<-tmp$V1
Kidney_chr5$hum_loci<-tmp$V3
FM_Kidney_hum$index<-paste0(FM_Kidney_hum$chr,"_",FM_Kidney_hum$variant_pos)
Kidney_chr5$index<-paste0("chr",Kidney_chr5$hum_chr,"_",Kidney_chr5$hum_loci)
chr5_Kidney_hum<-NULL
for(i in 1:length(same_Kidney_chr5)){
  a<-Kidney_chr5[match(same_Kidney_chr5[i],Kidney_chr5$pig_pos),]
  b<-FM_Kidney_hum[match(intersect(a$index,FM_Kidney_hum$index),FM_Kidney_hum$index),]
  chr5_Kidney_hum<-rbind(chr5_Kidney_hum,b)
}
chr5_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr5_Kidney_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr5_Kidney_hum_one2one<-chr5_Kidney_hum[match(intersect(chr5_pig2hum_one2one,chr5_Kidney_hum$gene_id),chr5_Kidney_hum$gene_id),]

chr6_Kidney_overloci_one2one<-chr6_Kidney_overloci[match(intersect(chr6_Kidney_genes,one2one_pig),chr6_Kidney_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Kidney_chr6$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Kidney_chr6$hum_chr<-tmp$V1
Kidney_chr6$hum_loci<-tmp$V3
FM_Kidney_hum$index<-paste0(FM_Kidney_hum$chr,"_",FM_Kidney_hum$variant_pos)
Kidney_chr6$index<-paste0("chr",Kidney_chr6$hum_chr,"_",Kidney_chr6$hum_loci)
chr6_Kidney_hum<-NULL
for(i in 1:length(same_Kidney_chr6)){
  a<-Kidney_chr6[match(same_Kidney_chr6[i],Kidney_chr6$pig_pos),]
  b<-FM_Kidney_hum[match(intersect(a$index,FM_Kidney_hum$index),FM_Kidney_hum$index),]
  chr6_Kidney_hum<-rbind(chr6_Kidney_hum,b)
}
chr6_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr6_Kidney_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr6_Kidney_hum_one2one<-chr6_Kidney_hum[match(intersect(chr6_pig2hum_one2one,chr6_Kidney_hum$gene_id),chr6_Kidney_hum$gene_id),]

chr7_Kidney_overloci_one2one<-chr7_Kidney_overloci[match(intersect(chr7_Kidney_genes,one2one_pig),chr7_Kidney_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Kidney_chr7$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Kidney_chr7$hum_chr<-tmp$V1
Kidney_chr7$hum_loci<-tmp$V3
FM_Kidney_hum$index<-paste0(FM_Kidney_hum$chr,"_",FM_Kidney_hum$variant_pos)
Kidney_chr7$index<-paste0("chr",Kidney_chr7$hum_chr,"_",Kidney_chr7$hum_loci)
chr7_Kidney_hum<-NULL
for(i in 1:length(same_Kidney_chr7)){
  a<-Kidney_chr7[match(same_Kidney_chr7[i],Kidney_chr7$pig_pos),]
  b<-FM_Kidney_hum[match(intersect(a$index,FM_Kidney_hum$index),FM_Kidney_hum$index),]
  chr7_Kidney_hum<-rbind(chr7_Kidney_hum,b)
}
chr7_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr7_Kidney_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr7_Kidney_hum_one2one<-chr7_Kidney_hum[match(intersect(chr7_pig2hum_one2one,chr7_Kidney_hum$gene_id),chr7_Kidney_hum$gene_id),]

chr8_Kidney_overloci_one2one<-chr8_Kidney_overloci[match(intersect(chr8_Kidney_genes,one2one_pig),chr8_Kidney_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Kidney_chr8$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Kidney_chr8$hum_chr<-tmp$V1
Kidney_chr8$hum_loci<-tmp$V3
FM_Kidney_hum$index<-paste0(FM_Kidney_hum$chr,"_",FM_Kidney_hum$variant_pos)
Kidney_chr8$index<-paste0("chr",Kidney_chr8$hum_chr,"_",Kidney_chr8$hum_loci)
chr8_Kidney_hum<-NULL
for(i in 1:length(same_Kidney_chr8)){
  a<-Kidney_chr8[match(same_Kidney_chr8[i],Kidney_chr8$pig_pos),]
  b<-FM_Kidney_hum[match(intersect(a$index,FM_Kidney_hum$index),FM_Kidney_hum$index),]
  chr8_Kidney_hum<-rbind(chr8_Kidney_hum,b)
}
chr8_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr8_Kidney_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr8_Kidney_hum_one2one<-chr8_Kidney_hum[match(intersect(chr8_pig2hum_one2one,chr8_Kidney_hum$gene_id),chr8_Kidney_hum$gene_id),]

chr9_Kidney_overloci_one2one<-chr9_Kidney_overloci[match(intersect(chr9_Kidney_genes,one2one_pig),chr9_Kidney_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Kidney_chr9$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Kidney_chr9$hum_chr<-tmp$V1
Kidney_chr9$hum_loci<-tmp$V3
FM_Kidney_hum$index<-paste0(FM_Kidney_hum$chr,"_",FM_Kidney_hum$variant_pos)
Kidney_chr9$index<-paste0("chr",Kidney_chr9$hum_chr,"_",Kidney_chr9$hum_loci)
chr9_Kidney_hum<-NULL
for(i in 1:length(same_Kidney_chr9)){
  a<-Kidney_chr9[match(same_Kidney_chr9[i],Kidney_chr9$pig_pos),]
  b<-FM_Kidney_hum[match(intersect(a$index,FM_Kidney_hum$index),FM_Kidney_hum$index),]
  chr9_Kidney_hum<-rbind(chr9_Kidney_hum,b)
}
chr9_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr9_Kidney_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr9_Kidney_hum_one2one<-chr9_Kidney_hum[match(intersect(chr9_pig2hum_one2one,chr9_Kidney_hum$gene_id),chr9_Kidney_hum$gene_id),]

chr10_Kidney_overloci_one2one<-chr10_Kidney_overloci[match(intersect(chr10_Kidney_genes,one2one_pig),chr10_Kidney_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Kidney_chr10$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Kidney_chr10$hum_chr<-tmp$V1
Kidney_chr10$hum_loci<-tmp$V3
FM_Kidney_hum$index<-paste0(FM_Kidney_hum$chr,"_",FM_Kidney_hum$variant_pos)
Kidney_chr10$index<-paste0("chr",Kidney_chr10$hum_chr,"_",Kidney_chr10$hum_loci)
chr10_Kidney_hum<-NULL
for(i in 1:length(same_Kidney_chr10)){
  a<-Kidney_chr10[match(same_Kidney_chr10[i],Kidney_chr10$pig_pos),]
  b<-FM_Kidney_hum[match(intersect(a$index,FM_Kidney_hum$index),FM_Kidney_hum$index),]
  chr10_Kidney_hum<-rbind(chr10_Kidney_hum,b)
}
chr10_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr10_Kidney_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr10_Kidney_hum_one2one<-chr10_Kidney_hum[match(intersect(chr10_pig2hum_one2one,chr10_Kidney_hum$gene_id),chr10_Kidney_hum$gene_id),]

chr11_Kidney_overloci_one2one<-chr11_Kidney_overloci[match(intersect(chr11_Kidney_genes,one2one_pig),chr11_Kidney_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Kidney_chr11$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Kidney_chr11$hum_chr<-tmp$V1
Kidney_chr11$hum_loci<-tmp$V3
FM_Kidney_hum$index<-paste0(FM_Kidney_hum$chr,"_",FM_Kidney_hum$variant_pos)
Kidney_chr11$index<-paste0("chr",Kidney_chr11$hum_chr,"_",Kidney_chr11$hum_loci)
chr11_Kidney_hum<-NULL
for(i in 1:length(same_Kidney_chr11)){
  a<-Kidney_chr11[match(same_Kidney_chr11[i],Kidney_chr11$pig_pos),]
  b<-FM_Kidney_hum[match(intersect(a$index,FM_Kidney_hum$index),FM_Kidney_hum$index),]
  chr11_Kidney_hum<-rbind(chr11_Kidney_hum,b)
}
chr11_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr11_Kidney_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr11_Kidney_hum_one2one<-chr11_Kidney_hum[match(intersect(chr11_pig2hum_one2one,chr11_Kidney_hum$gene_id),chr11_Kidney_hum$gene_id),]

chr12_Kidney_overloci_one2one<-chr12_Kidney_overloci[match(intersect(chr12_Kidney_genes,one2one_pig),chr12_Kidney_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Kidney_chr12$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Kidney_chr12$hum_chr<-tmp$V1
Kidney_chr12$hum_loci<-tmp$V3
FM_Kidney_hum$index<-paste0(FM_Kidney_hum$chr,"_",FM_Kidney_hum$variant_pos)
Kidney_chr12$index<-paste0("chr",Kidney_chr12$hum_chr,"_",Kidney_chr12$hum_loci)
chr12_Kidney_hum<-NULL
for(i in 1:length(same_Kidney_chr12)){
  a<-Kidney_chr12[match(same_Kidney_chr12[i],Kidney_chr12$pig_pos),]
  b<-FM_Kidney_hum[match(intersect(a$index,FM_Kidney_hum$index),FM_Kidney_hum$index),]
  chr12_Kidney_hum<-rbind(chr12_Kidney_hum,b)
}
chr12_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr12_Kidney_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr12_Kidney_hum_one2one<-chr12_Kidney_hum[match(intersect(chr12_pig2hum_one2one,chr12_Kidney_hum$gene_id),chr12_Kidney_hum$gene_id),]

chr13_Kidney_overloci_one2one<-chr13_Kidney_overloci[match(intersect(chr13_Kidney_genes,one2one_pig),chr13_Kidney_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Kidney_chr13$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Kidney_chr13$hum_chr<-tmp$V1
Kidney_chr13$hum_loci<-tmp$V3
FM_Kidney_hum$index<-paste0(FM_Kidney_hum$chr,"_",FM_Kidney_hum$variant_pos)
Kidney_chr13$index<-paste0("chr",Kidney_chr13$hum_chr,"_",Kidney_chr13$hum_loci)
chr13_Kidney_hum<-NULL
for(i in 1:length(same_Kidney_chr13)){
  a<-Kidney_chr13[match(same_Kidney_chr13[i],Kidney_chr13$pig_pos),]
  b<-FM_Kidney_hum[match(intersect(a$index,FM_Kidney_hum$index),FM_Kidney_hum$index),]
  chr13_Kidney_hum<-rbind(chr13_Kidney_hum,b)
}
chr13_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr13_Kidney_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr13_Kidney_hum_one2one<-chr13_Kidney_hum[match(intersect(chr13_pig2hum_one2one,chr13_Kidney_hum$gene_id),chr13_Kidney_hum$gene_id),]

chr14_Kidney_overloci_one2one<-chr14_Kidney_overloci[match(intersect(chr14_Kidney_genes,one2one_pig),chr14_Kidney_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Kidney_chr14$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Kidney_chr14$hum_chr<-tmp$V1
Kidney_chr14$hum_loci<-tmp$V3
FM_Kidney_hum$index<-paste0(FM_Kidney_hum$chr,"_",FM_Kidney_hum$variant_pos)
Kidney_chr14$index<-paste0("chr",Kidney_chr14$hum_chr,"_",Kidney_chr14$hum_loci)
chr14_Kidney_hum<-NULL
for(i in 1:length(same_Kidney_chr14)){
  a<-Kidney_chr14[match(same_Kidney_chr14[i],Kidney_chr14$pig_pos),]
  b<-FM_Kidney_hum[match(intersect(a$index,FM_Kidney_hum$index),FM_Kidney_hum$index),]
  chr14_Kidney_hum<-rbind(chr14_Kidney_hum,b)
}
chr14_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr14_Kidney_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr14_Kidney_hum_one2one<-chr14_Kidney_hum[match(intersect(chr14_pig2hum_one2one,chr14_Kidney_hum$gene_id),chr14_Kidney_hum$gene_id),]

chr15_Kidney_overloci_one2one<-chr15_Kidney_overloci[match(intersect(chr15_Kidney_genes,one2one_pig),chr15_Kidney_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Kidney_chr15$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Kidney_chr15$hum_chr<-tmp$V1
Kidney_chr15$hum_loci<-tmp$V3
FM_Kidney_hum$index<-paste0(FM_Kidney_hum$chr,"_",FM_Kidney_hum$variant_pos)
Kidney_chr15$index<-paste0("chr",Kidney_chr15$hum_chr,"_",Kidney_chr15$hum_loci)
chr15_Kidney_hum<-NULL
for(i in 1:length(same_Kidney_chr15)){
  a<-Kidney_chr15[match(same_Kidney_chr15[i],Kidney_chr15$pig_pos),]
  b<-FM_Kidney_hum[match(intersect(a$index,FM_Kidney_hum$index),FM_Kidney_hum$index),]
  chr15_Kidney_hum<-rbind(chr15_Kidney_hum,b)
}
chr15_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr15_Kidney_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr15_Kidney_hum_one2one<-chr15_Kidney_hum[match(intersect(chr15_pig2hum_one2one,chr15_Kidney_hum$gene_id),chr15_Kidney_hum$gene_id),]

chr16_Kidney_overloci_one2one<-chr16_Kidney_overloci[match(intersect(chr16_Kidney_genes,one2one_pig),chr16_Kidney_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Kidney_chr16$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Kidney_chr16$hum_chr<-tmp$V1
Kidney_chr16$hum_loci<-tmp$V3
FM_Kidney_hum$index<-paste0(FM_Kidney_hum$chr,"_",FM_Kidney_hum$variant_pos)
Kidney_chr16$index<-paste0("chr",Kidney_chr16$hum_chr,"_",Kidney_chr16$hum_loci)
chr16_Kidney_hum<-NULL
for(i in 1:length(same_Kidney_chr16)){
  a<-Kidney_chr16[match(same_Kidney_chr16[i],Kidney_chr16$pig_pos),]
  b<-FM_Kidney_hum[match(intersect(a$index,FM_Kidney_hum$index),FM_Kidney_hum$index),]
  chr16_Kidney_hum<-rbind(chr16_Kidney_hum,b)
}
chr16_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr16_Kidney_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr16_Kidney_hum_one2one<-chr16_Kidney_hum[match(intersect(chr16_pig2hum_one2one,chr16_Kidney_hum$gene_id),chr16_Kidney_hum$gene_id),]

chr17_Kidney_overloci_one2one<-chr17_Kidney_overloci[match(intersect(chr17_Kidney_genes,one2one_pig),chr17_Kidney_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Kidney_chr17$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Kidney_chr17$hum_chr<-tmp$V1
Kidney_chr17$hum_loci<-tmp$V3
FM_Kidney_hum$index<-paste0(FM_Kidney_hum$chr,"_",FM_Kidney_hum$variant_pos)
Kidney_chr17$index<-paste0("chr",Kidney_chr17$hum_chr,"_",Kidney_chr17$hum_loci)
chr17_Kidney_hum<-NULL
for(i in 1:length(same_Kidney_chr17)){
  a<-Kidney_chr17[match(same_Kidney_chr17[i],Kidney_chr17$pig_pos),]
  b<-FM_Kidney_hum[match(intersect(a$index,FM_Kidney_hum$index),FM_Kidney_hum$index),]
  chr17_Kidney_hum<-rbind(chr17_Kidney_hum,b)
}
chr17_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr17_Kidney_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr17_Kidney_hum_one2one<-chr17_Kidney_hum[match(intersect(chr17_pig2hum_one2one,chr17_Kidney_hum$gene_id),chr17_Kidney_hum$gene_id),]

chr18_Kidney_overloci_one2one<-chr18_Kidney_overloci[match(intersect(chr18_Kidney_genes,one2one_pig),chr18_Kidney_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Kidney_chr18$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Kidney_chr18$hum_chr<-tmp$V1
Kidney_chr18$hum_loci<-tmp$V3
FM_Kidney_hum$index<-paste0(FM_Kidney_hum$chr,"_",FM_Kidney_hum$variant_pos)
Kidney_chr18$index<-paste0("chr",Kidney_chr18$hum_chr,"_",Kidney_chr18$hum_loci)
chr18_Kidney_hum<-NULL
for(i in 1:length(same_Kidney_chr18)){
  a<-Kidney_chr18[match(same_Kidney_chr18[i],Kidney_chr18$pig_pos),]
  b<-FM_Kidney_hum[match(intersect(a$index,FM_Kidney_hum$index),FM_Kidney_hum$index),]
  chr18_Kidney_hum<-rbind(chr18_Kidney_hum,b)
}
chr18_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr18_Kidney_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr18_Kidney_hum_one2one<-chr18_Kidney_hum[match(intersect(chr18_pig2hum_one2one,chr18_Kidney_hum$gene_id),chr18_Kidney_hum$gene_id),]

Kidney_one2one_SNP_hum<-rbind(chr1_Kidney_hum_one2one,chr2_Kidney_hum_one2one,chr3_Kidney_hum_one2one,chr4_Kidney_hum_one2one,chr5_Kidney_hum_one2one,
                              chr6_Kidney_hum_one2one,chr7_Kidney_hum_one2one,chr8_Kidney_hum_one2one,chr9_Kidney_hum_one2one,chr10_Kidney_hum_one2one,
                              chr11_Kidney_hum_one2one,chr12_Kidney_hum_one2one,chr13_Kidney_hum_one2one,chr14_Kidney_hum_one2one,chr15_Kidney_hum_one2one,
                              chr16_Kidney_hum_one2one,chr17_Kidney_hum_one2one,chr18_Kidney_hum_one2one)

chr1_Kidney_pig_one2one<-chr1_Kidney_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr1_Kidney_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr1_Kidney_overloci_one2one$phenotype_id),1:9]
chr2_Kidney_pig_one2one<-chr2_Kidney_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr2_Kidney_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr2_Kidney_overloci_one2one$phenotype_id),1:9]
chr3_Kidney_pig_one2one<-chr3_Kidney_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr3_Kidney_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr3_Kidney_overloci_one2one$phenotype_id),1:9]
chr4_Kidney_pig_one2one<-chr4_Kidney_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr4_Kidney_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr4_Kidney_overloci_one2one$phenotype_id),1:9]
chr5_Kidney_pig_one2one<-chr5_Kidney_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr5_Kidney_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr5_Kidney_overloci_one2one$phenotype_id),1:9]
chr6_Kidney_pig_one2one<-chr6_Kidney_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr6_Kidney_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr6_Kidney_overloci_one2one$phenotype_id),1:9]
chr7_Kidney_pig_one2one<-chr7_Kidney_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr7_Kidney_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr7_Kidney_overloci_one2one$phenotype_id),1:9]
chr8_Kidney_pig_one2one<-chr8_Kidney_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr8_Kidney_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr8_Kidney_overloci_one2one$phenotype_id),1:9]
chr9_Kidney_pig_one2one<-chr9_Kidney_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr9_Kidney_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr9_Kidney_overloci_one2one$phenotype_id),1:9]
chr10_Kidney_pig_one2one<-chr10_Kidney_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr10_Kidney_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr10_Kidney_overloci_one2one$phenotype_id),1:9]
chr11_Kidney_pig_one2one<-chr11_Kidney_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr11_Kidney_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr11_Kidney_overloci_one2one$phenotype_id),1:9]
chr12_Kidney_pig_one2one<-chr12_Kidney_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr12_Kidney_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr12_Kidney_overloci_one2one$phenotype_id),1:9]
chr13_Kidney_pig_one2one<-chr13_Kidney_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr13_Kidney_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr13_Kidney_overloci_one2one$phenotype_id),1:9]
chr14_Kidney_pig_one2one<-chr14_Kidney_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr14_Kidney_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr14_Kidney_overloci_one2one$phenotype_id),1:9]
chr15_Kidney_pig_one2one<-chr15_Kidney_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr15_Kidney_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr15_Kidney_overloci_one2one$phenotype_id),1:9]
chr16_Kidney_pig_one2one<-chr16_Kidney_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr16_Kidney_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr16_Kidney_overloci_one2one$phenotype_id),1:9]
chr17_Kidney_pig_one2one<-chr17_Kidney_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr17_Kidney_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr17_Kidney_overloci_one2one$phenotype_id),1:9]
chr18_Kidney_pig_one2one<-chr18_Kidney_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr18_Kidney_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr18_Kidney_overloci_one2one$phenotype_id),1:9]

Kidney_one2one_SNP_pig<-rbind(chr1_Kidney_pig_one2one,chr2_Kidney_pig_one2one,chr3_Kidney_pig_one2one,chr4_Kidney_pig_one2one,chr5_Kidney_pig_one2one,
                              chr6_Kidney_pig_one2one,chr7_Kidney_pig_one2one,chr8_Kidney_pig_one2one,chr9_Kidney_pig_one2one,chr10_Kidney_pig_one2one,
                              chr11_Kidney_pig_one2one,chr12_Kidney_pig_one2one,chr13_Kidney_pig_one2one,chr14_Kidney_pig_one2one,chr15_Kidney_pig_one2one,
                              chr16_Kidney_pig_one2one,chr17_Kidney_pig_one2one,chr18_Kidney_pig_one2one)

Kidney_SNP_sum<-array(NA,dim=c(nrow(Kidney_one2one_SNP_hum),2))
colnames(Kidney_SNP_sum)<-c("Human","Pig")
Kidney_SNP_sum<-as.data.frame(Kidney_SNP_sum)
Kidney_SNP_sum$Human<-Kidney_one2one_SNP_hum$slope / Kidney_one2one_SNP_hum$slope_se
Kidney_SNP_sum$Pig<-Kidney_one2one_SNP_pig$slope / Kidney_one2one_SNP_pig$slope_se
cor<-cor(abs(Kidney_SNP_sum$Human),abs(Kidney_SNP_sum$Pig))
p_val<-t.test(abs(Kidney_SNP_sum$Human),abs(Kidney_SNP_sum$Pig))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/Kidney_SNP_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(Kidney_SNP_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The number of eGenes expressed in common tissues",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

save(Kidney_one2one_SNP_hum,Kidney_one2one_SNP_pig,file="/Users/baizhonghao/Downloads/human-pig GTEx/Kidney_SNP.Rdata")

#Liver_SNP_overlaploci#
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

eqtl_Liver_chr1<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Liver/Liver.cis_qtl_pairs.1.txt"))
eqtl_Liver_chr2<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Liver/Liver.cis_qtl_pairs.2.txt"))
eqtl_Liver_chr3<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Liver/Liver.cis_qtl_pairs.3.txt"))
eqtl_Liver_chr4<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Liver/Liver.cis_qtl_pairs.4.txt"))
eqtl_Liver_chr5<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Liver/Liver.cis_qtl_pairs.5.txt"))
eqtl_Liver_chr6<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Liver/Liver.cis_qtl_pairs.6.txt"))
eqtl_Liver_chr7<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Liver/Liver.cis_qtl_pairs.7.txt"))
eqtl_Liver_chr8<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Liver/Liver.cis_qtl_pairs.8.txt"))
eqtl_Liver_chr9<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Liver/Liver.cis_qtl_pairs.9.txt"))
eqtl_Liver_chr10<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Liver/Liver.cis_qtl_pairs.10.txt"))
eqtl_Liver_chr11<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Liver/Liver.cis_qtl_pairs.11.txt"))
eqtl_Liver_chr12<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Liver/Liver.cis_qtl_pairs.12.txt"))
eqtl_Liver_chr13<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Liver/Liver.cis_qtl_pairs.13.txt"))
eqtl_Liver_chr14<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Liver/Liver.cis_qtl_pairs.14.txt"))
eqtl_Liver_chr15<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Liver/Liver.cis_qtl_pairs.15.txt"))
eqtl_Liver_chr16<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Liver/Liver.cis_qtl_pairs.16.txt"))
eqtl_Liver_chr17<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Liver/Liver.cis_qtl_pairs.17.txt"))
eqtl_Liver_chr18<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Liver/Liver.cis_qtl_pairs.18.txt"))

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver_chr1$variant_id,split="_"))))
eqtl_Liver_chr1$chr<-tmp$V1
eqtl_Liver_chr1$loci<-tmp$V2

eqtl_Liver_chr1$index<-paste0(eqtl_Liver_chr1$chr,"-",eqtl_Liver_chr1$loci)
Liver_chr1$index<-paste0(Liver_chr1$chr,"-",Liver_chr1$pig_pos)
same_Liver_chr1<-intersect(Liver_chr1$pig_pos,eqtl_Liver_chr1$loci)


tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver_chr2$variant_id,split="_"))))
eqtl_Liver_chr2$chr<-tmp$V1
eqtl_Liver_chr2$loci<-tmp$V2

eqtl_Liver_chr2$index<-paste0(eqtl_Liver_chr2$chr,"-",eqtl_Liver_chr2$loci)
Liver_chr2$index<-paste0(Liver_chr2$chr,"-",Liver_chr2$pig_pos)
same_Liver_chr2<-intersect(Liver_chr2$pig_pos,eqtl_Liver_chr2$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver_chr3$variant_id,split="_"))))
eqtl_Liver_chr3$chr<-tmp$V1
eqtl_Liver_chr3$loci<-tmp$V2

eqtl_Liver_chr3$index<-paste0(eqtl_Liver_chr3$chr,"-",eqtl_Liver_chr3$loci)
Liver_chr3$index<-paste0(Liver_chr3$chr,"-",Liver_chr3$pig_pos)
same_Liver_chr3<-intersect(Liver_chr3$pig_pos,eqtl_Liver_chr3$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver_chr4$variant_id,split="_"))))
eqtl_Liver_chr4$chr<-tmp$V1
eqtl_Liver_chr4$loci<-tmp$V2

eqtl_Liver_chr4$index<-paste0(eqtl_Liver_chr4$chr,"-",eqtl_Liver_chr4$loci)
Liver_chr4$index<-paste0(Liver_chr4$chr,"-",Liver_chr4$pig_pos)
same_Liver_chr4<-intersect(Liver_chr4$pig_pos,eqtl_Liver_chr4$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver_chr5$variant_id,split="_"))))
eqtl_Liver_chr5$loci<-tmp$V2
same_Liver_chr5<-intersect(Liver_chr5$pig_pos,eqtl_Liver_chr5$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver_chr6$variant_id,split="_"))))
eqtl_Liver_chr6$loci<-tmp$V2
same_Liver_chr6<-intersect(Liver_chr6$pig_pos,eqtl_Liver_chr6$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver_chr7$variant_id,split="_"))))
eqtl_Liver_chr7$loci<-tmp$V2
same_Liver_chr7<-intersect(Liver_chr7$pig_pos,eqtl_Liver_chr7$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver_chr8$variant_id,split="_"))))
eqtl_Liver_chr8$loci<-tmp$V2
same_Liver_chr8<-intersect(Liver_chr8$pig_pos,eqtl_Liver_chr8$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver_chr9$variant_id,split="_"))))
eqtl_Liver_chr9$loci<-tmp$V2
same_Liver_chr9<-intersect(Liver_chr9$pig_pos,eqtl_Liver_chr9$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver_chr10$variant_id,split="_"))))
eqtl_Liver_chr10$loci<-tmp$V2
same_Liver_chr10<-intersect(Liver_chr10$pig_pos,eqtl_Liver_chr10$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver_chr11$variant_id,split="_"))))
eqtl_Liver_chr11$loci<-tmp$V2
same_Liver_chr11<-intersect(Liver_chr11$pig_pos,eqtl_Liver_chr11$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver_chr12$variant_id,split="_"))))
eqtl_Liver_chr12$loci<-tmp$V2
same_Liver_chr12<-intersect(Liver_chr12$pig_pos,eqtl_Liver_chr12$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver_chr13$variant_id,split="_"))))
eqtl_Liver_chr13$loci<-tmp$V2
same_Liver_chr13<-intersect(Liver_chr13$pig_pos,eqtl_Liver_chr13$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver_chr14$variant_id,split="_"))))
eqtl_Liver_chr14$loci<-tmp$V2
same_Liver_chr14<-intersect(Liver_chr14$pig_pos,eqtl_Liver_chr14$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver_chr15$variant_id,split="_"))))
eqtl_Liver_chr15$loci<-tmp$V2
same_Liver_chr15<-intersect(Liver_chr15$pig_pos,eqtl_Liver_chr15$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver_chr16$variant_id,split="_"))))
eqtl_Liver_chr16$loci<-tmp$V2
same_Liver_chr16<-intersect(Liver_chr16$pig_pos,eqtl_Liver_chr16$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver_chr17$variant_id,split="_"))))
eqtl_Liver_chr17$loci<-tmp$V2
same_Liver_chr17<-intersect(Liver_chr17$pig_pos,eqtl_Liver_chr17$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Liver_chr18$variant_id,split="_"))))
eqtl_Liver_chr18$loci<-tmp$V2
same_Liver_chr18<-intersect(Liver_chr18$pig_pos,eqtl_Liver_chr18$loci)

chr1_Liver_genes<-NULL
chr1_Liver_overloci<-NULL
if(length(same_Liver_chr1!=0)){
  for(i in 1:length(same_Liver_chr1)){
    a<-eqtl_Liver_chr1$phenotype_id[grep(same_Liver_chr1[i],eqtl_Liver_chr1$loci)]
    b<-eqtl_Liver_chr1[grep(same_Liver_chr1[i],eqtl_Liver_chr1$loci),]
    chr1_Liver_genes<-c(chr1_Liver_genes,a)
    chr1_Liver_overloci<-rbind(chr1_Liver_overloci,b)
  }
}

chr2_Liver_genes<-NULL
chr2_Liver_overloci<-NULL
if(length(same_Liver_chr2!=0)){
  for(i in 1:length(same_Liver_chr2)){
    a<-eqtl_Liver_chr2$phenotype_id[grep(same_Liver_chr2[i],eqtl_Liver_chr2$loci)]
    b<-eqtl_Liver_chr2[grep(same_Liver_chr2[i],eqtl_Liver_chr2$loci),]
    chr2_Liver_genes<-c(chr2_Liver_genes,a)
    chr2_Liver_overloci<-rbind(chr2_Liver_overloci,b)
  }
}

chr3_Liver_genes<-NULL
chr3_Liver_overloci<-NULL
if(length(same_Liver_chr3!=0)){
  for(i in 1:length(same_Liver_chr3)){
    a<-eqtl_Liver_chr3$phenotype_id[grep(same_Liver_chr3[i],eqtl_Liver_chr3$loci)]
    b<-eqtl_Liver_chr3[grep(same_Liver_chr3[i],eqtl_Liver_chr3$loci),]
    chr3_Liver_genes<-c(chr3_Liver_genes,a)
    chr3_Liver_overloci<-rbind(chr3_Liver_overloci,b)
  }
}

chr4_Liver_genes<-NULL
chr4_Liver_overloci<-NULL
if(length(same_Liver_chr4!=0)){
  for(i in 1:length(same_Liver_chr4)){
    a<-eqtl_Liver_chr4$phenotype_id[grep(same_Liver_chr4[i],eqtl_Liver_chr4$loci)]
    b<-eqtl_Liver_chr4[grep(same_Liver_chr4[i],eqtl_Liver_chr4$loci),]
    chr4_Liver_genes<-c(chr4_Liver_genes,a)
    chr4_Liver_overloci<-rbind(chr4_Liver_overloci,b)
  }
}

chr5_Liver_genes<-NULL
chr5_Liver_overloci<-NULL
if(length(same_Liver_chr5!=0)){
  for(i in 1:length(same_Liver_chr5)){
    a<-eqtl_Liver_chr5$phenotype_id[grep(same_Liver_chr5[i],eqtl_Liver_chr5$loci)]
    b<-eqtl_Liver_chr5[grep(same_Liver_chr5[i],eqtl_Liver_chr5$loci),]
    chr5_Liver_genes<-c(chr5_Liver_genes,a)
    chr5_Liver_overloci<-rbind(chr5_Liver_overloci,b)
  }
}

chr6_Liver_genes<-NULL
chr6_Liver_overloci<-NULL
if(length(same_Liver_chr6!=0)){
  for(i in 1:length(same_Liver_chr6)){
    a<-eqtl_Liver_chr6$phenotype_id[grep(same_Liver_chr6[i],eqtl_Liver_chr6$loci)]
    b<-eqtl_Liver_chr6[grep(same_Liver_chr6[i],eqtl_Liver_chr6$loci),]
    chr6_Liver_genes<-c(chr6_Liver_genes,a)
    chr6_Liver_overloci<-rbind(chr6_Liver_overloci,b)
  }
}

chr7_Liver_genes<-NULL
chr7_Liver_overloci<-NULL
if(length(same_Liver_chr7!=0)){
  for(i in 1:length(same_Liver_chr7)){
    a<-eqtl_Liver_chr7$phenotype_id[grep(same_Liver_chr7[i],eqtl_Liver_chr7$loci)]
    b<-eqtl_Liver_chr7[grep(same_Liver_chr7[i],eqtl_Liver_chr7$loci),]
    chr7_Liver_genes<-c(chr7_Liver_genes,a)
    chr7_Liver_overloci<-rbind(chr7_Liver_overloci,b)
  }
}

chr8_Liver_genes<-NULL
chr8_Liver_overloci<-NULL
if(length(same_Liver_chr8!=0)){
  for(i in 1:length(same_Liver_chr8)){
    a<-eqtl_Liver_chr8$phenotype_id[grep(same_Liver_chr8[i],eqtl_Liver_chr8$loci)]
    b<-eqtl_Liver_chr8[grep(same_Liver_chr8[i],eqtl_Liver_chr8$loci),]
    chr8_Liver_genes<-c(chr8_Liver_genes,a)
    chr8_Liver_overloci<-rbind(chr8_Liver_overloci,b)
  }
}

if(length(same_Liver_chr9!=0)){
  chr9_Liver_genes<-NULL
  chr9_Liver_overloci<-NULL
  for(i in 1:length(same_Liver_chr9)){
    a<-eqtl_Liver_chr9$phenotype_id[grep(same_Liver_chr9[i],eqtl_Liver_chr9$loci)]
    b<-eqtl_Liver_chr9[grep(same_Liver_chr9[i],eqtl_Liver_chr9$loci),]
    chr9_Liver_genes<-c(chr9_Liver_genes,a)
    chr9_Liver_overloci<-rbind(chr9_Liver_overloci,b)
  }
}

chr10_Liver_genes<-NULL
chr10_Liver_overloci<-NULL
if(length(same_Liver_chr10!=0)){
  for(i in 1:length(same_Liver_chr10)){
    a<-eqtl_Liver_chr10$phenotype_id[grep(same_Liver_chr10[i],eqtl_Liver_chr10$loci)]
    b<-eqtl_Liver_chr10[grep(same_Liver_chr10[i],eqtl_Liver_chr10$loci),]
    chr10_Liver_genes<-c(chr10_Liver_genes,a)
    chr10_Liver_overloci<-rbind(chr10_Liver_overloci,b)
  }
}

chr11_Liver_genes<-NULL
chr11_Liver_overloci<-NULL
if(length(same_Liver_chr11!=0)){
  for(i in 1:length(same_Liver_chr11)){
    a<-eqtl_Liver_chr11$phenotype_id[grep(same_Liver_chr11[i],eqtl_Liver_chr11$loci)]
    b<-eqtl_Liver_chr11[grep(same_Liver_chr11[i],eqtl_Liver_chr11$loci),]
    chr11_Liver_genes<-c(chr11_Liver_genes,a)
    chr11_Liver_overloci<-rbind(chr11_Liver_overloci,b)
  }
}

chr12_Liver_genes<-NULL
chr12_Liver_overloci<-NULL
if(length(same_Liver_chr12!=0)){
  for(i in 1:length(same_Liver_chr12)){
    a<-eqtl_Liver_chr12$phenotype_id[grep(same_Liver_chr12[i],eqtl_Liver_chr12$loci)]
    b<-eqtl_Liver_chr12[grep(same_Liver_chr12[i],eqtl_Liver_chr12$loci),]
    chr12_Liver_genes<-c(chr12_Liver_genes,a)
    chr12_Liver_overloci<-rbind(chr12_Liver_overloci,b)
  }
}

chr13_Liver_genes<-NULL
chr13_Liver_overloci<-NULL
if(length(same_Liver_chr13!=0)){
  for(i in 1:length(same_Liver_chr13)){
    a<-eqtl_Liver_chr13$phenotype_id[grep(same_Liver_chr13[i],eqtl_Liver_chr13$loci)]
    b<-eqtl_Liver_chr13[grep(same_Liver_chr13[i],eqtl_Liver_chr13$loci),]
    chr13_Liver_genes<-c(chr13_Liver_genes,a)
    chr13_Liver_overloci<-rbind(chr13_Liver_overloci,b)
  }
}

chr14_Liver_genes<-NULL
chr14_Liver_overloci<-NULL
if(length(same_Liver_chr14!=0)){
  for(i in 1:length(same_Liver_chr14)){
    a<-eqtl_Liver_chr14$phenotype_id[grep(same_Liver_chr14[i],eqtl_Liver_chr14$loci)]
    b<-eqtl_Liver_chr14[grep(same_Liver_chr14[i],eqtl_Liver_chr14$loci),]
    chr14_Liver_genes<-c(chr14_Liver_genes,a)
    chr14_Liver_overloci<-rbind(chr14_Liver_overloci,b)
  }
}

chr15_Liver_genes<-NULL
chr15_Liver_overloci<-NULL
if(length(same_Liver_chr15!=0)){
  for(i in 1:length(same_Liver_chr15)){
    a<-eqtl_Liver_chr15$phenotype_id[grep(same_Liver_chr15[i],eqtl_Liver_chr15$loci)]
    b<-eqtl_Liver_chr15[grep(same_Liver_chr15[i],eqtl_Liver_chr15$loci),]
    chr15_Liver_genes<-c(chr15_Liver_genes,a)
    chr15_Liver_overloci<-rbind(chr15_Liver_overloci,b)
  }
}

chr16_Liver_genes<-NULL
chr16_Liver_overloci<-NULL
if(length(same_Liver_chr16!=0)){
  for(i in 1:length(same_Liver_chr16)){
    a<-eqtl_Liver_chr16$phenotype_id[grep(same_Liver_chr16[i],eqtl_Liver_chr16$loci)]
    b<-eqtl_Liver_chr16[grep(same_Liver_chr16[i],eqtl_Liver_chr16$loci),]
    chr16_Liver_genes<-c(chr16_Liver_genes,a)
    chr16_Liver_overloci<-rbind(chr16_Liver_overloci,b)
  }
}

chr17_Liver_genes<-NULL
chr17_Liver_overloci<-NULL
if(length(same_Liver_chr17!=0)){
  for(i in 1:length(same_Liver_chr17)){
    a<-eqtl_Liver_chr17$phenotype_id[grep(same_Liver_chr17[i],eqtl_Liver_chr17$loci)]
    b<-eqtl_Liver_chr17[grep(same_Liver_chr17[i],eqtl_Liver_chr17$loci),]
    chr17_Liver_genes<-c(chr17_Liver_genes,a)
    chr17_Liver_overloci<-rbind(chr17_Liver_overloci,b)
  }
}

chr18_Liver_genes<-NULL
chr18_Liver_overloci<-NULL
if(length(same_Liver_chr18!=0)){
  for(i in 1:length(same_Liver_chr18)){
    a<-eqtl_Liver_chr18$phenotype_id[grep(same_Liver_chr18[i],eqtl_Liver_chr18$loci)]
    b<-eqtl_Liver_chr18[grep(same_Liver_chr18[i],eqtl_Liver_chr18$loci),]
    chr18_Liver_genes<-c(chr18_Liver_genes,a)
    chr18_Liver_overloci<-rbind(chr18_Liver_overloci,b)
  }
}
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
one2one_pig<-annotation$`Pig gene stable ID`[which(annotation$`Pig homology type`=="ortholog_one2one")]

chr1_Liver_overloci_one2one<-chr1_Liver_overloci[match(intersect(chr1_Liver_genes,one2one_pig),chr1_Liver_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Liver_chr1$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Liver_chr1$hum_chr<-tmp$V1
Liver_chr1$hum_loci<-tmp$V3
FM_Liver_hum$index<-paste0(FM_Liver_hum$chr,"_",FM_Liver_hum$variant_pos)
Liver_chr1$index<-paste0("chr",Liver_chr1$hum_chr,"_",Liver_chr1$hum_loci)
chr1_Liver_hum<-NULL
for(i in 1:length(same_Liver_chr1)){
  a<-Liver_chr1[match(same_Liver_chr1[i],Liver_chr1$pig_pos),]
  b<-FM_Liver_hum[match(intersect(a$index,FM_Liver_hum$index),FM_Liver_hum$index),]
  chr1_Liver_hum<-rbind(chr1_Liver_hum,b)
}
chr1_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr1_Liver_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr1_Liver_hum_one2one<-chr1_Liver_hum[match(intersect(chr1_pig2hum_one2one,chr1_Liver_hum$gene_id),chr1_Liver_hum$gene_id),]

chr2_Liver_overloci_one2one<-chr2_Liver_overloci[match(intersect(chr2_Liver_genes,one2one_pig),chr2_Liver_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Liver_chr2$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Liver_chr2$hum_chr<-tmp$V1
Liver_chr2$hum_loci<-tmp$V3
FM_Liver_hum$index<-paste0(FM_Liver_hum$chr,"_",FM_Liver_hum$variant_pos)
Liver_chr2$index<-paste0("chr",Liver_chr2$hum_chr,"_",Liver_chr2$hum_loci)
chr2_Liver_hum<-NULL
for(i in 1:length(same_Liver_chr2)){
  a<-Liver_chr2[match(same_Liver_chr2[i],Liver_chr2$pig_pos),]
  b<-FM_Liver_hum[match(intersect(a$index,FM_Liver_hum$index),FM_Liver_hum$index),]
  chr2_Liver_hum<-rbind(chr2_Liver_hum,b)
}
chr2_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr2_Liver_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr2_Liver_hum_one2one<-chr2_Liver_hum[match(intersect(chr2_pig2hum_one2one,chr2_Liver_hum$gene_id),chr2_Liver_hum$gene_id),]

chr3_Liver_overloci_one2one<-chr3_Liver_overloci[match(intersect(chr3_Liver_genes,one2one_pig),chr3_Liver_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Liver_chr3$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Liver_chr3$hum_chr<-tmp$V1
Liver_chr3$hum_loci<-tmp$V3
FM_Liver_hum$index<-paste0(FM_Liver_hum$chr,"_",FM_Liver_hum$variant_pos)
Liver_chr3$index<-paste0("chr",Liver_chr3$hum_chr,"_",Liver_chr3$hum_loci)
chr3_Liver_hum<-NULL
for(i in 1:length(same_Liver_chr3)){
  a<-Liver_chr3[match(same_Liver_chr3[i],Liver_chr3$pig_pos),]
  b<-FM_Liver_hum[match(intersect(a$index,FM_Liver_hum$index),FM_Liver_hum$index),]
  chr3_Liver_hum<-rbind(chr3_Liver_hum,b)
}
chr3_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr3_Liver_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr3_Liver_hum_one2one<-chr3_Liver_hum[match(intersect(chr3_pig2hum_one2one,chr3_Liver_hum$gene_id),chr3_Liver_hum$gene_id),]

chr4_Liver_overloci_one2one<-chr4_Liver_overloci[match(intersect(chr4_Liver_genes,one2one_pig),chr4_Liver_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Liver_chr4$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Liver_chr4$hum_chr<-tmp$V1
Liver_chr4$hum_loci<-tmp$V3
FM_Liver_hum$index<-paste0(FM_Liver_hum$chr,"_",FM_Liver_hum$variant_pos)
Liver_chr4$index<-paste0("chr",Liver_chr4$hum_chr,"_",Liver_chr4$hum_loci)
chr4_Liver_hum<-NULL
for(i in 1:length(same_Liver_chr4)){
  a<-Liver_chr4[match(same_Liver_chr4[i],Liver_chr4$pig_pos),]
  b<-FM_Liver_hum[match(intersect(a$index,FM_Liver_hum$index),FM_Liver_hum$index),]
  chr4_Liver_hum<-rbind(chr4_Liver_hum,b)
}
chr4_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr4_Liver_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr4_Liver_hum_one2one<-chr4_Liver_hum[match(intersect(chr4_pig2hum_one2one,chr4_Liver_hum$gene_id),chr4_Liver_hum$gene_id),]

chr5_Liver_overloci_one2one<-chr5_Liver_overloci[match(intersect(chr5_Liver_genes,one2one_pig),chr5_Liver_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Liver_chr5$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Liver_chr5$hum_chr<-tmp$V1
Liver_chr5$hum_loci<-tmp$V3
FM_Liver_hum$index<-paste0(FM_Liver_hum$chr,"_",FM_Liver_hum$variant_pos)
Liver_chr5$index<-paste0("chr",Liver_chr5$hum_chr,"_",Liver_chr5$hum_loci)
chr5_Liver_hum<-NULL
for(i in 1:length(same_Liver_chr5)){
  a<-Liver_chr5[match(same_Liver_chr5[i],Liver_chr5$pig_pos),]
  b<-FM_Liver_hum[match(intersect(a$index,FM_Liver_hum$index),FM_Liver_hum$index),]
  chr5_Liver_hum<-rbind(chr5_Liver_hum,b)
}
chr5_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr5_Liver_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr5_Liver_hum_one2one<-chr5_Liver_hum[match(intersect(chr5_pig2hum_one2one,chr5_Liver_hum$gene_id),chr5_Liver_hum$gene_id),]

chr6_Liver_overloci_one2one<-chr6_Liver_overloci[match(intersect(chr6_Liver_genes,one2one_pig),chr6_Liver_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Liver_chr6$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Liver_chr6$hum_chr<-tmp$V1
Liver_chr6$hum_loci<-tmp$V3
FM_Liver_hum$index<-paste0(FM_Liver_hum$chr,"_",FM_Liver_hum$variant_pos)
Liver_chr6$index<-paste0("chr",Liver_chr6$hum_chr,"_",Liver_chr6$hum_loci)
chr6_Liver_hum<-NULL
for(i in 1:length(same_Liver_chr6)){
  a<-Liver_chr6[match(same_Liver_chr6[i],Liver_chr6$pig_pos),]
  b<-FM_Liver_hum[match(intersect(a$index,FM_Liver_hum$index),FM_Liver_hum$index),]
  chr6_Liver_hum<-rbind(chr6_Liver_hum,b)
}
chr6_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr6_Liver_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr6_Liver_hum_one2one<-chr6_Liver_hum[match(intersect(chr6_pig2hum_one2one,chr6_Liver_hum$gene_id),chr6_Liver_hum$gene_id),]

chr7_Liver_overloci_one2one<-chr7_Liver_overloci[match(intersect(chr7_Liver_genes,one2one_pig),chr7_Liver_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Liver_chr7$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Liver_chr7$hum_chr<-tmp$V1
Liver_chr7$hum_loci<-tmp$V3
FM_Liver_hum$index<-paste0(FM_Liver_hum$chr,"_",FM_Liver_hum$variant_pos)
Liver_chr7$index<-paste0("chr",Liver_chr7$hum_chr,"_",Liver_chr7$hum_loci)
chr7_Liver_hum<-NULL
for(i in 1:length(same_Liver_chr7)){
  a<-Liver_chr7[match(same_Liver_chr7[i],Liver_chr7$pig_pos),]
  b<-FM_Liver_hum[match(intersect(a$index,FM_Liver_hum$index),FM_Liver_hum$index),]
  chr7_Liver_hum<-rbind(chr7_Liver_hum,b)
}
chr7_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr7_Liver_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr7_Liver_hum_one2one<-chr7_Liver_hum[match(intersect(chr7_pig2hum_one2one,chr7_Liver_hum$gene_id),chr7_Liver_hum$gene_id),]

chr8_Liver_overloci_one2one<-chr8_Liver_overloci[match(intersect(chr8_Liver_genes,one2one_pig),chr8_Liver_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Liver_chr8$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Liver_chr8$hum_chr<-tmp$V1
Liver_chr8$hum_loci<-tmp$V3
FM_Liver_hum$index<-paste0(FM_Liver_hum$chr,"_",FM_Liver_hum$variant_pos)
Liver_chr8$index<-paste0("chr",Liver_chr8$hum_chr,"_",Liver_chr8$hum_loci)
chr8_Liver_hum<-NULL
for(i in 1:length(same_Liver_chr8)){
  a<-Liver_chr8[match(same_Liver_chr8[i],Liver_chr8$pig_pos),]
  b<-FM_Liver_hum[match(intersect(a$index,FM_Liver_hum$index),FM_Liver_hum$index),]
  chr8_Liver_hum<-rbind(chr8_Liver_hum,b)
}
chr8_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr8_Liver_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr8_Liver_hum_one2one<-chr8_Liver_hum[match(intersect(chr8_pig2hum_one2one,chr8_Liver_hum$gene_id),chr8_Liver_hum$gene_id),]

chr9_Liver_overloci_one2one<-chr9_Liver_overloci[match(intersect(chr9_Liver_genes,one2one_pig),chr9_Liver_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Liver_chr9$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Liver_chr9$hum_chr<-tmp$V1
Liver_chr9$hum_loci<-tmp$V3
FM_Liver_hum$index<-paste0(FM_Liver_hum$chr,"_",FM_Liver_hum$variant_pos)
Liver_chr9$index<-paste0("chr",Liver_chr9$hum_chr,"_",Liver_chr9$hum_loci)
chr9_Liver_hum<-NULL
for(i in 1:length(same_Liver_chr9)){
  a<-Liver_chr9[match(same_Liver_chr9[i],Liver_chr9$pig_pos),]
  b<-FM_Liver_hum[match(intersect(a$index,FM_Liver_hum$index),FM_Liver_hum$index),]
  chr9_Liver_hum<-rbind(chr9_Liver_hum,b)
}
chr9_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr9_Liver_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr9_Liver_hum_one2one<-chr9_Liver_hum[match(intersect(chr9_pig2hum_one2one,chr9_Liver_hum$gene_id),chr9_Liver_hum$gene_id),]

chr10_Liver_overloci_one2one<-chr10_Liver_overloci[match(intersect(chr10_Liver_genes,one2one_pig),chr10_Liver_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Liver_chr10$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Liver_chr10$hum_chr<-tmp$V1
Liver_chr10$hum_loci<-tmp$V3
FM_Liver_hum$index<-paste0(FM_Liver_hum$chr,"_",FM_Liver_hum$variant_pos)
Liver_chr10$index<-paste0("chr",Liver_chr10$hum_chr,"_",Liver_chr10$hum_loci)
chr10_Liver_hum<-NULL
for(i in 1:length(same_Liver_chr10)){
  a<-Liver_chr10[match(same_Liver_chr10[i],Liver_chr10$pig_pos),]
  b<-FM_Liver_hum[match(intersect(a$index,FM_Liver_hum$index),FM_Liver_hum$index),]
  chr10_Liver_hum<-rbind(chr10_Liver_hum,b)
}
chr10_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr10_Liver_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr10_Liver_hum_one2one<-chr10_Liver_hum[match(intersect(chr10_pig2hum_one2one,chr10_Liver_hum$gene_id),chr10_Liver_hum$gene_id),]

chr11_Liver_overloci_one2one<-chr11_Liver_overloci[match(intersect(chr11_Liver_genes,one2one_pig),chr11_Liver_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Liver_chr11$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Liver_chr11$hum_chr<-tmp$V1
Liver_chr11$hum_loci<-tmp$V3
FM_Liver_hum$index<-paste0(FM_Liver_hum$chr,"_",FM_Liver_hum$variant_pos)
Liver_chr11$index<-paste0("chr",Liver_chr11$hum_chr,"_",Liver_chr11$hum_loci)
chr11_Liver_hum<-NULL
for(i in 1:length(same_Liver_chr11)){
  a<-Liver_chr11[match(same_Liver_chr11[i],Liver_chr11$pig_pos),]
  b<-FM_Liver_hum[match(intersect(a$index,FM_Liver_hum$index),FM_Liver_hum$index),]
  chr11_Liver_hum<-rbind(chr11_Liver_hum,b)
}
chr11_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr11_Liver_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr11_Liver_hum_one2one<-chr11_Liver_hum[match(intersect(chr11_pig2hum_one2one,chr11_Liver_hum$gene_id),chr11_Liver_hum$gene_id),]

chr12_Liver_overloci_one2one<-chr12_Liver_overloci[match(intersect(chr12_Liver_genes,one2one_pig),chr12_Liver_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Liver_chr12$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Liver_chr12$hum_chr<-tmp$V1
Liver_chr12$hum_loci<-tmp$V3
FM_Liver_hum$index<-paste0(FM_Liver_hum$chr,"_",FM_Liver_hum$variant_pos)
Liver_chr12$index<-paste0("chr",Liver_chr12$hum_chr,"_",Liver_chr12$hum_loci)
chr12_Liver_hum<-NULL
for(i in 1:length(same_Liver_chr12)){
  a<-Liver_chr12[match(same_Liver_chr12[i],Liver_chr12$pig_pos),]
  b<-FM_Liver_hum[match(intersect(a$index,FM_Liver_hum$index),FM_Liver_hum$index),]
  chr12_Liver_hum<-rbind(chr12_Liver_hum,b)
}
chr12_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr12_Liver_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr12_Liver_hum_one2one<-chr12_Liver_hum[match(intersect(chr12_pig2hum_one2one,chr12_Liver_hum$gene_id),chr12_Liver_hum$gene_id),]

chr13_Liver_overloci_one2one<-chr13_Liver_overloci[match(intersect(chr13_Liver_genes,one2one_pig),chr13_Liver_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Liver_chr13$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Liver_chr13$hum_chr<-tmp$V1
Liver_chr13$hum_loci<-tmp$V3
FM_Liver_hum$index<-paste0(FM_Liver_hum$chr,"_",FM_Liver_hum$variant_pos)
Liver_chr13$index<-paste0("chr",Liver_chr13$hum_chr,"_",Liver_chr13$hum_loci)
chr13_Liver_hum<-NULL
for(i in 1:length(same_Liver_chr13)){
  a<-Liver_chr13[match(same_Liver_chr13[i],Liver_chr13$pig_pos),]
  b<-FM_Liver_hum[match(intersect(a$index,FM_Liver_hum$index),FM_Liver_hum$index),]
  chr13_Liver_hum<-rbind(chr13_Liver_hum,b)
}
chr13_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr13_Liver_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr13_Liver_hum_one2one<-chr13_Liver_hum[match(intersect(chr13_pig2hum_one2one,chr13_Liver_hum$gene_id),chr13_Liver_hum$gene_id),]

chr14_Liver_overloci_one2one<-chr14_Liver_overloci[match(intersect(chr14_Liver_genes,one2one_pig),chr14_Liver_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Liver_chr14$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Liver_chr14$hum_chr<-tmp$V1
Liver_chr14$hum_loci<-tmp$V3
FM_Liver_hum$index<-paste0(FM_Liver_hum$chr,"_",FM_Liver_hum$variant_pos)
Liver_chr14$index<-paste0("chr",Liver_chr14$hum_chr,"_",Liver_chr14$hum_loci)
chr14_Liver_hum<-NULL
for(i in 1:length(same_Liver_chr14)){
  a<-Liver_chr14[match(same_Liver_chr14[i],Liver_chr14$pig_pos),]
  b<-FM_Liver_hum[match(intersect(a$index,FM_Liver_hum$index),FM_Liver_hum$index),]
  chr14_Liver_hum<-rbind(chr14_Liver_hum,b)
}
chr14_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr14_Liver_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr14_Liver_hum_one2one<-chr14_Liver_hum[match(intersect(chr14_pig2hum_one2one,chr14_Liver_hum$gene_id),chr14_Liver_hum$gene_id),]

chr15_Liver_overloci_one2one<-chr15_Liver_overloci[match(intersect(chr15_Liver_genes,one2one_pig),chr15_Liver_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Liver_chr15$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Liver_chr15$hum_chr<-tmp$V1
Liver_chr15$hum_loci<-tmp$V3
FM_Liver_hum$index<-paste0(FM_Liver_hum$chr,"_",FM_Liver_hum$variant_pos)
Liver_chr15$index<-paste0("chr",Liver_chr15$hum_chr,"_",Liver_chr15$hum_loci)
chr15_Liver_hum<-NULL
for(i in 1:length(same_Liver_chr15)){
  a<-Liver_chr15[match(same_Liver_chr15[i],Liver_chr15$pig_pos),]
  b<-FM_Liver_hum[match(intersect(a$index,FM_Liver_hum$index),FM_Liver_hum$index),]
  chr15_Liver_hum<-rbind(chr15_Liver_hum,b)
}
chr15_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr15_Liver_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr15_Liver_hum_one2one<-chr15_Liver_hum[match(intersect(chr15_pig2hum_one2one,chr15_Liver_hum$gene_id),chr15_Liver_hum$gene_id),]

chr16_Liver_overloci_one2one<-chr16_Liver_overloci[match(intersect(chr16_Liver_genes,one2one_pig),chr16_Liver_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Liver_chr16$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Liver_chr16$hum_chr<-tmp$V1
Liver_chr16$hum_loci<-tmp$V3
FM_Liver_hum$index<-paste0(FM_Liver_hum$chr,"_",FM_Liver_hum$variant_pos)
Liver_chr16$index<-paste0("chr",Liver_chr16$hum_chr,"_",Liver_chr16$hum_loci)
chr16_Liver_hum<-NULL
for(i in 1:length(same_Liver_chr16)){
  a<-Liver_chr16[match(same_Liver_chr16[i],Liver_chr16$pig_pos),]
  b<-FM_Liver_hum[match(intersect(a$index,FM_Liver_hum$index),FM_Liver_hum$index),]
  chr16_Liver_hum<-rbind(chr16_Liver_hum,b)
}
chr16_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr16_Liver_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr16_Liver_hum_one2one<-chr16_Liver_hum[match(intersect(chr16_pig2hum_one2one,chr16_Liver_hum$gene_id),chr16_Liver_hum$gene_id),]

chr17_Liver_overloci_one2one<-chr17_Liver_overloci[match(intersect(chr17_Liver_genes,one2one_pig),chr17_Liver_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Liver_chr17$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Liver_chr17$hum_chr<-tmp$V1
Liver_chr17$hum_loci<-tmp$V3
FM_Liver_hum$index<-paste0(FM_Liver_hum$chr,"_",FM_Liver_hum$variant_pos)
Liver_chr17$index<-paste0("chr",Liver_chr17$hum_chr,"_",Liver_chr17$hum_loci)
chr17_Liver_hum<-NULL
for(i in 1:length(same_Liver_chr17)){
  a<-Liver_chr17[match(same_Liver_chr17[i],Liver_chr17$pig_pos),]
  b<-FM_Liver_hum[match(intersect(a$index,FM_Liver_hum$index),FM_Liver_hum$index),]
  chr17_Liver_hum<-rbind(chr17_Liver_hum,b)
}
chr17_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr17_Liver_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr17_Liver_hum_one2one<-chr17_Liver_hum[match(intersect(chr17_pig2hum_one2one,chr17_Liver_hum$gene_id),chr17_Liver_hum$gene_id),]

chr18_Liver_overloci_one2one<-chr18_Liver_overloci[match(intersect(chr18_Liver_genes,one2one_pig),chr18_Liver_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Liver_chr18$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Liver_chr18$hum_chr<-tmp$V1
Liver_chr18$hum_loci<-tmp$V3
FM_Liver_hum$index<-paste0(FM_Liver_hum$chr,"_",FM_Liver_hum$variant_pos)
Liver_chr18$index<-paste0("chr",Liver_chr18$hum_chr,"_",Liver_chr18$hum_loci)
chr18_Liver_hum<-NULL
for(i in 1:length(same_Liver_chr18)){
  a<-Liver_chr18[match(same_Liver_chr18[i],Liver_chr18$pig_pos),]
  b<-FM_Liver_hum[match(intersect(a$index,FM_Liver_hum$index),FM_Liver_hum$index),]
  chr18_Liver_hum<-rbind(chr18_Liver_hum,b)
}
chr18_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr18_Liver_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr18_Liver_hum_one2one<-chr18_Liver_hum[match(intersect(chr18_pig2hum_one2one,chr18_Liver_hum$gene_id),chr18_Liver_hum$gene_id),]

Liver_one2one_SNP_hum<-rbind(chr1_Liver_hum_one2one,chr2_Liver_hum_one2one,chr3_Liver_hum_one2one,chr4_Liver_hum_one2one,chr5_Liver_hum_one2one,
                             chr6_Liver_hum_one2one,chr7_Liver_hum_one2one,chr8_Liver_hum_one2one,chr9_Liver_hum_one2one,chr10_Liver_hum_one2one,
                             chr11_Liver_hum_one2one,chr12_Liver_hum_one2one,chr13_Liver_hum_one2one,chr14_Liver_hum_one2one,chr15_Liver_hum_one2one,
                             chr16_Liver_hum_one2one,chr17_Liver_hum_one2one,chr18_Liver_hum_one2one)

chr1_Liver_pig_one2one<-chr1_Liver_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr1_Liver_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr1_Liver_overloci_one2one$phenotype_id),1:9]
chr2_Liver_pig_one2one<-chr2_Liver_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr2_Liver_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr2_Liver_overloci_one2one$phenotype_id),1:9]
chr3_Liver_pig_one2one<-chr3_Liver_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr3_Liver_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr3_Liver_overloci_one2one$phenotype_id),1:9]
chr4_Liver_pig_one2one<-chr4_Liver_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr4_Liver_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr4_Liver_overloci_one2one$phenotype_id),1:9]
chr5_Liver_pig_one2one<-chr5_Liver_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr5_Liver_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr5_Liver_overloci_one2one$phenotype_id),1:9]
chr6_Liver_pig_one2one<-chr6_Liver_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr6_Liver_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr6_Liver_overloci_one2one$phenotype_id),1:9]
chr7_Liver_pig_one2one<-chr7_Liver_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr7_Liver_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr7_Liver_overloci_one2one$phenotype_id),1:9]
chr8_Liver_pig_one2one<-chr8_Liver_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr8_Liver_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr8_Liver_overloci_one2one$phenotype_id),1:9]
chr9_Liver_pig_one2one<-chr9_Liver_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr9_Liver_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr9_Liver_overloci_one2one$phenotype_id),1:9]
chr10_Liver_pig_one2one<-chr10_Liver_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr10_Liver_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr10_Liver_overloci_one2one$phenotype_id),1:9]
chr11_Liver_pig_one2one<-chr11_Liver_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr11_Liver_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr11_Liver_overloci_one2one$phenotype_id),1:9]
chr12_Liver_pig_one2one<-chr12_Liver_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr12_Liver_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr12_Liver_overloci_one2one$phenotype_id),1:9]
chr13_Liver_pig_one2one<-chr13_Liver_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr13_Liver_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr13_Liver_overloci_one2one$phenotype_id),1:9]
chr14_Liver_pig_one2one<-chr14_Liver_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr14_Liver_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr14_Liver_overloci_one2one$phenotype_id),1:9]
chr15_Liver_pig_one2one<-chr15_Liver_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr15_Liver_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr15_Liver_overloci_one2one$phenotype_id),1:9]
chr16_Liver_pig_one2one<-chr16_Liver_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr16_Liver_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr16_Liver_overloci_one2one$phenotype_id),1:9]
chr17_Liver_pig_one2one<-chr17_Liver_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr17_Liver_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr17_Liver_overloci_one2one$phenotype_id),1:9]
chr18_Liver_pig_one2one<-chr18_Liver_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr18_Liver_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr18_Liver_overloci_one2one$phenotype_id),1:9]

Liver_one2one_SNP_pig<-rbind(chr1_Liver_pig_one2one,chr2_Liver_pig_one2one,chr3_Liver_pig_one2one,chr4_Liver_pig_one2one,chr5_Liver_pig_one2one,
                             chr6_Liver_pig_one2one,chr7_Liver_pig_one2one,chr8_Liver_pig_one2one,chr9_Liver_pig_one2one,chr10_Liver_pig_one2one,
                             chr11_Liver_pig_one2one,chr12_Liver_pig_one2one,chr13_Liver_pig_one2one,chr14_Liver_pig_one2one,chr15_Liver_pig_one2one,
                             chr16_Liver_pig_one2one,chr17_Liver_pig_one2one,chr18_Liver_pig_one2one)

Liver_SNP_sum<-array(NA,dim=c(nrow(Liver_one2one_SNP_hum),2))
colnames(Liver_SNP_sum)<-c("Human","Pig")
Liver_SNP_sum<-as.data.frame(Liver_SNP_sum)
Liver_SNP_sum$Human<-Liver_one2one_SNP_hum$slope / Liver_one2one_SNP_hum$slope_se
Liver_SNP_sum$Pig<-Liver_one2one_SNP_pig$slope / Liver_one2one_SNP_pig$slope_se
cor<-cor(abs(Liver_SNP_sum$Human),abs(Liver_SNP_sum$Pig))
p_val<-t.test(abs(Liver_SNP_sum$Human),abs(Liver_SNP_sum$Pig))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/Liver_SNP_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(Liver_SNP_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The number of eGenes expressed in common tissues",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

save(Liver_one2one_SNP_hum,Liver_one2one_SNP_pig,file="/Users/baizhonghao/Downloads/human-pig GTEx/Liver_SNP.Rdata")

#Lung_SNP_overlaploci#
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

eqtl_Lung_chr1<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Lung/Lung.cis_qtl_pairs.1.txt"))
eqtl_Lung_chr2<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Lung/Lung.cis_qtl_pairs.2.txt"))
eqtl_Lung_chr3<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Lung/Lung.cis_qtl_pairs.3.txt"))
eqtl_Lung_chr4<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Lung/Lung.cis_qtl_pairs.4.txt"))
eqtl_Lung_chr5<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Lung/Lung.cis_qtl_pairs.5.txt"))
eqtl_Lung_chr6<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Lung/Lung.cis_qtl_pairs.6.txt"))
eqtl_Lung_chr7<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Lung/Lung.cis_qtl_pairs.7.txt"))
eqtl_Lung_chr8<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Lung/Lung.cis_qtl_pairs.8.txt"))
eqtl_Lung_chr9<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Lung/Lung.cis_qtl_pairs.9.txt"))
eqtl_Lung_chr10<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Lung/Lung.cis_qtl_pairs.10.txt"))
eqtl_Lung_chr11<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Lung/Lung.cis_qtl_pairs.11.txt"))
eqtl_Lung_chr12<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Lung/Lung.cis_qtl_pairs.12.txt"))
eqtl_Lung_chr13<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Lung/Lung.cis_qtl_pairs.13.txt"))
eqtl_Lung_chr14<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Lung/Lung.cis_qtl_pairs.14.txt"))
eqtl_Lung_chr15<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Lung/Lung.cis_qtl_pairs.15.txt"))
eqtl_Lung_chr16<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Lung/Lung.cis_qtl_pairs.16.txt"))
eqtl_Lung_chr17<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Lung/Lung.cis_qtl_pairs.17.txt"))
eqtl_Lung_chr18<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Lung/Lung.cis_qtl_pairs.18.txt"))

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung_chr1$variant_id,split="_"))))
eqtl_Lung_chr1$chr<-tmp$V1
eqtl_Lung_chr1$loci<-tmp$V2

eqtl_Lung_chr1$index<-paste0(eqtl_Lung_chr1$chr,"-",eqtl_Lung_chr1$loci)
Lung_chr1$index<-paste0(Lung_chr1$chr,"-",Lung_chr1$pig_pos)
same_Lung_chr1<-intersect(Lung_chr1$pig_pos,eqtl_Lung_chr1$loci)


tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung_chr2$variant_id,split="_"))))
eqtl_Lung_chr2$chr<-tmp$V1
eqtl_Lung_chr2$loci<-tmp$V2

eqtl_Lung_chr2$index<-paste0(eqtl_Lung_chr2$chr,"-",eqtl_Lung_chr2$loci)
Lung_chr2$index<-paste0(Lung_chr2$chr,"-",Lung_chr2$pig_pos)
same_Lung_chr2<-intersect(Lung_chr2$pig_pos,eqtl_Lung_chr2$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung_chr3$variant_id,split="_"))))
eqtl_Lung_chr3$chr<-tmp$V1
eqtl_Lung_chr3$loci<-tmp$V2

eqtl_Lung_chr3$index<-paste0(eqtl_Lung_chr3$chr,"-",eqtl_Lung_chr3$loci)
Lung_chr3$index<-paste0(Lung_chr3$chr,"-",Lung_chr3$pig_pos)
same_Lung_chr3<-intersect(Lung_chr3$pig_pos,eqtl_Lung_chr3$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung_chr4$variant_id,split="_"))))
eqtl_Lung_chr4$chr<-tmp$V1
eqtl_Lung_chr4$loci<-tmp$V2

eqtl_Lung_chr4$index<-paste0(eqtl_Lung_chr4$chr,"-",eqtl_Lung_chr4$loci)
Lung_chr4$index<-paste0(Lung_chr4$chr,"-",Lung_chr4$pig_pos)
same_Lung_chr4<-intersect(Lung_chr4$pig_pos,eqtl_Lung_chr4$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung_chr5$variant_id,split="_"))))
eqtl_Lung_chr5$loci<-tmp$V2
same_Lung_chr5<-intersect(Lung_chr5$pig_pos,eqtl_Lung_chr5$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung_chr6$variant_id,split="_"))))
eqtl_Lung_chr6$loci<-tmp$V2
same_Lung_chr6<-intersect(Lung_chr6$pig_pos,eqtl_Lung_chr6$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung_chr7$variant_id,split="_"))))
eqtl_Lung_chr7$loci<-tmp$V2
same_Lung_chr7<-intersect(Lung_chr7$pig_pos,eqtl_Lung_chr7$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung_chr8$variant_id,split="_"))))
eqtl_Lung_chr8$loci<-tmp$V2
same_Lung_chr8<-intersect(Lung_chr8$pig_pos,eqtl_Lung_chr8$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung_chr9$variant_id,split="_"))))
eqtl_Lung_chr9$loci<-tmp$V2
same_Lung_chr9<-intersect(Lung_chr9$pig_pos,eqtl_Lung_chr9$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung_chr10$variant_id,split="_"))))
eqtl_Lung_chr10$loci<-tmp$V2
same_Lung_chr10<-intersect(Lung_chr10$pig_pos,eqtl_Lung_chr10$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung_chr11$variant_id,split="_"))))
eqtl_Lung_chr11$loci<-tmp$V2
same_Lung_chr11<-intersect(Lung_chr11$pig_pos,eqtl_Lung_chr11$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung_chr12$variant_id,split="_"))))
eqtl_Lung_chr12$loci<-tmp$V2
same_Lung_chr12<-intersect(Lung_chr12$pig_pos,eqtl_Lung_chr12$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung_chr13$variant_id,split="_"))))
eqtl_Lung_chr13$loci<-tmp$V2
same_Lung_chr13<-intersect(Lung_chr13$pig_pos,eqtl_Lung_chr13$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung_chr14$variant_id,split="_"))))
eqtl_Lung_chr14$loci<-tmp$V2
same_Lung_chr14<-intersect(Lung_chr14$pig_pos,eqtl_Lung_chr14$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung_chr15$variant_id,split="_"))))
eqtl_Lung_chr15$loci<-tmp$V2
same_Lung_chr15<-intersect(Lung_chr15$pig_pos,eqtl_Lung_chr15$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung_chr16$variant_id,split="_"))))
eqtl_Lung_chr16$loci<-tmp$V2
same_Lung_chr16<-intersect(Lung_chr16$pig_pos,eqtl_Lung_chr16$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung_chr17$variant_id,split="_"))))
eqtl_Lung_chr17$loci<-tmp$V2
same_Lung_chr17<-intersect(Lung_chr17$pig_pos,eqtl_Lung_chr17$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Lung_chr18$variant_id,split="_"))))
eqtl_Lung_chr18$loci<-tmp$V2
same_Lung_chr18<-intersect(Lung_chr18$pig_pos,eqtl_Lung_chr18$loci)

chr1_Lung_genes<-NULL
chr1_Lung_overloci<-NULL
if(length(same_Lung_chr1!=0)){
  for(i in 1:length(same_Lung_chr1)){
    a<-eqtl_Lung_chr1$phenotype_id[grep(same_Lung_chr1[i],eqtl_Lung_chr1$loci)]
    b<-eqtl_Lung_chr1[grep(same_Lung_chr1[i],eqtl_Lung_chr1$loci),]
    chr1_Lung_genes<-c(chr1_Lung_genes,a)
    chr1_Lung_overloci<-rbind(chr1_Lung_overloci,b)
  }
}

chr2_Lung_genes<-NULL
chr2_Lung_overloci<-NULL
if(length(same_Lung_chr2!=0)){
  for(i in 1:length(same_Lung_chr2)){
    a<-eqtl_Lung_chr2$phenotype_id[grep(same_Lung_chr2[i],eqtl_Lung_chr2$loci)]
    b<-eqtl_Lung_chr2[grep(same_Lung_chr2[i],eqtl_Lung_chr2$loci),]
    chr2_Lung_genes<-c(chr2_Lung_genes,a)
    chr2_Lung_overloci<-rbind(chr2_Lung_overloci,b)
  }
}

chr3_Lung_genes<-NULL
chr3_Lung_overloci<-NULL
if(length(same_Lung_chr3!=0)){
  for(i in 1:length(same_Lung_chr3)){
    a<-eqtl_Lung_chr3$phenotype_id[grep(same_Lung_chr3[i],eqtl_Lung_chr3$loci)]
    b<-eqtl_Lung_chr3[grep(same_Lung_chr3[i],eqtl_Lung_chr3$loci),]
    chr3_Lung_genes<-c(chr3_Lung_genes,a)
    chr3_Lung_overloci<-rbind(chr3_Lung_overloci,b)
  }
}

chr4_Lung_genes<-NULL
chr4_Lung_overloci<-NULL
if(length(same_Lung_chr4!=0)){
  for(i in 1:length(same_Lung_chr4)){
    a<-eqtl_Lung_chr4$phenotype_id[grep(same_Lung_chr4[i],eqtl_Lung_chr4$loci)]
    b<-eqtl_Lung_chr4[grep(same_Lung_chr4[i],eqtl_Lung_chr4$loci),]
    chr4_Lung_genes<-c(chr4_Lung_genes,a)
    chr4_Lung_overloci<-rbind(chr4_Lung_overloci,b)
  }
}

chr5_Lung_genes<-NULL
chr5_Lung_overloci<-NULL
if(length(same_Lung_chr5!=0)){
  for(i in 1:length(same_Lung_chr5)){
    a<-eqtl_Lung_chr5$phenotype_id[grep(same_Lung_chr5[i],eqtl_Lung_chr5$loci)]
    b<-eqtl_Lung_chr5[grep(same_Lung_chr5[i],eqtl_Lung_chr5$loci),]
    chr5_Lung_genes<-c(chr5_Lung_genes,a)
    chr5_Lung_overloci<-rbind(chr5_Lung_overloci,b)
  }
}

chr6_Lung_genes<-NULL
chr6_Lung_overloci<-NULL
if(length(same_Lung_chr6!=0)){
  for(i in 1:length(same_Lung_chr6)){
    a<-eqtl_Lung_chr6$phenotype_id[grep(same_Lung_chr6[i],eqtl_Lung_chr6$loci)]
    b<-eqtl_Lung_chr6[grep(same_Lung_chr6[i],eqtl_Lung_chr6$loci),]
    chr6_Lung_genes<-c(chr6_Lung_genes,a)
    chr6_Lung_overloci<-rbind(chr6_Lung_overloci,b)
  }
}

chr7_Lung_genes<-NULL
chr7_Lung_overloci<-NULL
if(length(same_Lung_chr7!=0)){
  for(i in 1:length(same_Lung_chr7)){
    a<-eqtl_Lung_chr7$phenotype_id[grep(same_Lung_chr7[i],eqtl_Lung_chr7$loci)]
    b<-eqtl_Lung_chr7[grep(same_Lung_chr7[i],eqtl_Lung_chr7$loci),]
    chr7_Lung_genes<-c(chr7_Lung_genes,a)
    chr7_Lung_overloci<-rbind(chr7_Lung_overloci,b)
  }
}

chr8_Lung_genes<-NULL
chr8_Lung_overloci<-NULL
if(length(same_Lung_chr8!=0)){
  for(i in 1:length(same_Lung_chr8)){
    a<-eqtl_Lung_chr8$phenotype_id[grep(same_Lung_chr8[i],eqtl_Lung_chr8$loci)]
    b<-eqtl_Lung_chr8[grep(same_Lung_chr8[i],eqtl_Lung_chr8$loci),]
    chr8_Lung_genes<-c(chr8_Lung_genes,a)
    chr8_Lung_overloci<-rbind(chr8_Lung_overloci,b)
  }
}


chr9_Lung_genes<-NULL
chr9_Lung_overloci<-NULL
if(length(same_Lung_chr9!=0)){
  for(i in 1:length(same_Lung_chr9)){
    a<-eqtl_Lung_chr9$phenotype_id[grep(same_Lung_chr9[i],eqtl_Lung_chr9$loci)]
    b<-eqtl_Lung_chr9[grep(same_Lung_chr9[i],eqtl_Lung_chr9$loci),]
    chr9_Lung_genes<-c(chr9_Lung_genes,a)
    chr9_Lung_overloci<-rbind(chr9_Lung_overloci,b)
  }
}

chr10_Lung_genes<-NULL
chr10_Lung_overloci<-NULL
if(length(same_Lung_chr10!=0)){
  for(i in 1:length(same_Lung_chr10)){
    a<-eqtl_Lung_chr10$phenotype_id[grep(same_Lung_chr10[i],eqtl_Lung_chr10$loci)]
    b<-eqtl_Lung_chr10[grep(same_Lung_chr10[i],eqtl_Lung_chr10$loci),]
    chr10_Lung_genes<-c(chr10_Lung_genes,a)
    chr10_Lung_overloci<-rbind(chr10_Lung_overloci,b)
  }
}

chr11_Lung_genes<-NULL
chr11_Lung_overloci<-NULL
if(length(same_Lung_chr11!=0)){
  for(i in 1:length(same_Lung_chr11)){
    a<-eqtl_Lung_chr11$phenotype_id[grep(same_Lung_chr11[i],eqtl_Lung_chr11$loci)]
    b<-eqtl_Lung_chr11[grep(same_Lung_chr11[i],eqtl_Lung_chr11$loci),]
    chr11_Lung_genes<-c(chr11_Lung_genes,a)
    chr11_Lung_overloci<-rbind(chr11_Lung_overloci,b)
  }
}

chr12_Lung_genes<-NULL
chr12_Lung_overloci<-NULL
if(length(same_Lung_chr12!=0)){
  for(i in 1:length(same_Lung_chr12)){
    a<-eqtl_Lung_chr12$phenotype_id[grep(same_Lung_chr12[i],eqtl_Lung_chr12$loci)]
    b<-eqtl_Lung_chr12[grep(same_Lung_chr12[i],eqtl_Lung_chr12$loci),]
    chr12_Lung_genes<-c(chr12_Lung_genes,a)
    chr12_Lung_overloci<-rbind(chr12_Lung_overloci,b)
  }
}

chr13_Lung_genes<-NULL
chr13_Lung_overloci<-NULL
if(length(same_Lung_chr13!=0)){
  for(i in 1:length(same_Lung_chr13)){
    a<-eqtl_Lung_chr13$phenotype_id[grep(same_Lung_chr13[i],eqtl_Lung_chr13$loci)]
    b<-eqtl_Lung_chr13[grep(same_Lung_chr13[i],eqtl_Lung_chr13$loci),]
    chr13_Lung_genes<-c(chr13_Lung_genes,a)
    chr13_Lung_overloci<-rbind(chr13_Lung_overloci,b)
  }
}

chr14_Lung_genes<-NULL
chr14_Lung_overloci<-NULL
if(length(same_Lung_chr14!=0)){
  for(i in 1:length(same_Lung_chr14)){
    a<-eqtl_Lung_chr14$phenotype_id[grep(same_Lung_chr14[i],eqtl_Lung_chr14$loci)]
    b<-eqtl_Lung_chr14[grep(same_Lung_chr14[i],eqtl_Lung_chr14$loci),]
    chr14_Lung_genes<-c(chr14_Lung_genes,a)
    chr14_Lung_overloci<-rbind(chr14_Lung_overloci,b)
  }
}

chr15_Lung_genes<-NULL
chr15_Lung_overloci<-NULL
if(length(same_Lung_chr15!=0)){
  for(i in 1:length(same_Lung_chr15)){
    a<-eqtl_Lung_chr15$phenotype_id[grep(same_Lung_chr15[i],eqtl_Lung_chr15$loci)]
    b<-eqtl_Lung_chr15[grep(same_Lung_chr15[i],eqtl_Lung_chr15$loci),]
    chr15_Lung_genes<-c(chr15_Lung_genes,a)
    chr15_Lung_overloci<-rbind(chr15_Lung_overloci,b)
  }
}

chr16_Lung_genes<-NULL
chr16_Lung_overloci<-NULL
if(length(same_Lung_chr16!=0)){
  for(i in 1:length(same_Lung_chr16)){
    a<-eqtl_Lung_chr16$phenotype_id[grep(same_Lung_chr16[i],eqtl_Lung_chr16$loci)]
    b<-eqtl_Lung_chr16[grep(same_Lung_chr16[i],eqtl_Lung_chr16$loci),]
    chr16_Lung_genes<-c(chr16_Lung_genes,a)
    chr16_Lung_overloci<-rbind(chr16_Lung_overloci,b)
  }
}

chr17_Lung_genes<-NULL
chr17_Lung_overloci<-NULL
if(length(same_Lung_chr17!=0)){
  for(i in 1:length(same_Lung_chr17)){
    a<-eqtl_Lung_chr17$phenotype_id[grep(same_Lung_chr17[i],eqtl_Lung_chr17$loci)]
    b<-eqtl_Lung_chr17[grep(same_Lung_chr17[i],eqtl_Lung_chr17$loci),]
    chr17_Lung_genes<-c(chr17_Lung_genes,a)
    chr17_Lung_overloci<-rbind(chr17_Lung_overloci,b)
  }
}

chr18_Lung_genes<-NULL
chr18_Lung_overloci<-NULL
if(length(same_Lung_chr18!=0)){
  for(i in 1:length(same_Lung_chr18)){
    a<-eqtl_Lung_chr18$phenotype_id[grep(same_Lung_chr18[i],eqtl_Lung_chr18$loci)]
    b<-eqtl_Lung_chr18[grep(same_Lung_chr18[i],eqtl_Lung_chr18$loci),]
    chr18_Lung_genes<-c(chr18_Lung_genes,a)
    chr18_Lung_overloci<-rbind(chr18_Lung_overloci,b)
  }
}
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
one2one_pig<-annotation$`Pig gene stable ID`[which(annotation$`Pig homology type`=="ortholog_one2one")]

chr1_Lung_overloci_one2one<-chr1_Lung_overloci[match(intersect(chr1_Lung_genes,one2one_pig),chr1_Lung_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Lung_chr1$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Lung_chr1$hum_chr<-tmp$V1
Lung_chr1$hum_loci<-tmp$V3
FM_Lung_hum$index<-paste0(FM_Lung_hum$chr,"_",FM_Lung_hum$variant_pos)
Lung_chr1$index<-paste0("chr",Lung_chr1$hum_chr,"_",Lung_chr1$hum_loci)
chr1_Lung_hum<-NULL
for(i in 1:length(same_Lung_chr1)){
  a<-Lung_chr1[match(same_Lung_chr1[i],Lung_chr1$pig_pos),]
  b<-FM_Lung_hum[match(intersect(a$index,FM_Lung_hum$index),FM_Lung_hum$index),]
  chr1_Lung_hum<-rbind(chr1_Lung_hum,b)
}
chr1_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr1_Lung_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr1_Lung_hum_one2one<-chr1_Lung_hum[match(intersect(chr1_pig2hum_one2one,chr1_Lung_hum$gene_id),chr1_Lung_hum$gene_id),]

chr2_Lung_overloci_one2one<-chr2_Lung_overloci[match(intersect(chr2_Lung_genes,one2one_pig),chr2_Lung_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Lung_chr2$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Lung_chr2$hum_chr<-tmp$V1
Lung_chr2$hum_loci<-tmp$V3
FM_Lung_hum$index<-paste0(FM_Lung_hum$chr,"_",FM_Lung_hum$variant_pos)
Lung_chr2$index<-paste0("chr",Lung_chr2$hum_chr,"_",Lung_chr2$hum_loci)
chr2_Lung_hum<-NULL
for(i in 1:length(same_Lung_chr2)){
  a<-Lung_chr2[match(same_Lung_chr2[i],Lung_chr2$pig_pos),]
  b<-FM_Lung_hum[match(intersect(a$index,FM_Lung_hum$index),FM_Lung_hum$index),]
  chr2_Lung_hum<-rbind(chr2_Lung_hum,b)
}
chr2_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr2_Lung_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr2_Lung_hum_one2one<-chr2_Lung_hum[match(intersect(chr2_pig2hum_one2one,chr2_Lung_hum$gene_id),chr2_Lung_hum$gene_id),]

chr3_Lung_overloci_one2one<-chr3_Lung_overloci[match(intersect(chr3_Lung_genes,one2one_pig),chr3_Lung_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Lung_chr3$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Lung_chr3$hum_chr<-tmp$V1
Lung_chr3$hum_loci<-tmp$V3
FM_Lung_hum$index<-paste0(FM_Lung_hum$chr,"_",FM_Lung_hum$variant_pos)
Lung_chr3$index<-paste0("chr",Lung_chr3$hum_chr,"_",Lung_chr3$hum_loci)
chr3_Lung_hum<-NULL
for(i in 1:length(same_Lung_chr3)){
  a<-Lung_chr3[match(same_Lung_chr3[i],Lung_chr3$pig_pos),]
  b<-FM_Lung_hum[match(intersect(a$index,FM_Lung_hum$index),FM_Lung_hum$index),]
  chr3_Lung_hum<-rbind(chr3_Lung_hum,b)
}
chr3_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr3_Lung_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr3_Lung_hum_one2one<-chr3_Lung_hum[match(intersect(chr3_pig2hum_one2one,chr3_Lung_hum$gene_id),chr3_Lung_hum$gene_id),]

chr4_Lung_overloci_one2one<-chr4_Lung_overloci[match(intersect(chr4_Lung_genes,one2one_pig),chr4_Lung_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Lung_chr4$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Lung_chr4$hum_chr<-tmp$V1
Lung_chr4$hum_loci<-tmp$V3
FM_Lung_hum$index<-paste0(FM_Lung_hum$chr,"_",FM_Lung_hum$variant_pos)
Lung_chr4$index<-paste0("chr",Lung_chr4$hum_chr,"_",Lung_chr4$hum_loci)
chr4_Lung_hum<-NULL
for(i in 1:length(same_Lung_chr4)){
  a<-Lung_chr4[match(same_Lung_chr4[i],Lung_chr4$pig_pos),]
  b<-FM_Lung_hum[match(intersect(a$index,FM_Lung_hum$index),FM_Lung_hum$index),]
  chr4_Lung_hum<-rbind(chr4_Lung_hum,b)
}
chr4_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr4_Lung_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr4_Lung_hum_one2one<-chr4_Lung_hum[match(intersect(chr4_pig2hum_one2one,chr4_Lung_hum$gene_id),chr4_Lung_hum$gene_id),]

chr5_Lung_overloci_one2one<-chr5_Lung_overloci[match(intersect(chr5_Lung_genes,one2one_pig),chr5_Lung_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Lung_chr5$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Lung_chr5$hum_chr<-tmp$V1
Lung_chr5$hum_loci<-tmp$V3
FM_Lung_hum$index<-paste0(FM_Lung_hum$chr,"_",FM_Lung_hum$variant_pos)
Lung_chr5$index<-paste0("chr",Lung_chr5$hum_chr,"_",Lung_chr5$hum_loci)
chr5_Lung_hum<-NULL
for(i in 1:length(same_Lung_chr5)){
  a<-Lung_chr5[match(same_Lung_chr5[i],Lung_chr5$pig_pos),]
  b<-FM_Lung_hum[match(intersect(a$index,FM_Lung_hum$index),FM_Lung_hum$index),]
  chr5_Lung_hum<-rbind(chr5_Lung_hum,b)
}
chr5_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr5_Lung_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr5_Lung_hum_one2one<-chr5_Lung_hum[match(intersect(chr5_pig2hum_one2one,chr5_Lung_hum$gene_id),chr5_Lung_hum$gene_id),]

chr6_Lung_overloci_one2one<-chr6_Lung_overloci[match(intersect(chr6_Lung_genes,one2one_pig),chr6_Lung_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Lung_chr6$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Lung_chr6$hum_chr<-tmp$V1
Lung_chr6$hum_loci<-tmp$V3
FM_Lung_hum$index<-paste0(FM_Lung_hum$chr,"_",FM_Lung_hum$variant_pos)
Lung_chr6$index<-paste0("chr",Lung_chr6$hum_chr,"_",Lung_chr6$hum_loci)
chr6_Lung_hum<-NULL
for(i in 1:length(same_Lung_chr6)){
  a<-Lung_chr6[match(same_Lung_chr6[i],Lung_chr6$pig_pos),]
  b<-FM_Lung_hum[match(intersect(a$index,FM_Lung_hum$index),FM_Lung_hum$index),]
  chr6_Lung_hum<-rbind(chr6_Lung_hum,b)
}
chr6_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr6_Lung_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr6_Lung_hum_one2one<-chr6_Lung_hum[match(intersect(chr6_pig2hum_one2one,chr6_Lung_hum$gene_id),chr6_Lung_hum$gene_id),]

chr7_Lung_overloci_one2one<-chr7_Lung_overloci[match(intersect(chr7_Lung_genes,one2one_pig),chr7_Lung_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Lung_chr7$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Lung_chr7$hum_chr<-tmp$V1
Lung_chr7$hum_loci<-tmp$V3
FM_Lung_hum$index<-paste0(FM_Lung_hum$chr,"_",FM_Lung_hum$variant_pos)
Lung_chr7$index<-paste0("chr",Lung_chr7$hum_chr,"_",Lung_chr7$hum_loci)
chr7_Lung_hum<-NULL
for(i in 1:length(same_Lung_chr7)){
  a<-Lung_chr7[match(same_Lung_chr7[i],Lung_chr7$pig_pos),]
  b<-FM_Lung_hum[match(intersect(a$index,FM_Lung_hum$index),FM_Lung_hum$index),]
  chr7_Lung_hum<-rbind(chr7_Lung_hum,b)
}
chr7_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr7_Lung_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr7_Lung_hum_one2one<-chr7_Lung_hum[match(intersect(chr7_pig2hum_one2one,chr7_Lung_hum$gene_id),chr7_Lung_hum$gene_id),]

chr8_Lung_overloci_one2one<-chr8_Lung_overloci[match(intersect(chr8_Lung_genes,one2one_pig),chr8_Lung_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Lung_chr8$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Lung_chr8$hum_chr<-tmp$V1
Lung_chr8$hum_loci<-tmp$V3
FM_Lung_hum$index<-paste0(FM_Lung_hum$chr,"_",FM_Lung_hum$variant_pos)
Lung_chr8$index<-paste0("chr",Lung_chr8$hum_chr,"_",Lung_chr8$hum_loci)
chr8_Lung_hum<-NULL
for(i in 1:length(same_Lung_chr8)){
  a<-Lung_chr8[match(same_Lung_chr8[i],Lung_chr8$pig_pos),]
  b<-FM_Lung_hum[match(intersect(a$index,FM_Lung_hum$index),FM_Lung_hum$index),]
  chr8_Lung_hum<-rbind(chr8_Lung_hum,b)
}
chr8_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr8_Lung_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr8_Lung_hum_one2one<-chr8_Lung_hum[match(intersect(chr8_pig2hum_one2one,chr8_Lung_hum$gene_id),chr8_Lung_hum$gene_id),]

chr9_Lung_overloci_one2one<-chr9_Lung_overloci[match(intersect(chr9_Lung_genes,one2one_pig),chr9_Lung_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Lung_chr9$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Lung_chr9$hum_chr<-tmp$V1
Lung_chr9$hum_loci<-tmp$V3
FM_Lung_hum$index<-paste0(FM_Lung_hum$chr,"_",FM_Lung_hum$variant_pos)
Lung_chr9$index<-paste0("chr",Lung_chr9$hum_chr,"_",Lung_chr9$hum_loci)
chr9_Lung_hum<-NULL
for(i in 1:length(same_Lung_chr9)){
  a<-Lung_chr9[match(same_Lung_chr9[i],Lung_chr9$pig_pos),]
  b<-FM_Lung_hum[match(intersect(a$index,FM_Lung_hum$index),FM_Lung_hum$index),]
  chr9_Lung_hum<-rbind(chr9_Lung_hum,b)
}
chr9_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr9_Lung_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr9_Lung_hum_one2one<-chr9_Lung_hum[match(intersect(chr9_pig2hum_one2one,chr9_Lung_hum$gene_id),chr9_Lung_hum$gene_id),]

chr10_Lung_overloci_one2one<-chr10_Lung_overloci[match(intersect(chr10_Lung_genes,one2one_pig),chr10_Lung_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Lung_chr10$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Lung_chr10$hum_chr<-tmp$V1
Lung_chr10$hum_loci<-tmp$V3
FM_Lung_hum$index<-paste0(FM_Lung_hum$chr,"_",FM_Lung_hum$variant_pos)
Lung_chr10$index<-paste0("chr",Lung_chr10$hum_chr,"_",Lung_chr10$hum_loci)
chr10_Lung_hum<-NULL
for(i in 1:length(same_Lung_chr10)){
  a<-Lung_chr10[match(same_Lung_chr10[i],Lung_chr10$pig_pos),]
  b<-FM_Lung_hum[match(intersect(a$index,FM_Lung_hum$index),FM_Lung_hum$index),]
  chr10_Lung_hum<-rbind(chr10_Lung_hum,b)
}
chr10_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr10_Lung_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr10_Lung_hum_one2one<-chr10_Lung_hum[match(intersect(chr10_pig2hum_one2one,chr10_Lung_hum$gene_id),chr10_Lung_hum$gene_id),]

chr11_Lung_overloci_one2one<-chr11_Lung_overloci[match(intersect(chr11_Lung_genes,one2one_pig),chr11_Lung_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Lung_chr11$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Lung_chr11$hum_chr<-tmp$V1
Lung_chr11$hum_loci<-tmp$V3
FM_Lung_hum$index<-paste0(FM_Lung_hum$chr,"_",FM_Lung_hum$variant_pos)
Lung_chr11$index<-paste0("chr",Lung_chr11$hum_chr,"_",Lung_chr11$hum_loci)
chr11_Lung_hum<-NULL
for(i in 1:length(same_Lung_chr11)){
  a<-Lung_chr11[match(same_Lung_chr11[i],Lung_chr11$pig_pos),]
  b<-FM_Lung_hum[match(intersect(a$index,FM_Lung_hum$index),FM_Lung_hum$index),]
  chr11_Lung_hum<-rbind(chr11_Lung_hum,b)
}
chr11_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr11_Lung_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr11_Lung_hum_one2one<-chr11_Lung_hum[match(intersect(chr11_pig2hum_one2one,chr11_Lung_hum$gene_id),chr11_Lung_hum$gene_id),]

chr12_Lung_overloci_one2one<-chr12_Lung_overloci[match(intersect(chr12_Lung_genes,one2one_pig),chr12_Lung_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Lung_chr12$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Lung_chr12$hum_chr<-tmp$V1
Lung_chr12$hum_loci<-tmp$V3
FM_Lung_hum$index<-paste0(FM_Lung_hum$chr,"_",FM_Lung_hum$variant_pos)
Lung_chr12$index<-paste0("chr",Lung_chr12$hum_chr,"_",Lung_chr12$hum_loci)
chr12_Lung_hum<-NULL
for(i in 1:length(same_Lung_chr12)){
  a<-Lung_chr12[match(same_Lung_chr12[i],Lung_chr12$pig_pos),]
  b<-FM_Lung_hum[match(intersect(a$index,FM_Lung_hum$index),FM_Lung_hum$index),]
  chr12_Lung_hum<-rbind(chr12_Lung_hum,b)
}
chr12_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr12_Lung_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr12_Lung_hum_one2one<-chr12_Lung_hum[match(intersect(chr12_pig2hum_one2one,chr12_Lung_hum$gene_id),chr12_Lung_hum$gene_id),]

chr13_Lung_overloci_one2one<-chr13_Lung_overloci[match(intersect(chr13_Lung_genes,one2one_pig),chr13_Lung_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Lung_chr13$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Lung_chr13$hum_chr<-tmp$V1
Lung_chr13$hum_loci<-tmp$V3
FM_Lung_hum$index<-paste0(FM_Lung_hum$chr,"_",FM_Lung_hum$variant_pos)
Lung_chr13$index<-paste0("chr",Lung_chr13$hum_chr,"_",Lung_chr13$hum_loci)
chr13_Lung_hum<-NULL
for(i in 1:length(same_Lung_chr13)){
  a<-Lung_chr13[match(same_Lung_chr13[i],Lung_chr13$pig_pos),]
  b<-FM_Lung_hum[match(intersect(a$index,FM_Lung_hum$index),FM_Lung_hum$index),]
  chr13_Lung_hum<-rbind(chr13_Lung_hum,b)
}
chr13_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr13_Lung_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr13_Lung_hum_one2one<-chr13_Lung_hum[match(intersect(chr13_pig2hum_one2one,chr13_Lung_hum$gene_id),chr13_Lung_hum$gene_id),]

chr14_Lung_overloci_one2one<-chr14_Lung_overloci[match(intersect(chr14_Lung_genes,one2one_pig),chr14_Lung_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Lung_chr14$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Lung_chr14$hum_chr<-tmp$V1
Lung_chr14$hum_loci<-tmp$V3
FM_Lung_hum$index<-paste0(FM_Lung_hum$chr,"_",FM_Lung_hum$variant_pos)
Lung_chr14$index<-paste0("chr",Lung_chr14$hum_chr,"_",Lung_chr14$hum_loci)
chr14_Lung_hum<-NULL
for(i in 1:length(same_Lung_chr14)){
  a<-Lung_chr14[match(same_Lung_chr14[i],Lung_chr14$pig_pos),]
  b<-FM_Lung_hum[match(intersect(a$index,FM_Lung_hum$index),FM_Lung_hum$index),]
  chr14_Lung_hum<-rbind(chr14_Lung_hum,b)
}
chr14_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr14_Lung_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr14_Lung_hum_one2one<-chr14_Lung_hum[match(intersect(chr14_pig2hum_one2one,chr14_Lung_hum$gene_id),chr14_Lung_hum$gene_id),]

chr15_Lung_overloci_one2one<-chr15_Lung_overloci[match(intersect(chr15_Lung_genes,one2one_pig),chr15_Lung_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Lung_chr15$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Lung_chr15$hum_chr<-tmp$V1
Lung_chr15$hum_loci<-tmp$V3
FM_Lung_hum$index<-paste0(FM_Lung_hum$chr,"_",FM_Lung_hum$variant_pos)
Lung_chr15$index<-paste0("chr",Lung_chr15$hum_chr,"_",Lung_chr15$hum_loci)
chr15_Lung_hum<-NULL
for(i in 1:length(same_Lung_chr15)){
  a<-Lung_chr15[match(same_Lung_chr15[i],Lung_chr15$pig_pos),]
  b<-FM_Lung_hum[match(intersect(a$index,FM_Lung_hum$index),FM_Lung_hum$index),]
  chr15_Lung_hum<-rbind(chr15_Lung_hum,b)
}
chr15_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr15_Lung_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr15_Lung_hum_one2one<-chr15_Lung_hum[match(intersect(chr15_pig2hum_one2one,chr15_Lung_hum$gene_id),chr15_Lung_hum$gene_id),]

chr16_Lung_overloci_one2one<-chr16_Lung_overloci[match(intersect(chr16_Lung_genes,one2one_pig),chr16_Lung_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Lung_chr16$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Lung_chr16$hum_chr<-tmp$V1
Lung_chr16$hum_loci<-tmp$V3
FM_Lung_hum$index<-paste0(FM_Lung_hum$chr,"_",FM_Lung_hum$variant_pos)
Lung_chr16$index<-paste0("chr",Lung_chr16$hum_chr,"_",Lung_chr16$hum_loci)
chr16_Lung_hum<-NULL
for(i in 1:length(same_Lung_chr16)){
  a<-Lung_chr16[match(same_Lung_chr16[i],Lung_chr16$pig_pos),]
  b<-FM_Lung_hum[match(intersect(a$index,FM_Lung_hum$index),FM_Lung_hum$index),]
  chr16_Lung_hum<-rbind(chr16_Lung_hum,b)
}
chr16_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr16_Lung_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr16_Lung_hum_one2one<-chr16_Lung_hum[match(intersect(chr16_pig2hum_one2one,chr16_Lung_hum$gene_id),chr16_Lung_hum$gene_id),]

chr17_Lung_overloci_one2one<-chr17_Lung_overloci[match(intersect(chr17_Lung_genes,one2one_pig),chr17_Lung_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Lung_chr17$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Lung_chr17$hum_chr<-tmp$V1
Lung_chr17$hum_loci<-tmp$V3
FM_Lung_hum$index<-paste0(FM_Lung_hum$chr,"_",FM_Lung_hum$variant_pos)
Lung_chr17$index<-paste0("chr",Lung_chr17$hum_chr,"_",Lung_chr17$hum_loci)
chr17_Lung_hum<-NULL
for(i in 1:length(same_Lung_chr17)){
  a<-Lung_chr17[match(same_Lung_chr17[i],Lung_chr17$pig_pos),]
  b<-FM_Lung_hum[match(intersect(a$index,FM_Lung_hum$index),FM_Lung_hum$index),]
  chr17_Lung_hum<-rbind(chr17_Lung_hum,b)
}
chr17_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr17_Lung_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr17_Lung_hum_one2one<-chr17_Lung_hum[match(intersect(chr17_pig2hum_one2one,chr17_Lung_hum$gene_id),chr17_Lung_hum$gene_id),]

chr18_Lung_overloci_one2one<-chr18_Lung_overloci[match(intersect(chr18_Lung_genes,one2one_pig),chr18_Lung_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Lung_chr18$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Lung_chr18$hum_chr<-tmp$V1
Lung_chr18$hum_loci<-tmp$V3
FM_Lung_hum$index<-paste0(FM_Lung_hum$chr,"_",FM_Lung_hum$variant_pos)
Lung_chr18$index<-paste0("chr",Lung_chr18$hum_chr,"_",Lung_chr18$hum_loci)
chr18_Lung_hum<-NULL
for(i in 1:length(same_Lung_chr18)){
  a<-Lung_chr18[match(same_Lung_chr18[i],Lung_chr18$pig_pos),]
  b<-FM_Lung_hum[match(intersect(a$index,FM_Lung_hum$index),FM_Lung_hum$index),]
  chr18_Lung_hum<-rbind(chr18_Lung_hum,b)
}
chr18_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr18_Lung_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr18_Lung_hum_one2one<-chr18_Lung_hum[match(intersect(chr18_pig2hum_one2one,chr18_Lung_hum$gene_id),chr18_Lung_hum$gene_id),]

Lung_one2one_SNP_hum<-rbind(chr1_Lung_hum_one2one,chr2_Lung_hum_one2one,chr3_Lung_hum_one2one,chr4_Lung_hum_one2one,chr5_Lung_hum_one2one,
                            chr6_Lung_hum_one2one,chr7_Lung_hum_one2one,chr8_Lung_hum_one2one,chr9_Lung_hum_one2one,chr10_Lung_hum_one2one,
                            chr11_Lung_hum_one2one,chr12_Lung_hum_one2one,chr13_Lung_hum_one2one,chr14_Lung_hum_one2one,chr15_Lung_hum_one2one,
                            chr16_Lung_hum_one2one,chr17_Lung_hum_one2one,chr18_Lung_hum_one2one)

chr1_Lung_pig_one2one<-chr1_Lung_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr1_Lung_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr1_Lung_overloci_one2one$phenotype_id),1:9]
chr2_Lung_pig_one2one<-chr2_Lung_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr2_Lung_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr2_Lung_overloci_one2one$phenotype_id),1:9]
chr3_Lung_pig_one2one<-chr3_Lung_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr3_Lung_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr3_Lung_overloci_one2one$phenotype_id),1:9]
chr4_Lung_pig_one2one<-chr4_Lung_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr4_Lung_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr4_Lung_overloci_one2one$phenotype_id),1:9]
chr5_Lung_pig_one2one<-chr5_Lung_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr5_Lung_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr5_Lung_overloci_one2one$phenotype_id),1:9]
chr6_Lung_pig_one2one<-chr6_Lung_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr6_Lung_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr6_Lung_overloci_one2one$phenotype_id),1:9]
chr7_Lung_pig_one2one<-chr7_Lung_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr7_Lung_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr7_Lung_overloci_one2one$phenotype_id),1:9]
chr8_Lung_pig_one2one<-chr8_Lung_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr8_Lung_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr8_Lung_overloci_one2one$phenotype_id),1:9]
chr9_Lung_pig_one2one<-chr9_Lung_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr9_Lung_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr9_Lung_overloci_one2one$phenotype_id),1:9]
chr10_Lung_pig_one2one<-chr10_Lung_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr10_Lung_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr10_Lung_overloci_one2one$phenotype_id),1:9]
chr11_Lung_pig_one2one<-chr11_Lung_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr11_Lung_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr11_Lung_overloci_one2one$phenotype_id),1:9]
chr12_Lung_pig_one2one<-chr12_Lung_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr12_Lung_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr12_Lung_overloci_one2one$phenotype_id),1:9]
chr13_Lung_pig_one2one<-chr13_Lung_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr13_Lung_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr13_Lung_overloci_one2one$phenotype_id),1:9]
chr14_Lung_pig_one2one<-chr14_Lung_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr14_Lung_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr14_Lung_overloci_one2one$phenotype_id),1:9]
chr15_Lung_pig_one2one<-chr15_Lung_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr15_Lung_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr15_Lung_overloci_one2one$phenotype_id),1:9]
chr16_Lung_pig_one2one<-chr16_Lung_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr16_Lung_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr16_Lung_overloci_one2one$phenotype_id),1:9]
chr17_Lung_pig_one2one<-chr17_Lung_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr17_Lung_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr17_Lung_overloci_one2one$phenotype_id),1:9]
chr18_Lung_pig_one2one<-chr18_Lung_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr18_Lung_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr18_Lung_overloci_one2one$phenotype_id),1:9]

Lung_one2one_SNP_pig<-rbind(chr1_Lung_pig_one2one,chr2_Lung_pig_one2one,chr3_Lung_pig_one2one,chr4_Lung_pig_one2one,chr5_Lung_pig_one2one,
                            chr6_Lung_pig_one2one,chr7_Lung_pig_one2one,chr8_Lung_pig_one2one,chr9_Lung_pig_one2one,chr10_Lung_pig_one2one,
                            chr11_Lung_pig_one2one,chr12_Lung_pig_one2one,chr13_Lung_pig_one2one,chr14_Lung_pig_one2one,chr15_Lung_pig_one2one,
                            chr16_Lung_pig_one2one,chr17_Lung_pig_one2one,chr18_Lung_pig_one2one)

Lung_SNP_sum<-array(NA,dim=c(nrow(Lung_one2one_SNP_hum),2))
colnames(Lung_SNP_sum)<-c("Human","Pig")
Lung_SNP_sum<-as.data.frame(Lung_SNP_sum)
Lung_SNP_sum$Human<-Lung_one2one_SNP_hum$slope / Lung_one2one_SNP_hum$slope_se
Lung_SNP_sum$Pig<-Lung_one2one_SNP_pig$slope / Lung_one2one_SNP_pig$slope_se
cor<-cor(abs(Lung_SNP_sum$Human),abs(Lung_SNP_sum$Pig))
p_val<-t.test(abs(Lung_SNP_sum$Human),abs(Lung_SNP_sum$Pig))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/Lung_SNP_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(Lung_SNP_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The number of eGenes expressed in common tissues",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

save(Lung_one2one_SNP_hum,Lung_one2one_SNP_pig,file="/Users/baizhonghao/Downloads/human-pig GTEx/Lung_SNP.Rdata")

#Ovary_SNP_overlaploci#
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

eqtl_Ovary_chr1<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ovary/Ovary.cis_qtl_pairs.1.txt"))
eqtl_Ovary_chr2<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ovary/Ovary.cis_qtl_pairs.2.txt"))
eqtl_Ovary_chr3<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ovary/Ovary.cis_qtl_pairs.3.txt"))
eqtl_Ovary_chr4<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ovary/Ovary.cis_qtl_pairs.4.txt"))
eqtl_Ovary_chr5<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ovary/Ovary.cis_qtl_pairs.5.txt"))
eqtl_Ovary_chr6<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ovary/Ovary.cis_qtl_pairs.6.txt"))
eqtl_Ovary_chr7<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ovary/Ovary.cis_qtl_pairs.7.txt"))
eqtl_Ovary_chr8<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ovary/Ovary.cis_qtl_pairs.8.txt"))
eqtl_Ovary_chr9<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ovary/Ovary.cis_qtl_pairs.9.txt"))
eqtl_Ovary_chr10<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ovary/Ovary.cis_qtl_pairs.10.txt"))
eqtl_Ovary_chr11<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ovary/Ovary.cis_qtl_pairs.11.txt"))
eqtl_Ovary_chr12<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ovary/Ovary.cis_qtl_pairs.12.txt"))
eqtl_Ovary_chr13<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ovary/Ovary.cis_qtl_pairs.13.txt"))
eqtl_Ovary_chr14<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ovary/Ovary.cis_qtl_pairs.14.txt"))
eqtl_Ovary_chr15<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ovary/Ovary.cis_qtl_pairs.15.txt"))
eqtl_Ovary_chr16<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ovary/Ovary.cis_qtl_pairs.16.txt"))
eqtl_Ovary_chr17<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ovary/Ovary.cis_qtl_pairs.17.txt"))
eqtl_Ovary_chr18<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Ovary/Ovary.cis_qtl_pairs.18.txt"))

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary_chr1$variant_id,split="_"))))
eqtl_Ovary_chr1$chr<-tmp$V1
eqtl_Ovary_chr1$loci<-tmp$V2

eqtl_Ovary_chr1$index<-paste0(eqtl_Ovary_chr1$chr,"-",eqtl_Ovary_chr1$loci)
Ovary_chr1$index<-paste0(Ovary_chr1$chr,"-",Ovary_chr1$pig_pos)
same_Ovary_chr1<-intersect(Ovary_chr1$pig_pos,eqtl_Ovary_chr1$loci)


tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary_chr2$variant_id,split="_"))))
eqtl_Ovary_chr2$chr<-tmp$V1
eqtl_Ovary_chr2$loci<-tmp$V2

eqtl_Ovary_chr2$index<-paste0(eqtl_Ovary_chr2$chr,"-",eqtl_Ovary_chr2$loci)
Ovary_chr2$index<-paste0(Ovary_chr2$chr,"-",Ovary_chr2$pig_pos)
same_Ovary_chr2<-intersect(Ovary_chr2$pig_pos,eqtl_Ovary_chr2$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary_chr3$variant_id,split="_"))))
eqtl_Ovary_chr3$chr<-tmp$V1
eqtl_Ovary_chr3$loci<-tmp$V2

eqtl_Ovary_chr3$index<-paste0(eqtl_Ovary_chr3$chr,"-",eqtl_Ovary_chr3$loci)
Ovary_chr3$index<-paste0(Ovary_chr3$chr,"-",Ovary_chr3$pig_pos)
same_Ovary_chr3<-intersect(Ovary_chr3$pig_pos,eqtl_Ovary_chr3$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary_chr4$variant_id,split="_"))))
eqtl_Ovary_chr4$chr<-tmp$V1
eqtl_Ovary_chr4$loci<-tmp$V2

eqtl_Ovary_chr4$index<-paste0(eqtl_Ovary_chr4$chr,"-",eqtl_Ovary_chr4$loci)
Ovary_chr4$index<-paste0(Ovary_chr4$chr,"-",Ovary_chr4$pig_pos)
same_Ovary_chr4<-intersect(Ovary_chr4$pig_pos,eqtl_Ovary_chr4$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary_chr5$variant_id,split="_"))))
eqtl_Ovary_chr5$loci<-tmp$V2
same_Ovary_chr5<-intersect(Ovary_chr5$pig_pos,eqtl_Ovary_chr5$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary_chr6$variant_id,split="_"))))
eqtl_Ovary_chr6$loci<-tmp$V2
same_Ovary_chr6<-intersect(Ovary_chr6$pig_pos,eqtl_Ovary_chr6$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary_chr7$variant_id,split="_"))))
eqtl_Ovary_chr7$loci<-tmp$V2
same_Ovary_chr7<-intersect(Ovary_chr7$pig_pos,eqtl_Ovary_chr7$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary_chr8$variant_id,split="_"))))
eqtl_Ovary_chr8$loci<-tmp$V2
same_Ovary_chr8<-intersect(Ovary_chr8$pig_pos,eqtl_Ovary_chr8$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary_chr9$variant_id,split="_"))))
eqtl_Ovary_chr9$loci<-tmp$V2
same_Ovary_chr9<-intersect(Ovary_chr9$pig_pos,eqtl_Ovary_chr9$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary_chr10$variant_id,split="_"))))
eqtl_Ovary_chr10$loci<-tmp$V2
same_Ovary_chr10<-intersect(Ovary_chr10$pig_pos,eqtl_Ovary_chr10$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary_chr11$variant_id,split="_"))))
eqtl_Ovary_chr11$loci<-tmp$V2
same_Ovary_chr11<-intersect(Ovary_chr11$pig_pos,eqtl_Ovary_chr11$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary_chr12$variant_id,split="_"))))
eqtl_Ovary_chr12$loci<-tmp$V2
same_Ovary_chr12<-intersect(Ovary_chr12$pig_pos,eqtl_Ovary_chr12$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary_chr13$variant_id,split="_"))))
eqtl_Ovary_chr13$loci<-tmp$V2
same_Ovary_chr13<-intersect(Ovary_chr13$pig_pos,eqtl_Ovary_chr13$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary_chr14$variant_id,split="_"))))
eqtl_Ovary_chr14$loci<-tmp$V2
same_Ovary_chr14<-intersect(Ovary_chr14$pig_pos,eqtl_Ovary_chr14$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary_chr15$variant_id,split="_"))))
eqtl_Ovary_chr15$loci<-tmp$V2
same_Ovary_chr15<-intersect(Ovary_chr15$pig_pos,eqtl_Ovary_chr15$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary_chr16$variant_id,split="_"))))
eqtl_Ovary_chr16$loci<-tmp$V2
same_Ovary_chr16<-intersect(Ovary_chr16$pig_pos,eqtl_Ovary_chr16$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary_chr17$variant_id,split="_"))))
eqtl_Ovary_chr17$loci<-tmp$V2
same_Ovary_chr17<-intersect(Ovary_chr17$pig_pos,eqtl_Ovary_chr17$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Ovary_chr18$variant_id,split="_"))))
eqtl_Ovary_chr18$loci<-tmp$V2
same_Ovary_chr18<-intersect(Ovary_chr18$pig_pos,eqtl_Ovary_chr18$loci)

chr1_Ovary_genes<-NULL
chr1_Ovary_overloci<-NULL
if(length(same_Ovary_chr1!=0)){
  for(i in 1:length(same_Ovary_chr1)){
    a<-eqtl_Ovary_chr1$phenotype_id[grep(same_Ovary_chr1[i],eqtl_Ovary_chr1$loci)]
    b<-eqtl_Ovary_chr1[grep(same_Ovary_chr1[i],eqtl_Ovary_chr1$loci),]
    chr1_Ovary_genes<-c(chr1_Ovary_genes,a)
    chr1_Ovary_overloci<-rbind(chr1_Ovary_overloci,b)
  }
}

chr2_Ovary_genes<-NULL
chr2_Ovary_overloci<-NULL
if(length(same_Ovary_chr2!=0)){
  for(i in 1:length(same_Ovary_chr2)){
    a<-eqtl_Ovary_chr2$phenotype_id[grep(same_Ovary_chr2[i],eqtl_Ovary_chr2$loci)]
    b<-eqtl_Ovary_chr2[grep(same_Ovary_chr2[i],eqtl_Ovary_chr2$loci),]
    chr2_Ovary_genes<-c(chr2_Ovary_genes,a)
    chr2_Ovary_overloci<-rbind(chr2_Ovary_overloci,b)
  }
}

chr3_Ovary_genes<-NULL
chr3_Ovary_overloci<-NULL
if(length(same_Ovary_chr3!=0)){
  for(i in 1:length(same_Ovary_chr3)){
    a<-eqtl_Ovary_chr3$phenotype_id[grep(same_Ovary_chr3[i],eqtl_Ovary_chr3$loci)]
    b<-eqtl_Ovary_chr3[grep(same_Ovary_chr3[i],eqtl_Ovary_chr3$loci),]
    chr3_Ovary_genes<-c(chr3_Ovary_genes,a)
    chr3_Ovary_overloci<-rbind(chr3_Ovary_overloci,b)
  }
}

chr4_Ovary_genes<-NULL
chr4_Ovary_overloci<-NULL
if(length(same_Ovary_chr4!=0)){
  for(i in 1:length(same_Ovary_chr4)){
    a<-eqtl_Ovary_chr4$phenotype_id[grep(same_Ovary_chr4[i],eqtl_Ovary_chr4$loci)]
    b<-eqtl_Ovary_chr4[grep(same_Ovary_chr4[i],eqtl_Ovary_chr4$loci),]
    chr4_Ovary_genes<-c(chr4_Ovary_genes,a)
    chr4_Ovary_overloci<-rbind(chr4_Ovary_overloci,b)
  }
}

chr5_Ovary_genes<-NULL
chr5_Ovary_overloci<-NULL
if(length(same_Ovary_chr5!=0)){
  for(i in 1:length(same_Ovary_chr5)){
    a<-eqtl_Ovary_chr5$phenotype_id[grep(same_Ovary_chr5[i],eqtl_Ovary_chr5$loci)]
    b<-eqtl_Ovary_chr5[grep(same_Ovary_chr5[i],eqtl_Ovary_chr5$loci),]
    chr5_Ovary_genes<-c(chr5_Ovary_genes,a)
    chr5_Ovary_overloci<-rbind(chr5_Ovary_overloci,b)
  }
}

chr6_Ovary_genes<-NULL
chr6_Ovary_overloci<-NULL
if(length(same_Ovary_chr6!=0)){
  for(i in 1:length(same_Ovary_chr6)){
    a<-eqtl_Ovary_chr6$phenotype_id[grep(same_Ovary_chr6[i],eqtl_Ovary_chr6$loci)]
    b<-eqtl_Ovary_chr6[grep(same_Ovary_chr6[i],eqtl_Ovary_chr6$loci),]
    chr6_Ovary_genes<-c(chr6_Ovary_genes,a)
    chr6_Ovary_overloci<-rbind(chr6_Ovary_overloci,b)
  }
}

chr7_Ovary_genes<-NULL
chr7_Ovary_overloci<-NULL
if(length(same_Ovary_chr7!=0)){
  for(i in 1:length(same_Ovary_chr7)){
    a<-eqtl_Ovary_chr7$phenotype_id[grep(same_Ovary_chr7[i],eqtl_Ovary_chr7$loci)]
    b<-eqtl_Ovary_chr7[grep(same_Ovary_chr7[i],eqtl_Ovary_chr7$loci),]
    chr7_Ovary_genes<-c(chr7_Ovary_genes,a)
    chr7_Ovary_overloci<-rbind(chr7_Ovary_overloci,b)
  }
}

chr8_Ovary_genes<-NULL
chr8_Ovary_overloci<-NULL
if(length(same_Ovary_chr8!=0)){
  for(i in 1:length(same_Ovary_chr8)){
    a<-eqtl_Ovary_chr8$phenotype_id[grep(same_Ovary_chr8[i],eqtl_Ovary_chr8$loci)]
    b<-eqtl_Ovary_chr8[grep(same_Ovary_chr8[i],eqtl_Ovary_chr8$loci),]
    chr8_Ovary_genes<-c(chr8_Ovary_genes,a)
    chr8_Ovary_overloci<-rbind(chr8_Ovary_overloci,b)
  }
}

chr9_Ovary_genes<-NULL
chr9_Ovary_overloci<-NULL
if(length(same_Ovary_chr9!=0)){
  for(i in 1:length(same_Ovary_chr9)){
    a<-eqtl_Ovary_chr9$phenotype_id[grep(same_Ovary_chr9[i],eqtl_Ovary_chr9$loci)]
    b<-eqtl_Ovary_chr9[grep(same_Ovary_chr9[i],eqtl_Ovary_chr9$loci),]
    chr9_Ovary_genes<-c(chr9_Ovary_genes,a)
    chr9_Ovary_overloci<-rbind(chr9_Ovary_overloci,b)
  }
}

chr10_Ovary_genes<-NULL
chr10_Ovary_overloci<-NULL
if(length(same_Ovary_chr10!=0)){
  for(i in 1:length(same_Ovary_chr10)){
    a<-eqtl_Ovary_chr10$phenotype_id[grep(same_Ovary_chr10[i],eqtl_Ovary_chr10$loci)]
    b<-eqtl_Ovary_chr10[grep(same_Ovary_chr10[i],eqtl_Ovary_chr10$loci),]
    chr10_Ovary_genes<-c(chr10_Ovary_genes,a)
    chr10_Ovary_overloci<-rbind(chr10_Ovary_overloci,b)
  }
}

chr11_Ovary_genes<-NULL
chr11_Ovary_overloci<-NULL
if(length(same_Ovary_chr11!=0)){
  for(i in 1:length(same_Ovary_chr11)){
    a<-eqtl_Ovary_chr11$phenotype_id[grep(same_Ovary_chr11[i],eqtl_Ovary_chr11$loci)]
    b<-eqtl_Ovary_chr11[grep(same_Ovary_chr11[i],eqtl_Ovary_chr11$loci),]
    chr11_Ovary_genes<-c(chr11_Ovary_genes,a)
    chr11_Ovary_overloci<-rbind(chr11_Ovary_overloci,b)
  }
}

chr12_Ovary_genes<-NULL
chr12_Ovary_overloci<-NULL
if(length(same_Ovary_chr12!=0)){
  for(i in 1:length(same_Ovary_chr12)){
    a<-eqtl_Ovary_chr12$phenotype_id[grep(same_Ovary_chr12[i],eqtl_Ovary_chr12$loci)]
    b<-eqtl_Ovary_chr12[grep(same_Ovary_chr12[i],eqtl_Ovary_chr12$loci),]
    chr12_Ovary_genes<-c(chr12_Ovary_genes,a)
    chr12_Ovary_overloci<-rbind(chr12_Ovary_overloci,b)
  }
}

chr13_Ovary_genes<-NULL
chr13_Ovary_overloci<-NULL
if(length(same_Ovary_chr13!=0)){
  for(i in 1:length(same_Ovary_chr13)){
    a<-eqtl_Ovary_chr13$phenotype_id[grep(same_Ovary_chr13[i],eqtl_Ovary_chr13$loci)]
    b<-eqtl_Ovary_chr13[grep(same_Ovary_chr13[i],eqtl_Ovary_chr13$loci),]
    chr13_Ovary_genes<-c(chr13_Ovary_genes,a)
    chr13_Ovary_overloci<-rbind(chr13_Ovary_overloci,b)
  }
}

chr14_Ovary_genes<-NULL
chr14_Ovary_overloci<-NULL
if(length(same_Ovary_chr14!=0)){
  for(i in 1:length(same_Ovary_chr14)){
    a<-eqtl_Ovary_chr14$phenotype_id[grep(same_Ovary_chr14[i],eqtl_Ovary_chr14$loci)]
    b<-eqtl_Ovary_chr14[grep(same_Ovary_chr14[i],eqtl_Ovary_chr14$loci),]
    chr14_Ovary_genes<-c(chr14_Ovary_genes,a)
    chr14_Ovary_overloci<-rbind(chr14_Ovary_overloci,b)
  }
}

chr15_Ovary_genes<-NULL
chr15_Ovary_overloci<-NULL
if(length(same_Ovary_chr15!=0)){
  for(i in 1:length(same_Ovary_chr15)){
    a<-eqtl_Ovary_chr15$phenotype_id[grep(same_Ovary_chr15[i],eqtl_Ovary_chr15$loci)]
    b<-eqtl_Ovary_chr15[grep(same_Ovary_chr15[i],eqtl_Ovary_chr15$loci),]
    chr15_Ovary_genes<-c(chr15_Ovary_genes,a)
    chr15_Ovary_overloci<-rbind(chr15_Ovary_overloci,b)
  }
}

chr16_Ovary_genes<-NULL
chr16_Ovary_overloci<-NULL
if(length(same_Ovary_chr16!=0)){
  for(i in 1:length(same_Ovary_chr16)){
    a<-eqtl_Ovary_chr16$phenotype_id[grep(same_Ovary_chr16[i],eqtl_Ovary_chr16$loci)]
    b<-eqtl_Ovary_chr16[grep(same_Ovary_chr16[i],eqtl_Ovary_chr16$loci),]
    chr16_Ovary_genes<-c(chr16_Ovary_genes,a)
    chr16_Ovary_overloci<-rbind(chr16_Ovary_overloci,b)
  }
}

chr17_Ovary_genes<-NULL
chr17_Ovary_overloci<-NULL
if(length(same_Ovary_chr17!=0)){
  for(i in 1:length(same_Ovary_chr17)){
    a<-eqtl_Ovary_chr17$phenotype_id[grep(same_Ovary_chr17[i],eqtl_Ovary_chr17$loci)]
    b<-eqtl_Ovary_chr17[grep(same_Ovary_chr17[i],eqtl_Ovary_chr17$loci),]
    chr17_Ovary_genes<-c(chr17_Ovary_genes,a)
    chr17_Ovary_overloci<-rbind(chr17_Ovary_overloci,b)
  }
}

chr18_Ovary_genes<-NULL
chr18_Ovary_overloci<-NULL
if(length(same_Ovary_chr18!=0)){
  for(i in 1:length(same_Ovary_chr18)){
    a<-eqtl_Ovary_chr18$phenotype_id[grep(same_Ovary_chr18[i],eqtl_Ovary_chr18$loci)]
    b<-eqtl_Ovary_chr18[grep(same_Ovary_chr18[i],eqtl_Ovary_chr18$loci),]
    chr18_Ovary_genes<-c(chr18_Ovary_genes,a)
    chr18_Ovary_overloci<-rbind(chr18_Ovary_overloci,b)
  }
}
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
one2one_pig<-annotation$`Pig gene stable ID`[which(annotation$`Pig homology type`=="ortholog_one2one")]

chr1_Ovary_overloci_one2one<-chr1_Ovary_overloci[match(intersect(chr1_Ovary_genes,one2one_pig),chr1_Ovary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ovary_chr1$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ovary_chr1$hum_chr<-tmp$V1
Ovary_chr1$hum_loci<-tmp$V3
FM_Ovary_hum$index<-paste0(FM_Ovary_hum$chr,"_",FM_Ovary_hum$variant_pos)
Ovary_chr1$index<-paste0("chr",Ovary_chr1$hum_chr,"_",Ovary_chr1$hum_loci)
chr1_Ovary_hum<-NULL
for(i in 1:length(same_Ovary_chr1)){
  a<-Ovary_chr1[match(same_Ovary_chr1[i],Ovary_chr1$pig_pos),]
  b<-FM_Ovary_hum[match(intersect(a$index,FM_Ovary_hum$index),FM_Ovary_hum$index),]
  chr1_Ovary_hum<-rbind(chr1_Ovary_hum,b)
}
chr1_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr1_Ovary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr1_Ovary_hum_one2one<-chr1_Ovary_hum[match(intersect(chr1_pig2hum_one2one,chr1_Ovary_hum$gene_id),chr1_Ovary_hum$gene_id),]

chr2_Ovary_overloci_one2one<-chr2_Ovary_overloci[match(intersect(chr2_Ovary_genes,one2one_pig),chr2_Ovary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ovary_chr2$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ovary_chr2$hum_chr<-tmp$V1
Ovary_chr2$hum_loci<-tmp$V3
FM_Ovary_hum$index<-paste0(FM_Ovary_hum$chr,"_",FM_Ovary_hum$variant_pos)
Ovary_chr2$index<-paste0("chr",Ovary_chr2$hum_chr,"_",Ovary_chr2$hum_loci)
chr2_Ovary_hum<-NULL
for(i in 1:length(same_Ovary_chr2)){
  a<-Ovary_chr2[match(same_Ovary_chr2[i],Ovary_chr2$pig_pos),]
  b<-FM_Ovary_hum[match(intersect(a$index,FM_Ovary_hum$index),FM_Ovary_hum$index),]
  chr2_Ovary_hum<-rbind(chr2_Ovary_hum,b)
}
chr2_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr2_Ovary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr2_Ovary_hum_one2one<-chr2_Ovary_hum[match(intersect(chr2_pig2hum_one2one,chr2_Ovary_hum$gene_id),chr2_Ovary_hum$gene_id),]

chr3_Ovary_overloci_one2one<-chr3_Ovary_overloci[match(intersect(chr3_Ovary_genes,one2one_pig),chr3_Ovary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ovary_chr3$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ovary_chr3$hum_chr<-tmp$V1
Ovary_chr3$hum_loci<-tmp$V3
FM_Ovary_hum$index<-paste0(FM_Ovary_hum$chr,"_",FM_Ovary_hum$variant_pos)
Ovary_chr3$index<-paste0("chr",Ovary_chr3$hum_chr,"_",Ovary_chr3$hum_loci)
chr3_Ovary_hum<-NULL
for(i in 1:length(same_Ovary_chr3)){
  a<-Ovary_chr3[match(same_Ovary_chr3[i],Ovary_chr3$pig_pos),]
  b<-FM_Ovary_hum[match(intersect(a$index,FM_Ovary_hum$index),FM_Ovary_hum$index),]
  chr3_Ovary_hum<-rbind(chr3_Ovary_hum,b)
}
chr3_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr3_Ovary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr3_Ovary_hum_one2one<-chr3_Ovary_hum[match(intersect(chr3_pig2hum_one2one,chr3_Ovary_hum$gene_id),chr3_Ovary_hum$gene_id),]

chr4_Ovary_overloci_one2one<-chr4_Ovary_overloci[match(intersect(chr4_Ovary_genes,one2one_pig),chr4_Ovary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ovary_chr4$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ovary_chr4$hum_chr<-tmp$V1
Ovary_chr4$hum_loci<-tmp$V3
FM_Ovary_hum$index<-paste0(FM_Ovary_hum$chr,"_",FM_Ovary_hum$variant_pos)
Ovary_chr4$index<-paste0("chr",Ovary_chr4$hum_chr,"_",Ovary_chr4$hum_loci)
chr4_Ovary_hum<-NULL
for(i in 1:length(same_Ovary_chr4)){
  a<-Ovary_chr4[match(same_Ovary_chr4[i],Ovary_chr4$pig_pos),]
  b<-FM_Ovary_hum[match(intersect(a$index,FM_Ovary_hum$index),FM_Ovary_hum$index),]
  chr4_Ovary_hum<-rbind(chr4_Ovary_hum,b)
}
chr4_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr4_Ovary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr4_Ovary_hum_one2one<-chr4_Ovary_hum[match(intersect(chr4_pig2hum_one2one,chr4_Ovary_hum$gene_id),chr4_Ovary_hum$gene_id),]

chr5_Ovary_overloci_one2one<-chr5_Ovary_overloci[match(intersect(chr5_Ovary_genes,one2one_pig),chr5_Ovary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ovary_chr5$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ovary_chr5$hum_chr<-tmp$V1
Ovary_chr5$hum_loci<-tmp$V3
FM_Ovary_hum$index<-paste0(FM_Ovary_hum$chr,"_",FM_Ovary_hum$variant_pos)
Ovary_chr5$index<-paste0("chr",Ovary_chr5$hum_chr,"_",Ovary_chr5$hum_loci)
chr5_Ovary_hum<-NULL
for(i in 1:length(same_Ovary_chr5)){
  a<-Ovary_chr5[match(same_Ovary_chr5[i],Ovary_chr5$pig_pos),]
  b<-FM_Ovary_hum[match(intersect(a$index,FM_Ovary_hum$index),FM_Ovary_hum$index),]
  chr5_Ovary_hum<-rbind(chr5_Ovary_hum,b)
}
chr5_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr5_Ovary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr5_Ovary_hum_one2one<-chr5_Ovary_hum[match(intersect(chr5_pig2hum_one2one,chr5_Ovary_hum$gene_id),chr5_Ovary_hum$gene_id),]

chr6_Ovary_overloci_one2one<-chr6_Ovary_overloci[match(intersect(chr6_Ovary_genes,one2one_pig),chr6_Ovary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ovary_chr6$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ovary_chr6$hum_chr<-tmp$V1
Ovary_chr6$hum_loci<-tmp$V3
FM_Ovary_hum$index<-paste0(FM_Ovary_hum$chr,"_",FM_Ovary_hum$variant_pos)
Ovary_chr6$index<-paste0("chr",Ovary_chr6$hum_chr,"_",Ovary_chr6$hum_loci)
chr6_Ovary_hum<-NULL
for(i in 1:length(same_Ovary_chr6)){
  a<-Ovary_chr6[match(same_Ovary_chr6[i],Ovary_chr6$pig_pos),]
  b<-FM_Ovary_hum[match(intersect(a$index,FM_Ovary_hum$index),FM_Ovary_hum$index),]
  chr6_Ovary_hum<-rbind(chr6_Ovary_hum,b)
}
chr6_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr6_Ovary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr6_Ovary_hum_one2one<-chr6_Ovary_hum[match(intersect(chr6_pig2hum_one2one,chr6_Ovary_hum$gene_id),chr6_Ovary_hum$gene_id),]

chr7_Ovary_overloci_one2one<-chr7_Ovary_overloci[match(intersect(chr7_Ovary_genes,one2one_pig),chr7_Ovary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ovary_chr7$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ovary_chr7$hum_chr<-tmp$V1
Ovary_chr7$hum_loci<-tmp$V3
FM_Ovary_hum$index<-paste0(FM_Ovary_hum$chr,"_",FM_Ovary_hum$variant_pos)
Ovary_chr7$index<-paste0("chr",Ovary_chr7$hum_chr,"_",Ovary_chr7$hum_loci)
chr7_Ovary_hum<-NULL
for(i in 1:length(same_Ovary_chr7)){
  a<-Ovary_chr7[match(same_Ovary_chr7[i],Ovary_chr7$pig_pos),]
  b<-FM_Ovary_hum[match(intersect(a$index,FM_Ovary_hum$index),FM_Ovary_hum$index),]
  chr7_Ovary_hum<-rbind(chr7_Ovary_hum,b)
}
chr7_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr7_Ovary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr7_Ovary_hum_one2one<-chr7_Ovary_hum[match(intersect(chr7_pig2hum_one2one,chr7_Ovary_hum$gene_id),chr7_Ovary_hum$gene_id),]

chr8_Ovary_overloci_one2one<-chr8_Ovary_overloci[match(intersect(chr8_Ovary_genes,one2one_pig),chr8_Ovary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ovary_chr8$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ovary_chr8$hum_chr<-tmp$V1
Ovary_chr8$hum_loci<-tmp$V3
FM_Ovary_hum$index<-paste0(FM_Ovary_hum$chr,"_",FM_Ovary_hum$variant_pos)
Ovary_chr8$index<-paste0("chr",Ovary_chr8$hum_chr,"_",Ovary_chr8$hum_loci)
chr8_Ovary_hum<-NULL
for(i in 1:length(same_Ovary_chr8)){
  a<-Ovary_chr8[match(same_Ovary_chr8[i],Ovary_chr8$pig_pos),]
  b<-FM_Ovary_hum[match(intersect(a$index,FM_Ovary_hum$index),FM_Ovary_hum$index),]
  chr8_Ovary_hum<-rbind(chr8_Ovary_hum,b)
}
chr8_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr8_Ovary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr8_Ovary_hum_one2one<-chr8_Ovary_hum[match(intersect(chr8_pig2hum_one2one,chr8_Ovary_hum$gene_id),chr8_Ovary_hum$gene_id),]

chr9_Ovary_overloci_one2one<-chr9_Ovary_overloci[match(intersect(chr9_Ovary_genes,one2one_pig),chr9_Ovary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ovary_chr9$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ovary_chr9$hum_chr<-tmp$V1
Ovary_chr9$hum_loci<-tmp$V3
FM_Ovary_hum$index<-paste0(FM_Ovary_hum$chr,"_",FM_Ovary_hum$variant_pos)
Ovary_chr9$index<-paste0("chr",Ovary_chr9$hum_chr,"_",Ovary_chr9$hum_loci)
chr9_Ovary_hum<-NULL
for(i in 1:length(same_Ovary_chr9)){
  a<-Ovary_chr9[match(same_Ovary_chr9[i],Ovary_chr9$pig_pos),]
  b<-FM_Ovary_hum[match(intersect(a$index,FM_Ovary_hum$index),FM_Ovary_hum$index),]
  chr9_Ovary_hum<-rbind(chr9_Ovary_hum,b)
}
chr9_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr9_Ovary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr9_Ovary_hum_one2one<-chr9_Ovary_hum[match(intersect(chr9_pig2hum_one2one,chr9_Ovary_hum$gene_id),chr9_Ovary_hum$gene_id),]

chr10_Ovary_overloci_one2one<-chr10_Ovary_overloci[match(intersect(chr10_Ovary_genes,one2one_pig),chr10_Ovary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ovary_chr10$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ovary_chr10$hum_chr<-tmp$V1
Ovary_chr10$hum_loci<-tmp$V3
FM_Ovary_hum$index<-paste0(FM_Ovary_hum$chr,"_",FM_Ovary_hum$variant_pos)
Ovary_chr10$index<-paste0("chr",Ovary_chr10$hum_chr,"_",Ovary_chr10$hum_loci)
chr10_Ovary_hum<-NULL
for(i in 1:length(same_Ovary_chr10)){
  a<-Ovary_chr10[match(same_Ovary_chr10[i],Ovary_chr10$pig_pos),]
  b<-FM_Ovary_hum[match(intersect(a$index,FM_Ovary_hum$index),FM_Ovary_hum$index),]
  chr10_Ovary_hum<-rbind(chr10_Ovary_hum,b)
}
chr10_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr10_Ovary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr10_Ovary_hum_one2one<-chr10_Ovary_hum[match(intersect(chr10_pig2hum_one2one,chr10_Ovary_hum$gene_id),chr10_Ovary_hum$gene_id),]

chr11_Ovary_overloci_one2one<-chr11_Ovary_overloci[match(intersect(chr11_Ovary_genes,one2one_pig),chr11_Ovary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ovary_chr11$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ovary_chr11$hum_chr<-tmp$V1
Ovary_chr11$hum_loci<-tmp$V3
FM_Ovary_hum$index<-paste0(FM_Ovary_hum$chr,"_",FM_Ovary_hum$variant_pos)
Ovary_chr11$index<-paste0("chr",Ovary_chr11$hum_chr,"_",Ovary_chr11$hum_loci)
chr11_Ovary_hum<-NULL
for(i in 1:length(same_Ovary_chr11)){
  a<-Ovary_chr11[match(same_Ovary_chr11[i],Ovary_chr11$pig_pos),]
  b<-FM_Ovary_hum[match(intersect(a$index,FM_Ovary_hum$index),FM_Ovary_hum$index),]
  chr11_Ovary_hum<-rbind(chr11_Ovary_hum,b)
}
chr11_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr11_Ovary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr11_Ovary_hum_one2one<-chr11_Ovary_hum[match(intersect(chr11_pig2hum_one2one,chr11_Ovary_hum$gene_id),chr11_Ovary_hum$gene_id),]

chr12_Ovary_overloci_one2one<-chr12_Ovary_overloci[match(intersect(chr12_Ovary_genes,one2one_pig),chr12_Ovary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ovary_chr12$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ovary_chr12$hum_chr<-tmp$V1
Ovary_chr12$hum_loci<-tmp$V3
FM_Ovary_hum$index<-paste0(FM_Ovary_hum$chr,"_",FM_Ovary_hum$variant_pos)
Ovary_chr12$index<-paste0("chr",Ovary_chr12$hum_chr,"_",Ovary_chr12$hum_loci)
chr12_Ovary_hum<-NULL
for(i in 1:length(same_Ovary_chr12)){
  a<-Ovary_chr12[match(same_Ovary_chr12[i],Ovary_chr12$pig_pos),]
  b<-FM_Ovary_hum[match(intersect(a$index,FM_Ovary_hum$index),FM_Ovary_hum$index),]
  chr12_Ovary_hum<-rbind(chr12_Ovary_hum,b)
}
chr12_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr12_Ovary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr12_Ovary_hum_one2one<-chr12_Ovary_hum[match(intersect(chr12_pig2hum_one2one,chr12_Ovary_hum$gene_id),chr12_Ovary_hum$gene_id),]

chr13_Ovary_overloci_one2one<-chr13_Ovary_overloci[match(intersect(chr13_Ovary_genes,one2one_pig),chr13_Ovary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ovary_chr13$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ovary_chr13$hum_chr<-tmp$V1
Ovary_chr13$hum_loci<-tmp$V3
FM_Ovary_hum$index<-paste0(FM_Ovary_hum$chr,"_",FM_Ovary_hum$variant_pos)
Ovary_chr13$index<-paste0("chr",Ovary_chr13$hum_chr,"_",Ovary_chr13$hum_loci)
chr13_Ovary_hum<-NULL
for(i in 1:length(same_Ovary_chr13)){
  a<-Ovary_chr13[match(same_Ovary_chr13[i],Ovary_chr13$pig_pos),]
  b<-FM_Ovary_hum[match(intersect(a$index,FM_Ovary_hum$index),FM_Ovary_hum$index),]
  chr13_Ovary_hum<-rbind(chr13_Ovary_hum,b)
}
chr13_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr13_Ovary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr13_Ovary_hum_one2one<-chr13_Ovary_hum[match(intersect(chr13_pig2hum_one2one,chr13_Ovary_hum$gene_id),chr13_Ovary_hum$gene_id),]

chr14_Ovary_overloci_one2one<-chr14_Ovary_overloci[match(intersect(chr14_Ovary_genes,one2one_pig),chr14_Ovary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ovary_chr14$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ovary_chr14$hum_chr<-tmp$V1
Ovary_chr14$hum_loci<-tmp$V3
FM_Ovary_hum$index<-paste0(FM_Ovary_hum$chr,"_",FM_Ovary_hum$variant_pos)
Ovary_chr14$index<-paste0("chr",Ovary_chr14$hum_chr,"_",Ovary_chr14$hum_loci)
chr14_Ovary_hum<-NULL
for(i in 1:length(same_Ovary_chr14)){
  a<-Ovary_chr14[match(same_Ovary_chr14[i],Ovary_chr14$pig_pos),]
  b<-FM_Ovary_hum[match(intersect(a$index,FM_Ovary_hum$index),FM_Ovary_hum$index),]
  chr14_Ovary_hum<-rbind(chr14_Ovary_hum,b)
}
chr14_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr14_Ovary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr14_Ovary_hum_one2one<-chr14_Ovary_hum[match(intersect(chr14_pig2hum_one2one,chr14_Ovary_hum$gene_id),chr14_Ovary_hum$gene_id),]

chr15_Ovary_overloci_one2one<-chr15_Ovary_overloci[match(intersect(chr15_Ovary_genes,one2one_pig),chr15_Ovary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ovary_chr15$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ovary_chr15$hum_chr<-tmp$V1
Ovary_chr15$hum_loci<-tmp$V3
FM_Ovary_hum$index<-paste0(FM_Ovary_hum$chr,"_",FM_Ovary_hum$variant_pos)
Ovary_chr15$index<-paste0("chr",Ovary_chr15$hum_chr,"_",Ovary_chr15$hum_loci)
chr15_Ovary_hum<-NULL
for(i in 1:length(same_Ovary_chr15)){
  a<-Ovary_chr15[match(same_Ovary_chr15[i],Ovary_chr15$pig_pos),]
  b<-FM_Ovary_hum[match(intersect(a$index,FM_Ovary_hum$index),FM_Ovary_hum$index),]
  chr15_Ovary_hum<-rbind(chr15_Ovary_hum,b)
}
chr15_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr15_Ovary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr15_Ovary_hum_one2one<-chr15_Ovary_hum[match(intersect(chr15_pig2hum_one2one,chr15_Ovary_hum$gene_id),chr15_Ovary_hum$gene_id),]

chr16_Ovary_overloci_one2one<-chr16_Ovary_overloci[match(intersect(chr16_Ovary_genes,one2one_pig),chr16_Ovary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ovary_chr16$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ovary_chr16$hum_chr<-tmp$V1
Ovary_chr16$hum_loci<-tmp$V3
FM_Ovary_hum$index<-paste0(FM_Ovary_hum$chr,"_",FM_Ovary_hum$variant_pos)
Ovary_chr16$index<-paste0("chr",Ovary_chr16$hum_chr,"_",Ovary_chr16$hum_loci)
chr16_Ovary_hum<-NULL
for(i in 1:length(same_Ovary_chr16)){
  a<-Ovary_chr16[match(same_Ovary_chr16[i],Ovary_chr16$pig_pos),]
  b<-FM_Ovary_hum[match(intersect(a$index,FM_Ovary_hum$index),FM_Ovary_hum$index),]
  chr16_Ovary_hum<-rbind(chr16_Ovary_hum,b)
}
chr16_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr16_Ovary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr16_Ovary_hum_one2one<-chr16_Ovary_hum[match(intersect(chr16_pig2hum_one2one,chr16_Ovary_hum$gene_id),chr16_Ovary_hum$gene_id),]

chr17_Ovary_overloci_one2one<-chr17_Ovary_overloci[match(intersect(chr17_Ovary_genes,one2one_pig),chr17_Ovary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ovary_chr17$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ovary_chr17$hum_chr<-tmp$V1
Ovary_chr17$hum_loci<-tmp$V3
FM_Ovary_hum$index<-paste0(FM_Ovary_hum$chr,"_",FM_Ovary_hum$variant_pos)
Ovary_chr17$index<-paste0("chr",Ovary_chr17$hum_chr,"_",Ovary_chr17$hum_loci)
chr17_Ovary_hum<-NULL
for(i in 1:length(same_Ovary_chr17)){
  a<-Ovary_chr17[match(same_Ovary_chr17[i],Ovary_chr17$pig_pos),]
  b<-FM_Ovary_hum[match(intersect(a$index,FM_Ovary_hum$index),FM_Ovary_hum$index),]
  chr17_Ovary_hum<-rbind(chr17_Ovary_hum,b)
}
chr17_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr17_Ovary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr17_Ovary_hum_one2one<-chr17_Ovary_hum[match(intersect(chr17_pig2hum_one2one,chr17_Ovary_hum$gene_id),chr17_Ovary_hum$gene_id),]

chr18_Ovary_overloci_one2one<-chr18_Ovary_overloci[match(intersect(chr18_Ovary_genes,one2one_pig),chr18_Ovary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Ovary_chr18$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Ovary_chr18$hum_chr<-tmp$V1
Ovary_chr18$hum_loci<-tmp$V3
FM_Ovary_hum$index<-paste0(FM_Ovary_hum$chr,"_",FM_Ovary_hum$variant_pos)
Ovary_chr18$index<-paste0("chr",Ovary_chr18$hum_chr,"_",Ovary_chr18$hum_loci)
chr18_Ovary_hum<-NULL
for(i in 1:length(same_Ovary_chr18)){
  a<-Ovary_chr18[match(same_Ovary_chr18[i],Ovary_chr18$pig_pos),]
  b<-FM_Ovary_hum[match(intersect(a$index,FM_Ovary_hum$index),FM_Ovary_hum$index),]
  chr18_Ovary_hum<-rbind(chr18_Ovary_hum,b)
}
chr18_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr18_Ovary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr18_Ovary_hum_one2one<-chr18_Ovary_hum[match(intersect(chr18_pig2hum_one2one,chr18_Ovary_hum$gene_id),chr18_Ovary_hum$gene_id),]

Ovary_one2one_SNP_hum<-rbind(chr1_Ovary_hum_one2one,chr2_Ovary_hum_one2one,chr3_Ovary_hum_one2one,chr4_Ovary_hum_one2one,chr5_Ovary_hum_one2one,
                             chr6_Ovary_hum_one2one,chr7_Ovary_hum_one2one,chr8_Ovary_hum_one2one,chr9_Ovary_hum_one2one,chr10_Ovary_hum_one2one,
                             chr11_Ovary_hum_one2one,chr12_Ovary_hum_one2one,chr13_Ovary_hum_one2one,chr14_Ovary_hum_one2one,chr15_Ovary_hum_one2one,
                             chr16_Ovary_hum_one2one,chr17_Ovary_hum_one2one,chr18_Ovary_hum_one2one)

chr1_Ovary_pig_one2one<-chr1_Ovary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr1_Ovary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr1_Ovary_overloci_one2one$phenotype_id),1:9]
chr2_Ovary_pig_one2one<-chr2_Ovary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr2_Ovary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr2_Ovary_overloci_one2one$phenotype_id),1:9]
chr3_Ovary_pig_one2one<-chr3_Ovary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr3_Ovary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr3_Ovary_overloci_one2one$phenotype_id),1:9]
chr4_Ovary_pig_one2one<-chr4_Ovary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr4_Ovary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr4_Ovary_overloci_one2one$phenotype_id),1:9]
chr5_Ovary_pig_one2one<-chr5_Ovary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr5_Ovary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr5_Ovary_overloci_one2one$phenotype_id),1:9]
chr6_Ovary_pig_one2one<-chr6_Ovary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr6_Ovary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr6_Ovary_overloci_one2one$phenotype_id),1:9]
chr7_Ovary_pig_one2one<-chr7_Ovary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr7_Ovary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr7_Ovary_overloci_one2one$phenotype_id),1:9]
chr8_Ovary_pig_one2one<-chr8_Ovary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr8_Ovary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr8_Ovary_overloci_one2one$phenotype_id),1:9]
chr9_Ovary_pig_one2one<-chr9_Ovary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr9_Ovary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr9_Ovary_overloci_one2one$phenotype_id),1:9]
chr10_Ovary_pig_one2one<-chr10_Ovary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr10_Ovary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr10_Ovary_overloci_one2one$phenotype_id),1:9]
chr11_Ovary_pig_one2one<-chr11_Ovary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr11_Ovary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr11_Ovary_overloci_one2one$phenotype_id),1:9]
chr12_Ovary_pig_one2one<-chr12_Ovary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr12_Ovary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr12_Ovary_overloci_one2one$phenotype_id),1:9]
chr13_Ovary_pig_one2one<-chr13_Ovary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr13_Ovary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr13_Ovary_overloci_one2one$phenotype_id),1:9]
chr14_Ovary_pig_one2one<-chr14_Ovary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr14_Ovary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr14_Ovary_overloci_one2one$phenotype_id),1:9]
chr15_Ovary_pig_one2one<-chr15_Ovary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr15_Ovary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr15_Ovary_overloci_one2one$phenotype_id),1:9]
chr16_Ovary_pig_one2one<-chr16_Ovary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr16_Ovary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr16_Ovary_overloci_one2one$phenotype_id),1:9]
chr17_Ovary_pig_one2one<-chr17_Ovary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr17_Ovary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr17_Ovary_overloci_one2one$phenotype_id),1:9]
chr18_Ovary_pig_one2one<-chr18_Ovary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr18_Ovary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr18_Ovary_overloci_one2one$phenotype_id),1:9]

Ovary_one2one_SNP_pig<-rbind(chr1_Ovary_pig_one2one,chr2_Ovary_pig_one2one,chr3_Ovary_pig_one2one,chr4_Ovary_pig_one2one,chr5_Ovary_pig_one2one,
                             chr6_Ovary_pig_one2one,chr7_Ovary_pig_one2one,chr8_Ovary_pig_one2one,chr9_Ovary_pig_one2one,chr10_Ovary_pig_one2one,
                             chr11_Ovary_pig_one2one,chr12_Ovary_pig_one2one,chr13_Ovary_pig_one2one,chr14_Ovary_pig_one2one,chr15_Ovary_pig_one2one,
                             chr16_Ovary_pig_one2one,chr17_Ovary_pig_one2one,chr18_Ovary_pig_one2one)

Ovary_SNP_sum<-array(NA,dim=c(nrow(Ovary_one2one_SNP_hum),2))
colnames(Ovary_SNP_sum)<-c("Human","Pig")
Ovary_SNP_sum<-as.data.frame(Ovary_SNP_sum)
Ovary_SNP_sum$Human<-Ovary_one2one_SNP_hum$slope / Ovary_one2one_SNP_hum$slope_se
Ovary_SNP_sum$Pig<-Ovary_one2one_SNP_pig$slope / Ovary_one2one_SNP_pig$slope_se
cor<-cor(abs(Ovary_SNP_sum$Human),abs(Ovary_SNP_sum$Pig))
p_val<-t.test(abs(Ovary_SNP_sum$Human),abs(Ovary_SNP_sum$Pig))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/Ovary_SNP_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(Ovary_SNP_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The number of eGenes expressed in common tissues",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

save(Ovary_one2one_SNP_hum,Ovary_one2one_SNP_pig,file="/Users/baizhonghao/Downloads/human-pig GTEx/Ovary_SNP.Rdata")

#Pituitary_SNP_overlaploci#
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

eqtl_Pituitary_chr1<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Pituitary/Pituitary.cis_qtl_pairs.1.txt"))
eqtl_Pituitary_chr2<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Pituitary/Pituitary.cis_qtl_pairs.2.txt"))
eqtl_Pituitary_chr3<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Pituitary/Pituitary.cis_qtl_pairs.3.txt"))
eqtl_Pituitary_chr4<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Pituitary/Pituitary.cis_qtl_pairs.4.txt"))
eqtl_Pituitary_chr5<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Pituitary/Pituitary.cis_qtl_pairs.5.txt"))
eqtl_Pituitary_chr6<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Pituitary/Pituitary.cis_qtl_pairs.6.txt"))
eqtl_Pituitary_chr7<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Pituitary/Pituitary.cis_qtl_pairs.7.txt"))
eqtl_Pituitary_chr8<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Pituitary/Pituitary.cis_qtl_pairs.8.txt"))
eqtl_Pituitary_chr9<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Pituitary/Pituitary.cis_qtl_pairs.9.txt"))
eqtl_Pituitary_chr10<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Pituitary/Pituitary.cis_qtl_pairs.10.txt"))
eqtl_Pituitary_chr11<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Pituitary/Pituitary.cis_qtl_pairs.11.txt"))
eqtl_Pituitary_chr12<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Pituitary/Pituitary.cis_qtl_pairs.12.txt"))
eqtl_Pituitary_chr13<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Pituitary/Pituitary.cis_qtl_pairs.13.txt"))
eqtl_Pituitary_chr14<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Pituitary/Pituitary.cis_qtl_pairs.14.txt"))
eqtl_Pituitary_chr15<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Pituitary/Pituitary.cis_qtl_pairs.15.txt"))
eqtl_Pituitary_chr16<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Pituitary/Pituitary.cis_qtl_pairs.16.txt"))
eqtl_Pituitary_chr17<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Pituitary/Pituitary.cis_qtl_pairs.17.txt"))
eqtl_Pituitary_chr18<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Pituitary/Pituitary.cis_qtl_pairs.18.txt"))

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary_chr1$variant_id,split="_"))))
eqtl_Pituitary_chr1$chr<-tmp$V1
eqtl_Pituitary_chr1$loci<-tmp$V2

eqtl_Pituitary_chr1$index<-paste0(eqtl_Pituitary_chr1$chr,"-",eqtl_Pituitary_chr1$loci)
Pituitary_chr1$index<-paste0(Pituitary_chr1$chr,"-",Pituitary_chr1$pig_pos)
same_Pituitary_chr1<-intersect(Pituitary_chr1$pig_pos,eqtl_Pituitary_chr1$loci)


tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary_chr2$variant_id,split="_"))))
eqtl_Pituitary_chr2$chr<-tmp$V1
eqtl_Pituitary_chr2$loci<-tmp$V2

eqtl_Pituitary_chr2$index<-paste0(eqtl_Pituitary_chr2$chr,"-",eqtl_Pituitary_chr2$loci)
Pituitary_chr2$index<-paste0(Pituitary_chr2$chr,"-",Pituitary_chr2$pig_pos)
same_Pituitary_chr2<-intersect(Pituitary_chr2$pig_pos,eqtl_Pituitary_chr2$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary_chr3$variant_id,split="_"))))
eqtl_Pituitary_chr3$chr<-tmp$V1
eqtl_Pituitary_chr3$loci<-tmp$V2

eqtl_Pituitary_chr3$index<-paste0(eqtl_Pituitary_chr3$chr,"-",eqtl_Pituitary_chr3$loci)
Pituitary_chr3$index<-paste0(Pituitary_chr3$chr,"-",Pituitary_chr3$pig_pos)
same_Pituitary_chr3<-intersect(Pituitary_chr3$pig_pos,eqtl_Pituitary_chr3$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary_chr4$variant_id,split="_"))))
eqtl_Pituitary_chr4$chr<-tmp$V1
eqtl_Pituitary_chr4$loci<-tmp$V2

eqtl_Pituitary_chr4$index<-paste0(eqtl_Pituitary_chr4$chr,"-",eqtl_Pituitary_chr4$loci)
Pituitary_chr4$index<-paste0(Pituitary_chr4$chr,"-",Pituitary_chr4$pig_pos)
same_Pituitary_chr4<-intersect(Pituitary_chr4$pig_pos,eqtl_Pituitary_chr4$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary_chr5$variant_id,split="_"))))
eqtl_Pituitary_chr5$loci<-tmp$V2
same_Pituitary_chr5<-intersect(Pituitary_chr5$pig_pos,eqtl_Pituitary_chr5$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary_chr6$variant_id,split="_"))))
eqtl_Pituitary_chr6$loci<-tmp$V2
same_Pituitary_chr6<-intersect(Pituitary_chr6$pig_pos,eqtl_Pituitary_chr6$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary_chr7$variant_id,split="_"))))
eqtl_Pituitary_chr7$loci<-tmp$V2
same_Pituitary_chr7<-intersect(Pituitary_chr7$pig_pos,eqtl_Pituitary_chr7$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary_chr8$variant_id,split="_"))))
eqtl_Pituitary_chr8$loci<-tmp$V2
same_Pituitary_chr8<-intersect(Pituitary_chr8$pig_pos,eqtl_Pituitary_chr8$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary_chr9$variant_id,split="_"))))
eqtl_Pituitary_chr9$loci<-tmp$V2
same_Pituitary_chr9<-intersect(Pituitary_chr9$pig_pos,eqtl_Pituitary_chr9$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary_chr10$variant_id,split="_"))))
eqtl_Pituitary_chr10$loci<-tmp$V2
same_Pituitary_chr10<-intersect(Pituitary_chr10$pig_pos,eqtl_Pituitary_chr10$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary_chr11$variant_id,split="_"))))
eqtl_Pituitary_chr11$loci<-tmp$V2
same_Pituitary_chr11<-intersect(Pituitary_chr11$pig_pos,eqtl_Pituitary_chr11$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary_chr12$variant_id,split="_"))))
eqtl_Pituitary_chr12$loci<-tmp$V2
same_Pituitary_chr12<-intersect(Pituitary_chr12$pig_pos,eqtl_Pituitary_chr12$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary_chr13$variant_id,split="_"))))
eqtl_Pituitary_chr13$loci<-tmp$V2
same_Pituitary_chr13<-intersect(Pituitary_chr13$pig_pos,eqtl_Pituitary_chr13$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary_chr14$variant_id,split="_"))))
eqtl_Pituitary_chr14$loci<-tmp$V2
same_Pituitary_chr14<-intersect(Pituitary_chr14$pig_pos,eqtl_Pituitary_chr14$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary_chr15$variant_id,split="_"))))
eqtl_Pituitary_chr15$loci<-tmp$V2
same_Pituitary_chr15<-intersect(Pituitary_chr15$pig_pos,eqtl_Pituitary_chr15$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary_chr16$variant_id,split="_"))))
eqtl_Pituitary_chr16$loci<-tmp$V2
same_Pituitary_chr16<-intersect(Pituitary_chr16$pig_pos,eqtl_Pituitary_chr16$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary_chr17$variant_id,split="_"))))
eqtl_Pituitary_chr17$loci<-tmp$V2
same_Pituitary_chr17<-intersect(Pituitary_chr17$pig_pos,eqtl_Pituitary_chr17$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Pituitary_chr18$variant_id,split="_"))))
eqtl_Pituitary_chr18$loci<-tmp$V2
same_Pituitary_chr18<-intersect(Pituitary_chr18$pig_pos,eqtl_Pituitary_chr18$loci)

chr1_Pituitary_genes<-NULL
chr1_Pituitary_overloci<-NULL
if(length(same_Pituitary_chr1!=0)){
  for(i in 1:length(same_Pituitary_chr1)){
    a<-eqtl_Pituitary_chr1$phenotype_id[grep(same_Pituitary_chr1[i],eqtl_Pituitary_chr1$loci)]
    b<-eqtl_Pituitary_chr1[grep(same_Pituitary_chr1[i],eqtl_Pituitary_chr1$loci),]
    chr1_Pituitary_genes<-c(chr1_Pituitary_genes,a)
    chr1_Pituitary_overloci<-rbind(chr1_Pituitary_overloci,b)
  }
}

chr2_Pituitary_genes<-NULL
chr2_Pituitary_overloci<-NULL
if(length(same_Pituitary_chr2!=0)){
  for(i in 1:length(same_Pituitary_chr2)){
    a<-eqtl_Pituitary_chr2$phenotype_id[grep(same_Pituitary_chr2[i],eqtl_Pituitary_chr2$loci)]
    b<-eqtl_Pituitary_chr2[grep(same_Pituitary_chr2[i],eqtl_Pituitary_chr2$loci),]
    chr2_Pituitary_genes<-c(chr2_Pituitary_genes,a)
    chr2_Pituitary_overloci<-rbind(chr2_Pituitary_overloci,b)
  }
}

chr3_Pituitary_genes<-NULL
chr3_Pituitary_overloci<-NULL
if(length(same_Pituitary_chr3!=0)){
  for(i in 1:length(same_Pituitary_chr3)){
    a<-eqtl_Pituitary_chr3$phenotype_id[grep(same_Pituitary_chr3[i],eqtl_Pituitary_chr3$loci)]
    b<-eqtl_Pituitary_chr3[grep(same_Pituitary_chr3[i],eqtl_Pituitary_chr3$loci),]
    chr3_Pituitary_genes<-c(chr3_Pituitary_genes,a)
    chr3_Pituitary_overloci<-rbind(chr3_Pituitary_overloci,b)
  }
}

chr4_Pituitary_genes<-NULL
chr4_Pituitary_overloci<-NULL
if(length(same_Pituitary_chr4!=0)){
  for(i in 1:length(same_Pituitary_chr4)){
    a<-eqtl_Pituitary_chr4$phenotype_id[grep(same_Pituitary_chr4[i],eqtl_Pituitary_chr4$loci)]
    b<-eqtl_Pituitary_chr4[grep(same_Pituitary_chr4[i],eqtl_Pituitary_chr4$loci),]
    chr4_Pituitary_genes<-c(chr4_Pituitary_genes,a)
    chr4_Pituitary_overloci<-rbind(chr4_Pituitary_overloci,b)
  }
}

chr5_Pituitary_genes<-NULL
chr5_Pituitary_overloci<-NULL
if(length(same_Pituitary_chr5!=0)){
  for(i in 1:length(same_Pituitary_chr5)){
    a<-eqtl_Pituitary_chr5$phenotype_id[grep(same_Pituitary_chr5[i],eqtl_Pituitary_chr5$loci)]
    b<-eqtl_Pituitary_chr5[grep(same_Pituitary_chr5[i],eqtl_Pituitary_chr5$loci),]
    chr5_Pituitary_genes<-c(chr5_Pituitary_genes,a)
    chr5_Pituitary_overloci<-rbind(chr5_Pituitary_overloci,b)
  }
}

chr6_Pituitary_genes<-NULL
chr6_Pituitary_overloci<-NULL
if(length(same_Pituitary_chr6!=0)){
  for(i in 1:length(same_Pituitary_chr6)){
    a<-eqtl_Pituitary_chr6$phenotype_id[grep(same_Pituitary_chr6[i],eqtl_Pituitary_chr6$loci)]
    b<-eqtl_Pituitary_chr6[grep(same_Pituitary_chr6[i],eqtl_Pituitary_chr6$loci),]
    chr6_Pituitary_genes<-c(chr6_Pituitary_genes,a)
    chr6_Pituitary_overloci<-rbind(chr6_Pituitary_overloci,b)
  }
}

chr7_Pituitary_genes<-NULL
chr7_Pituitary_overloci<-NULL
if(length(same_Pituitary_chr7!=0)){
  for(i in 1:length(same_Pituitary_chr7)){
    a<-eqtl_Pituitary_chr7$phenotype_id[grep(same_Pituitary_chr7[i],eqtl_Pituitary_chr7$loci)]
    b<-eqtl_Pituitary_chr7[grep(same_Pituitary_chr7[i],eqtl_Pituitary_chr7$loci),]
    chr7_Pituitary_genes<-c(chr7_Pituitary_genes,a)
    chr7_Pituitary_overloci<-rbind(chr7_Pituitary_overloci,b)
  }
}

chr8_Pituitary_genes<-NULL
chr8_Pituitary_overloci<-NULL
if(length(same_Pituitary_chr8!=0)){
  for(i in 1:length(same_Pituitary_chr8)){
    a<-eqtl_Pituitary_chr8$phenotype_id[grep(same_Pituitary_chr8[i],eqtl_Pituitary_chr8$loci)]
    b<-eqtl_Pituitary_chr8[grep(same_Pituitary_chr8[i],eqtl_Pituitary_chr8$loci),]
    chr8_Pituitary_genes<-c(chr8_Pituitary_genes,a)
    chr8_Pituitary_overloci<-rbind(chr8_Pituitary_overloci,b)
  }
}

chr9_Pituitary_genes<-NULL
chr9_Pituitary_overloci<-NULL
if(length(same_Pituitary_chr9!=0)){
  for(i in 1:length(same_Pituitary_chr9)){
    a<-eqtl_Pituitary_chr9$phenotype_id[grep(same_Pituitary_chr9[i],eqtl_Pituitary_chr9$loci)]
    b<-eqtl_Pituitary_chr9[grep(same_Pituitary_chr9[i],eqtl_Pituitary_chr9$loci),]
    chr9_Pituitary_genes<-c(chr9_Pituitary_genes,a)
    chr9_Pituitary_overloci<-rbind(chr9_Pituitary_overloci,b)
  }
}

chr10_Pituitary_genes<-NULL
chr10_Pituitary_overloci<-NULL
if(length(same_Pituitary_chr10!=0)){
  for(i in 1:length(same_Pituitary_chr10)){
    a<-eqtl_Pituitary_chr10$phenotype_id[grep(same_Pituitary_chr10[i],eqtl_Pituitary_chr10$loci)]
    b<-eqtl_Pituitary_chr10[grep(same_Pituitary_chr10[i],eqtl_Pituitary_chr10$loci),]
    chr10_Pituitary_genes<-c(chr10_Pituitary_genes,a)
    chr10_Pituitary_overloci<-rbind(chr10_Pituitary_overloci,b)
  }
}

chr11_Pituitary_genes<-NULL
chr11_Pituitary_overloci<-NULL
if(length(same_Pituitary_chr11!=0)){
  for(i in 1:length(same_Pituitary_chr11)){
    a<-eqtl_Pituitary_chr11$phenotype_id[grep(same_Pituitary_chr11[i],eqtl_Pituitary_chr11$loci)]
    b<-eqtl_Pituitary_chr11[grep(same_Pituitary_chr11[i],eqtl_Pituitary_chr11$loci),]
    chr11_Pituitary_genes<-c(chr11_Pituitary_genes,a)
    chr11_Pituitary_overloci<-rbind(chr11_Pituitary_overloci,b)
  }
}

chr12_Pituitary_genes<-NULL
chr12_Pituitary_overloci<-NULL
if(length(same_Pituitary_chr12!=0)){
  for(i in 1:length(same_Pituitary_chr12)){
    a<-eqtl_Pituitary_chr12$phenotype_id[grep(same_Pituitary_chr12[i],eqtl_Pituitary_chr12$loci)]
    b<-eqtl_Pituitary_chr12[grep(same_Pituitary_chr12[i],eqtl_Pituitary_chr12$loci),]
    chr12_Pituitary_genes<-c(chr12_Pituitary_genes,a)
    chr12_Pituitary_overloci<-rbind(chr12_Pituitary_overloci,b)
  }
}

chr13_Pituitary_genes<-NULL
chr13_Pituitary_overloci<-NULL
if(length(same_Pituitary_chr13!=0)){
  for(i in 1:length(same_Pituitary_chr13)){
    a<-eqtl_Pituitary_chr13$phenotype_id[grep(same_Pituitary_chr13[i],eqtl_Pituitary_chr13$loci)]
    b<-eqtl_Pituitary_chr13[grep(same_Pituitary_chr13[i],eqtl_Pituitary_chr13$loci),]
    chr13_Pituitary_genes<-c(chr13_Pituitary_genes,a)
    chr13_Pituitary_overloci<-rbind(chr13_Pituitary_overloci,b)
  }
}

chr14_Pituitary_genes<-NULL
chr14_Pituitary_overloci<-NULL
if(length(same_Pituitary_chr14!=0)){
  for(i in 1:length(same_Pituitary_chr14)){
    a<-eqtl_Pituitary_chr14$phenotype_id[grep(same_Pituitary_chr14[i],eqtl_Pituitary_chr14$loci)]
    b<-eqtl_Pituitary_chr14[grep(same_Pituitary_chr14[i],eqtl_Pituitary_chr14$loci),]
    chr14_Pituitary_genes<-c(chr14_Pituitary_genes,a)
    chr14_Pituitary_overloci<-rbind(chr14_Pituitary_overloci,b)
  }
}

chr15_Pituitary_genes<-NULL
chr15_Pituitary_overloci<-NULL
if(length(same_Pituitary_chr15!=0)){
  for(i in 1:length(same_Pituitary_chr15)){
    a<-eqtl_Pituitary_chr15$phenotype_id[grep(same_Pituitary_chr15[i],eqtl_Pituitary_chr15$loci)]
    b<-eqtl_Pituitary_chr15[grep(same_Pituitary_chr15[i],eqtl_Pituitary_chr15$loci),]
    chr15_Pituitary_genes<-c(chr15_Pituitary_genes,a)
    chr15_Pituitary_overloci<-rbind(chr15_Pituitary_overloci,b)
  }
}

chr16_Pituitary_genes<-NULL
chr16_Pituitary_overloci<-NULL
if(length(same_Pituitary_chr16!=0)){
  for(i in 1:length(same_Pituitary_chr16)){
    a<-eqtl_Pituitary_chr16$phenotype_id[grep(same_Pituitary_chr16[i],eqtl_Pituitary_chr16$loci)]
    b<-eqtl_Pituitary_chr16[grep(same_Pituitary_chr16[i],eqtl_Pituitary_chr16$loci),]
    chr16_Pituitary_genes<-c(chr16_Pituitary_genes,a)
    chr16_Pituitary_overloci<-rbind(chr16_Pituitary_overloci,b)
  }
}

chr17_Pituitary_genes<-NULL
chr17_Pituitary_overloci<-NULL
if(length(same_Pituitary_chr17!=0)){
  for(i in 1:length(same_Pituitary_chr17)){
    a<-eqtl_Pituitary_chr17$phenotype_id[grep(same_Pituitary_chr17[i],eqtl_Pituitary_chr17$loci)]
    b<-eqtl_Pituitary_chr17[grep(same_Pituitary_chr17[i],eqtl_Pituitary_chr17$loci),]
    chr17_Pituitary_genes<-c(chr17_Pituitary_genes,a)
    chr17_Pituitary_overloci<-rbind(chr17_Pituitary_overloci,b)
  }
}

chr18_Pituitary_genes<-NULL
chr18_Pituitary_overloci<-NULL
if(length(same_Pituitary_chr18!=0)){
  for(i in 1:length(same_Pituitary_chr18)){
    a<-eqtl_Pituitary_chr18$phenotype_id[grep(same_Pituitary_chr18[i],eqtl_Pituitary_chr18$loci)]
    b<-eqtl_Pituitary_chr18[grep(same_Pituitary_chr18[i],eqtl_Pituitary_chr18$loci),]
    chr18_Pituitary_genes<-c(chr18_Pituitary_genes,a)
    chr18_Pituitary_overloci<-rbind(chr18_Pituitary_overloci,b)
  }
}
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
one2one_pig<-annotation$`Pig gene stable ID`[which(annotation$`Pig homology type`=="ortholog_one2one")]

chr1_Pituitary_overloci_one2one<-chr1_Pituitary_overloci[match(intersect(chr1_Pituitary_genes,one2one_pig),chr1_Pituitary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Pituitary_chr1$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Pituitary_chr1$hum_chr<-tmp$V1
Pituitary_chr1$hum_loci<-tmp$V3
FM_Pituitary_hum$index<-paste0(FM_Pituitary_hum$chr,"_",FM_Pituitary_hum$variant_pos)
Pituitary_chr1$index<-paste0("chr",Pituitary_chr1$hum_chr,"_",Pituitary_chr1$hum_loci)
chr1_Pituitary_hum<-NULL
for(i in 1:length(same_Pituitary_chr1)){
  a<-Pituitary_chr1[match(same_Pituitary_chr1[i],Pituitary_chr1$pig_pos),]
  b<-FM_Pituitary_hum[match(intersect(a$index,FM_Pituitary_hum$index),FM_Pituitary_hum$index),]
  chr1_Pituitary_hum<-rbind(chr1_Pituitary_hum,b)
}
chr1_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr1_Pituitary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr1_Pituitary_hum_one2one<-chr1_Pituitary_hum[match(intersect(chr1_pig2hum_one2one,chr1_Pituitary_hum$gene_id),chr1_Pituitary_hum$gene_id),]

chr2_Pituitary_overloci_one2one<-chr2_Pituitary_overloci[match(intersect(chr2_Pituitary_genes,one2one_pig),chr2_Pituitary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Pituitary_chr2$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Pituitary_chr2$hum_chr<-tmp$V1
Pituitary_chr2$hum_loci<-tmp$V3
FM_Pituitary_hum$index<-paste0(FM_Pituitary_hum$chr,"_",FM_Pituitary_hum$variant_pos)
Pituitary_chr2$index<-paste0("chr",Pituitary_chr2$hum_chr,"_",Pituitary_chr2$hum_loci)
chr2_Pituitary_hum<-NULL
for(i in 1:length(same_Pituitary_chr2)){
  a<-Pituitary_chr2[match(same_Pituitary_chr2[i],Pituitary_chr2$pig_pos),]
  b<-FM_Pituitary_hum[match(intersect(a$index,FM_Pituitary_hum$index),FM_Pituitary_hum$index),]
  chr2_Pituitary_hum<-rbind(chr2_Pituitary_hum,b)
}
chr2_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr2_Pituitary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr2_Pituitary_hum_one2one<-chr2_Pituitary_hum[match(intersect(chr2_pig2hum_one2one,chr2_Pituitary_hum$gene_id),chr2_Pituitary_hum$gene_id),]

chr3_Pituitary_overloci_one2one<-chr3_Pituitary_overloci[match(intersect(chr3_Pituitary_genes,one2one_pig),chr3_Pituitary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Pituitary_chr3$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Pituitary_chr3$hum_chr<-tmp$V1
Pituitary_chr3$hum_loci<-tmp$V3
FM_Pituitary_hum$index<-paste0(FM_Pituitary_hum$chr,"_",FM_Pituitary_hum$variant_pos)
Pituitary_chr3$index<-paste0("chr",Pituitary_chr3$hum_chr,"_",Pituitary_chr3$hum_loci)
chr3_Pituitary_hum<-NULL
for(i in 1:length(same_Pituitary_chr3)){
  a<-Pituitary_chr3[match(same_Pituitary_chr3[i],Pituitary_chr3$pig_pos),]
  b<-FM_Pituitary_hum[match(intersect(a$index,FM_Pituitary_hum$index),FM_Pituitary_hum$index),]
  chr3_Pituitary_hum<-rbind(chr3_Pituitary_hum,b)
}
chr3_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr3_Pituitary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr3_Pituitary_hum_one2one<-chr3_Pituitary_hum[match(intersect(chr3_pig2hum_one2one,chr3_Pituitary_hum$gene_id),chr3_Pituitary_hum$gene_id),]

chr4_Pituitary_overloci_one2one<-chr4_Pituitary_overloci[match(intersect(chr4_Pituitary_genes,one2one_pig),chr4_Pituitary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Pituitary_chr4$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Pituitary_chr4$hum_chr<-tmp$V1
Pituitary_chr4$hum_loci<-tmp$V3
FM_Pituitary_hum$index<-paste0(FM_Pituitary_hum$chr,"_",FM_Pituitary_hum$variant_pos)
Pituitary_chr4$index<-paste0("chr",Pituitary_chr4$hum_chr,"_",Pituitary_chr4$hum_loci)
chr4_Pituitary_hum<-NULL
for(i in 1:length(same_Pituitary_chr4)){
  a<-Pituitary_chr4[match(same_Pituitary_chr4[i],Pituitary_chr4$pig_pos),]
  b<-FM_Pituitary_hum[match(intersect(a$index,FM_Pituitary_hum$index),FM_Pituitary_hum$index),]
  chr4_Pituitary_hum<-rbind(chr4_Pituitary_hum,b)
}
chr4_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr4_Pituitary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr4_Pituitary_hum_one2one<-chr4_Pituitary_hum[match(intersect(chr4_pig2hum_one2one,chr4_Pituitary_hum$gene_id),chr4_Pituitary_hum$gene_id),]

chr5_Pituitary_overloci_one2one<-chr5_Pituitary_overloci[match(intersect(chr5_Pituitary_genes,one2one_pig),chr5_Pituitary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Pituitary_chr5$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Pituitary_chr5$hum_chr<-tmp$V1
Pituitary_chr5$hum_loci<-tmp$V3
FM_Pituitary_hum$index<-paste0(FM_Pituitary_hum$chr,"_",FM_Pituitary_hum$variant_pos)
Pituitary_chr5$index<-paste0("chr",Pituitary_chr5$hum_chr,"_",Pituitary_chr5$hum_loci)
chr5_Pituitary_hum<-NULL
for(i in 1:length(same_Pituitary_chr5)){
  a<-Pituitary_chr5[match(same_Pituitary_chr5[i],Pituitary_chr5$pig_pos),]
  b<-FM_Pituitary_hum[match(intersect(a$index,FM_Pituitary_hum$index),FM_Pituitary_hum$index),]
  chr5_Pituitary_hum<-rbind(chr5_Pituitary_hum,b)
}
chr5_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr5_Pituitary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr5_Pituitary_hum_one2one<-chr5_Pituitary_hum[match(intersect(chr5_pig2hum_one2one,chr5_Pituitary_hum$gene_id),chr5_Pituitary_hum$gene_id),]

chr6_Pituitary_overloci_one2one<-chr6_Pituitary_overloci[match(intersect(chr6_Pituitary_genes,one2one_pig),chr6_Pituitary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Pituitary_chr6$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Pituitary_chr6$hum_chr<-tmp$V1
Pituitary_chr6$hum_loci<-tmp$V3
FM_Pituitary_hum$index<-paste0(FM_Pituitary_hum$chr,"_",FM_Pituitary_hum$variant_pos)
Pituitary_chr6$index<-paste0("chr",Pituitary_chr6$hum_chr,"_",Pituitary_chr6$hum_loci)
chr6_Pituitary_hum<-NULL
for(i in 1:length(same_Pituitary_chr6)){
  a<-Pituitary_chr6[match(same_Pituitary_chr6[i],Pituitary_chr6$pig_pos),]
  b<-FM_Pituitary_hum[match(intersect(a$index,FM_Pituitary_hum$index),FM_Pituitary_hum$index),]
  chr6_Pituitary_hum<-rbind(chr6_Pituitary_hum,b)
}
chr6_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr6_Pituitary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr6_Pituitary_hum_one2one<-chr6_Pituitary_hum[match(intersect(chr6_pig2hum_one2one,chr6_Pituitary_hum$gene_id),chr6_Pituitary_hum$gene_id),]

chr7_Pituitary_overloci_one2one<-chr7_Pituitary_overloci[match(intersect(chr7_Pituitary_genes,one2one_pig),chr7_Pituitary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Pituitary_chr7$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Pituitary_chr7$hum_chr<-tmp$V1
Pituitary_chr7$hum_loci<-tmp$V3
FM_Pituitary_hum$index<-paste0(FM_Pituitary_hum$chr,"_",FM_Pituitary_hum$variant_pos)
Pituitary_chr7$index<-paste0("chr",Pituitary_chr7$hum_chr,"_",Pituitary_chr7$hum_loci)
chr7_Pituitary_hum<-NULL
for(i in 1:length(same_Pituitary_chr7)){
  a<-Pituitary_chr7[match(same_Pituitary_chr7[i],Pituitary_chr7$pig_pos),]
  b<-FM_Pituitary_hum[match(intersect(a$index,FM_Pituitary_hum$index),FM_Pituitary_hum$index),]
  chr7_Pituitary_hum<-rbind(chr7_Pituitary_hum,b)
}
chr7_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr7_Pituitary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr7_Pituitary_hum_one2one<-chr7_Pituitary_hum[match(intersect(chr7_pig2hum_one2one,chr7_Pituitary_hum$gene_id),chr7_Pituitary_hum$gene_id),]

chr8_Pituitary_overloci_one2one<-chr8_Pituitary_overloci[match(intersect(chr8_Pituitary_genes,one2one_pig),chr8_Pituitary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Pituitary_chr8$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Pituitary_chr8$hum_chr<-tmp$V1
Pituitary_chr8$hum_loci<-tmp$V3
FM_Pituitary_hum$index<-paste0(FM_Pituitary_hum$chr,"_",FM_Pituitary_hum$variant_pos)
Pituitary_chr8$index<-paste0("chr",Pituitary_chr8$hum_chr,"_",Pituitary_chr8$hum_loci)
chr8_Pituitary_hum<-NULL
for(i in 1:length(same_Pituitary_chr8)){
  a<-Pituitary_chr8[match(same_Pituitary_chr8[i],Pituitary_chr8$pig_pos),]
  b<-FM_Pituitary_hum[match(intersect(a$index,FM_Pituitary_hum$index),FM_Pituitary_hum$index),]
  chr8_Pituitary_hum<-rbind(chr8_Pituitary_hum,b)
}
chr8_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr8_Pituitary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr8_Pituitary_hum_one2one<-chr8_Pituitary_hum[match(intersect(chr8_pig2hum_one2one,chr8_Pituitary_hum$gene_id),chr8_Pituitary_hum$gene_id),]

chr9_Pituitary_overloci_one2one<-chr9_Pituitary_overloci[match(intersect(chr9_Pituitary_genes,one2one_pig),chr9_Pituitary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Pituitary_chr9$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Pituitary_chr9$hum_chr<-tmp$V1
Pituitary_chr9$hum_loci<-tmp$V3
FM_Pituitary_hum$index<-paste0(FM_Pituitary_hum$chr,"_",FM_Pituitary_hum$variant_pos)
Pituitary_chr9$index<-paste0("chr",Pituitary_chr9$hum_chr,"_",Pituitary_chr9$hum_loci)
chr9_Pituitary_hum<-NULL
for(i in 1:length(same_Pituitary_chr9)){
  a<-Pituitary_chr9[match(same_Pituitary_chr9[i],Pituitary_chr9$pig_pos),]
  b<-FM_Pituitary_hum[match(intersect(a$index,FM_Pituitary_hum$index),FM_Pituitary_hum$index),]
  chr9_Pituitary_hum<-rbind(chr9_Pituitary_hum,b)
}
chr9_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr9_Pituitary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr9_Pituitary_hum_one2one<-chr9_Pituitary_hum[match(intersect(chr9_pig2hum_one2one,chr9_Pituitary_hum$gene_id),chr9_Pituitary_hum$gene_id),]

chr10_Pituitary_overloci_one2one<-chr10_Pituitary_overloci[match(intersect(chr10_Pituitary_genes,one2one_pig),chr10_Pituitary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Pituitary_chr10$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Pituitary_chr10$hum_chr<-tmp$V1
Pituitary_chr10$hum_loci<-tmp$V3
FM_Pituitary_hum$index<-paste0(FM_Pituitary_hum$chr,"_",FM_Pituitary_hum$variant_pos)
Pituitary_chr10$index<-paste0("chr",Pituitary_chr10$hum_chr,"_",Pituitary_chr10$hum_loci)
chr10_Pituitary_hum<-NULL
for(i in 1:length(same_Pituitary_chr10)){
  a<-Pituitary_chr10[match(same_Pituitary_chr10[i],Pituitary_chr10$pig_pos),]
  b<-FM_Pituitary_hum[match(intersect(a$index,FM_Pituitary_hum$index),FM_Pituitary_hum$index),]
  chr10_Pituitary_hum<-rbind(chr10_Pituitary_hum,b)
}
chr10_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr10_Pituitary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr10_Pituitary_hum_one2one<-chr10_Pituitary_hum[match(intersect(chr10_pig2hum_one2one,chr10_Pituitary_hum$gene_id),chr10_Pituitary_hum$gene_id),]

chr11_Pituitary_overloci_one2one<-chr11_Pituitary_overloci[match(intersect(chr11_Pituitary_genes,one2one_pig),chr11_Pituitary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Pituitary_chr11$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Pituitary_chr11$hum_chr<-tmp$V1
Pituitary_chr11$hum_loci<-tmp$V3
FM_Pituitary_hum$index<-paste0(FM_Pituitary_hum$chr,"_",FM_Pituitary_hum$variant_pos)
Pituitary_chr11$index<-paste0("chr",Pituitary_chr11$hum_chr,"_",Pituitary_chr11$hum_loci)
chr11_Pituitary_hum<-NULL
for(i in 1:length(same_Pituitary_chr11)){
  a<-Pituitary_chr11[match(same_Pituitary_chr11[i],Pituitary_chr11$pig_pos),]
  b<-FM_Pituitary_hum[match(intersect(a$index,FM_Pituitary_hum$index),FM_Pituitary_hum$index),]
  chr11_Pituitary_hum<-rbind(chr11_Pituitary_hum,b)
}
chr11_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr11_Pituitary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr11_Pituitary_hum_one2one<-chr11_Pituitary_hum[match(intersect(chr11_pig2hum_one2one,chr11_Pituitary_hum$gene_id),chr11_Pituitary_hum$gene_id),]

chr12_Pituitary_overloci_one2one<-chr12_Pituitary_overloci[match(intersect(chr12_Pituitary_genes,one2one_pig),chr12_Pituitary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Pituitary_chr12$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Pituitary_chr12$hum_chr<-tmp$V1
Pituitary_chr12$hum_loci<-tmp$V3
FM_Pituitary_hum$index<-paste0(FM_Pituitary_hum$chr,"_",FM_Pituitary_hum$variant_pos)
Pituitary_chr12$index<-paste0("chr",Pituitary_chr12$hum_chr,"_",Pituitary_chr12$hum_loci)
chr12_Pituitary_hum<-NULL
for(i in 1:length(same_Pituitary_chr12)){
  a<-Pituitary_chr12[match(same_Pituitary_chr12[i],Pituitary_chr12$pig_pos),]
  b<-FM_Pituitary_hum[match(intersect(a$index,FM_Pituitary_hum$index),FM_Pituitary_hum$index),]
  chr12_Pituitary_hum<-rbind(chr12_Pituitary_hum,b)
}
chr12_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr12_Pituitary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr12_Pituitary_hum_one2one<-chr12_Pituitary_hum[match(intersect(chr12_pig2hum_one2one,chr12_Pituitary_hum$gene_id),chr12_Pituitary_hum$gene_id),]

chr13_Pituitary_overloci_one2one<-chr13_Pituitary_overloci[match(intersect(chr13_Pituitary_genes,one2one_pig),chr13_Pituitary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Pituitary_chr13$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Pituitary_chr13$hum_chr<-tmp$V1
Pituitary_chr13$hum_loci<-tmp$V3
FM_Pituitary_hum$index<-paste0(FM_Pituitary_hum$chr,"_",FM_Pituitary_hum$variant_pos)
Pituitary_chr13$index<-paste0("chr",Pituitary_chr13$hum_chr,"_",Pituitary_chr13$hum_loci)
chr13_Pituitary_hum<-NULL
for(i in 1:length(same_Pituitary_chr13)){
  a<-Pituitary_chr13[match(same_Pituitary_chr13[i],Pituitary_chr13$pig_pos),]
  b<-FM_Pituitary_hum[match(intersect(a$index,FM_Pituitary_hum$index),FM_Pituitary_hum$index),]
  chr13_Pituitary_hum<-rbind(chr13_Pituitary_hum,b)
}
chr13_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr13_Pituitary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr13_Pituitary_hum_one2one<-chr13_Pituitary_hum[match(intersect(chr13_pig2hum_one2one,chr13_Pituitary_hum$gene_id),chr13_Pituitary_hum$gene_id),]

chr14_Pituitary_overloci_one2one<-chr14_Pituitary_overloci[match(intersect(chr14_Pituitary_genes,one2one_pig),chr14_Pituitary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Pituitary_chr14$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Pituitary_chr14$hum_chr<-tmp$V1
Pituitary_chr14$hum_loci<-tmp$V3
FM_Pituitary_hum$index<-paste0(FM_Pituitary_hum$chr,"_",FM_Pituitary_hum$variant_pos)
Pituitary_chr14$index<-paste0("chr",Pituitary_chr14$hum_chr,"_",Pituitary_chr14$hum_loci)
chr14_Pituitary_hum<-NULL
for(i in 1:length(same_Pituitary_chr14)){
  a<-Pituitary_chr14[match(same_Pituitary_chr14[i],Pituitary_chr14$pig_pos),]
  b<-FM_Pituitary_hum[match(intersect(a$index,FM_Pituitary_hum$index),FM_Pituitary_hum$index),]
  chr14_Pituitary_hum<-rbind(chr14_Pituitary_hum,b)
}
chr14_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr14_Pituitary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr14_Pituitary_hum_one2one<-chr14_Pituitary_hum[match(intersect(chr14_pig2hum_one2one,chr14_Pituitary_hum$gene_id),chr14_Pituitary_hum$gene_id),]

chr15_Pituitary_overloci_one2one<-chr15_Pituitary_overloci[match(intersect(chr15_Pituitary_genes,one2one_pig),chr15_Pituitary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Pituitary_chr15$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Pituitary_chr15$hum_chr<-tmp$V1
Pituitary_chr15$hum_loci<-tmp$V3
FM_Pituitary_hum$index<-paste0(FM_Pituitary_hum$chr,"_",FM_Pituitary_hum$variant_pos)
Pituitary_chr15$index<-paste0("chr",Pituitary_chr15$hum_chr,"_",Pituitary_chr15$hum_loci)
chr15_Pituitary_hum<-NULL
for(i in 1:length(same_Pituitary_chr15)){
  a<-Pituitary_chr15[match(same_Pituitary_chr15[i],Pituitary_chr15$pig_pos),]
  b<-FM_Pituitary_hum[match(intersect(a$index,FM_Pituitary_hum$index),FM_Pituitary_hum$index),]
  chr15_Pituitary_hum<-rbind(chr15_Pituitary_hum,b)
}
chr15_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr15_Pituitary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr15_Pituitary_hum_one2one<-chr15_Pituitary_hum[match(intersect(chr15_pig2hum_one2one,chr15_Pituitary_hum$gene_id),chr15_Pituitary_hum$gene_id),]

chr16_Pituitary_overloci_one2one<-chr16_Pituitary_overloci[match(intersect(chr16_Pituitary_genes,one2one_pig),chr16_Pituitary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Pituitary_chr16$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Pituitary_chr16$hum_chr<-tmp$V1
Pituitary_chr16$hum_loci<-tmp$V3
FM_Pituitary_hum$index<-paste0(FM_Pituitary_hum$chr,"_",FM_Pituitary_hum$variant_pos)
Pituitary_chr16$index<-paste0("chr",Pituitary_chr16$hum_chr,"_",Pituitary_chr16$hum_loci)
chr16_Pituitary_hum<-NULL
for(i in 1:length(same_Pituitary_chr16)){
  a<-Pituitary_chr16[match(same_Pituitary_chr16[i],Pituitary_chr16$pig_pos),]
  b<-FM_Pituitary_hum[match(intersect(a$index,FM_Pituitary_hum$index),FM_Pituitary_hum$index),]
  chr16_Pituitary_hum<-rbind(chr16_Pituitary_hum,b)
}
chr16_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr16_Pituitary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr16_Pituitary_hum_one2one<-chr16_Pituitary_hum[match(intersect(chr16_pig2hum_one2one,chr16_Pituitary_hum$gene_id),chr16_Pituitary_hum$gene_id),]

chr17_Pituitary_overloci_one2one<-chr17_Pituitary_overloci[match(intersect(chr17_Pituitary_genes,one2one_pig),chr17_Pituitary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Pituitary_chr17$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Pituitary_chr17$hum_chr<-tmp$V1
Pituitary_chr17$hum_loci<-tmp$V3
FM_Pituitary_hum$index<-paste0(FM_Pituitary_hum$chr,"_",FM_Pituitary_hum$variant_pos)
Pituitary_chr17$index<-paste0("chr",Pituitary_chr17$hum_chr,"_",Pituitary_chr17$hum_loci)
chr17_Pituitary_hum<-NULL
for(i in 1:length(same_Pituitary_chr17)){
  a<-Pituitary_chr17[match(same_Pituitary_chr17[i],Pituitary_chr17$pig_pos),]
  b<-FM_Pituitary_hum[match(intersect(a$index,FM_Pituitary_hum$index),FM_Pituitary_hum$index),]
  chr17_Pituitary_hum<-rbind(chr17_Pituitary_hum,b)
}
chr17_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr17_Pituitary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr17_Pituitary_hum_one2one<-chr17_Pituitary_hum[match(intersect(chr17_pig2hum_one2one,chr17_Pituitary_hum$gene_id),chr17_Pituitary_hum$gene_id),]

chr18_Pituitary_overloci_one2one<-chr18_Pituitary_overloci[match(intersect(chr18_Pituitary_genes,one2one_pig),chr18_Pituitary_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Pituitary_chr18$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Pituitary_chr18$hum_chr<-tmp$V1
Pituitary_chr18$hum_loci<-tmp$V3
FM_Pituitary_hum$index<-paste0(FM_Pituitary_hum$chr,"_",FM_Pituitary_hum$variant_pos)
Pituitary_chr18$index<-paste0("chr",Pituitary_chr18$hum_chr,"_",Pituitary_chr18$hum_loci)
chr18_Pituitary_hum<-NULL
for(i in 1:length(same_Pituitary_chr18)){
  a<-Pituitary_chr18[match(same_Pituitary_chr18[i],Pituitary_chr18$pig_pos),]
  b<-FM_Pituitary_hum[match(intersect(a$index,FM_Pituitary_hum$index),FM_Pituitary_hum$index),]
  chr18_Pituitary_hum<-rbind(chr18_Pituitary_hum,b)
}
chr18_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr18_Pituitary_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr18_Pituitary_hum_one2one<-chr18_Pituitary_hum[match(intersect(chr18_pig2hum_one2one,chr18_Pituitary_hum$gene_id),chr18_Pituitary_hum$gene_id),]

Pituitary_one2one_SNP_hum<-rbind(chr1_Pituitary_hum_one2one,chr2_Pituitary_hum_one2one,chr3_Pituitary_hum_one2one,chr4_Pituitary_hum_one2one,chr5_Pituitary_hum_one2one,
                                 chr6_Pituitary_hum_one2one,chr7_Pituitary_hum_one2one,chr8_Pituitary_hum_one2one,chr9_Pituitary_hum_one2one,chr10_Pituitary_hum_one2one,
                                 chr11_Pituitary_hum_one2one,chr12_Pituitary_hum_one2one,chr13_Pituitary_hum_one2one,chr14_Pituitary_hum_one2one,chr15_Pituitary_hum_one2one,
                                 chr16_Pituitary_hum_one2one,chr17_Pituitary_hum_one2one,chr18_Pituitary_hum_one2one)

chr1_Pituitary_pig_one2one<-chr1_Pituitary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr1_Pituitary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr1_Pituitary_overloci_one2one$phenotype_id),1:9]
chr2_Pituitary_pig_one2one<-chr2_Pituitary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr2_Pituitary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr2_Pituitary_overloci_one2one$phenotype_id),1:9]
chr3_Pituitary_pig_one2one<-chr3_Pituitary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr3_Pituitary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr3_Pituitary_overloci_one2one$phenotype_id),1:9]
chr4_Pituitary_pig_one2one<-chr4_Pituitary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr4_Pituitary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr4_Pituitary_overloci_one2one$phenotype_id),1:9]
chr5_Pituitary_pig_one2one<-chr5_Pituitary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr5_Pituitary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr5_Pituitary_overloci_one2one$phenotype_id),1:9]
chr6_Pituitary_pig_one2one<-chr6_Pituitary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr6_Pituitary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr6_Pituitary_overloci_one2one$phenotype_id),1:9]
chr7_Pituitary_pig_one2one<-chr7_Pituitary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr7_Pituitary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr7_Pituitary_overloci_one2one$phenotype_id),1:9]
chr8_Pituitary_pig_one2one<-chr8_Pituitary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr8_Pituitary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr8_Pituitary_overloci_one2one$phenotype_id),1:9]
chr9_Pituitary_pig_one2one<-chr9_Pituitary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr9_Pituitary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr9_Pituitary_overloci_one2one$phenotype_id),1:9]
chr10_Pituitary_pig_one2one<-chr10_Pituitary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr10_Pituitary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr10_Pituitary_overloci_one2one$phenotype_id),1:9]
chr11_Pituitary_pig_one2one<-chr11_Pituitary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr11_Pituitary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr11_Pituitary_overloci_one2one$phenotype_id),1:9]
chr12_Pituitary_pig_one2one<-chr12_Pituitary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr12_Pituitary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr12_Pituitary_overloci_one2one$phenotype_id),1:9]
chr13_Pituitary_pig_one2one<-chr13_Pituitary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr13_Pituitary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr13_Pituitary_overloci_one2one$phenotype_id),1:9]
chr14_Pituitary_pig_one2one<-chr14_Pituitary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr14_Pituitary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr14_Pituitary_overloci_one2one$phenotype_id),1:9]
chr15_Pituitary_pig_one2one<-chr15_Pituitary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr15_Pituitary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr15_Pituitary_overloci_one2one$phenotype_id),1:9]
chr16_Pituitary_pig_one2one<-chr16_Pituitary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr16_Pituitary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr16_Pituitary_overloci_one2one$phenotype_id),1:9]
chr17_Pituitary_pig_one2one<-chr17_Pituitary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr17_Pituitary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr17_Pituitary_overloci_one2one$phenotype_id),1:9]
chr18_Pituitary_pig_one2one<-chr18_Pituitary_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr18_Pituitary_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr18_Pituitary_overloci_one2one$phenotype_id),1:9]

Pituitary_one2one_SNP_pig<-rbind(chr1_Pituitary_pig_one2one,chr2_Pituitary_pig_one2one,chr3_Pituitary_pig_one2one,chr4_Pituitary_pig_one2one,chr5_Pituitary_pig_one2one,
                                 chr6_Pituitary_pig_one2one,chr7_Pituitary_pig_one2one,chr8_Pituitary_pig_one2one,chr9_Pituitary_pig_one2one,chr10_Pituitary_pig_one2one,
                                 chr11_Pituitary_pig_one2one,chr12_Pituitary_pig_one2one,chr13_Pituitary_pig_one2one,chr14_Pituitary_pig_one2one,chr15_Pituitary_pig_one2one,
                                 chr16_Pituitary_pig_one2one,chr17_Pituitary_pig_one2one,chr18_Pituitary_pig_one2one)

Pituitary_SNP_sum<-array(NA,dim=c(nrow(Pituitary_one2one_SNP_hum),2))
colnames(Pituitary_SNP_sum)<-c("Human","Pig")
Pituitary_SNP_sum<-as.data.frame(Pituitary_SNP_sum)
Pituitary_SNP_sum$Human<-Pituitary_one2one_SNP_hum$slope / Pituitary_one2one_SNP_hum$slope_se
Pituitary_SNP_sum$Pig<-Pituitary_one2one_SNP_pig$slope / Pituitary_one2one_SNP_pig$slope_se
cor<-cor(abs(Pituitary_SNP_sum$Human),abs(Pituitary_SNP_sum$Pig))
p_val<-t.test(abs(Pituitary_SNP_sum$Human),abs(Pituitary_SNP_sum$Pig))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/Pituitary_SNP_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(Pituitary_SNP_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The number of eGenes expressed in common tissues",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

save(Pituitary_one2one_SNP_hum,Pituitary_one2one_SNP_pig,file="/Users/baizhonghao/Downloads/human-pig GTEx/Pituitary_SNP.Rdata")

#Spleen_SNP_overlaploci#
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

eqtl_Spleen_chr1<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Spleen/Spleen.cis_qtl_pairs.1.txt"))
eqtl_Spleen_chr2<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Spleen/Spleen.cis_qtl_pairs.2.txt"))
eqtl_Spleen_chr3<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Spleen/Spleen.cis_qtl_pairs.3.txt"))
eqtl_Spleen_chr4<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Spleen/Spleen.cis_qtl_pairs.4.txt"))
eqtl_Spleen_chr5<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Spleen/Spleen.cis_qtl_pairs.5.txt"))
eqtl_Spleen_chr6<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Spleen/Spleen.cis_qtl_pairs.6.txt"))
eqtl_Spleen_chr7<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Spleen/Spleen.cis_qtl_pairs.7.txt"))
eqtl_Spleen_chr8<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Spleen/Spleen.cis_qtl_pairs.8.txt"))
eqtl_Spleen_chr9<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Spleen/Spleen.cis_qtl_pairs.9.txt"))
eqtl_Spleen_chr10<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Spleen/Spleen.cis_qtl_pairs.10.txt"))
eqtl_Spleen_chr11<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Spleen/Spleen.cis_qtl_pairs.11.txt"))
eqtl_Spleen_chr12<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Spleen/Spleen.cis_qtl_pairs.12.txt"))
eqtl_Spleen_chr13<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Spleen/Spleen.cis_qtl_pairs.13.txt"))
eqtl_Spleen_chr14<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Spleen/Spleen.cis_qtl_pairs.14.txt"))
eqtl_Spleen_chr15<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Spleen/Spleen.cis_qtl_pairs.15.txt"))
eqtl_Spleen_chr16<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Spleen/Spleen.cis_qtl_pairs.16.txt"))
eqtl_Spleen_chr17<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Spleen/Spleen.cis_qtl_pairs.17.txt"))
eqtl_Spleen_chr18<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Spleen/Spleen.cis_qtl_pairs.18.txt"))

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen_chr1$variant_id,split="_"))))
eqtl_Spleen_chr1$chr<-tmp$V1
eqtl_Spleen_chr1$loci<-tmp$V2

eqtl_Spleen_chr1$index<-paste0(eqtl_Spleen_chr1$chr,"-",eqtl_Spleen_chr1$loci)
Spleen_chr1$index<-paste0(Spleen_chr1$chr,"-",Spleen_chr1$pig_pos)
same_Spleen_chr1<-intersect(Spleen_chr1$pig_pos,eqtl_Spleen_chr1$loci)


tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen_chr2$variant_id,split="_"))))
eqtl_Spleen_chr2$chr<-tmp$V1
eqtl_Spleen_chr2$loci<-tmp$V2

eqtl_Spleen_chr2$index<-paste0(eqtl_Spleen_chr2$chr,"-",eqtl_Spleen_chr2$loci)
Spleen_chr2$index<-paste0(Spleen_chr2$chr,"-",Spleen_chr2$pig_pos)
same_Spleen_chr2<-intersect(Spleen_chr2$pig_pos,eqtl_Spleen_chr2$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen_chr3$variant_id,split="_"))))
eqtl_Spleen_chr3$chr<-tmp$V1
eqtl_Spleen_chr3$loci<-tmp$V2

eqtl_Spleen_chr3$index<-paste0(eqtl_Spleen_chr3$chr,"-",eqtl_Spleen_chr3$loci)
Spleen_chr3$index<-paste0(Spleen_chr3$chr,"-",Spleen_chr3$pig_pos)
same_Spleen_chr3<-intersect(Spleen_chr3$pig_pos,eqtl_Spleen_chr3$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen_chr4$variant_id,split="_"))))
eqtl_Spleen_chr4$chr<-tmp$V1
eqtl_Spleen_chr4$loci<-tmp$V2

eqtl_Spleen_chr4$index<-paste0(eqtl_Spleen_chr4$chr,"-",eqtl_Spleen_chr4$loci)
Spleen_chr4$index<-paste0(Spleen_chr4$chr,"-",Spleen_chr4$pig_pos)
same_Spleen_chr4<-intersect(Spleen_chr4$pig_pos,eqtl_Spleen_chr4$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen_chr5$variant_id,split="_"))))
eqtl_Spleen_chr5$loci<-tmp$V2
same_Spleen_chr5<-intersect(Spleen_chr5$pig_pos,eqtl_Spleen_chr5$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen_chr6$variant_id,split="_"))))
eqtl_Spleen_chr6$loci<-tmp$V2
same_Spleen_chr6<-intersect(Spleen_chr6$pig_pos,eqtl_Spleen_chr6$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen_chr7$variant_id,split="_"))))
eqtl_Spleen_chr7$loci<-tmp$V2
same_Spleen_chr7<-intersect(Spleen_chr7$pig_pos,eqtl_Spleen_chr7$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen_chr8$variant_id,split="_"))))
eqtl_Spleen_chr8$loci<-tmp$V2
same_Spleen_chr8<-intersect(Spleen_chr8$pig_pos,eqtl_Spleen_chr8$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen_chr9$variant_id,split="_"))))
eqtl_Spleen_chr9$loci<-tmp$V2
same_Spleen_chr9<-intersect(Spleen_chr9$pig_pos,eqtl_Spleen_chr9$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen_chr10$variant_id,split="_"))))
eqtl_Spleen_chr10$loci<-tmp$V2
same_Spleen_chr10<-intersect(Spleen_chr10$pig_pos,eqtl_Spleen_chr10$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen_chr11$variant_id,split="_"))))
eqtl_Spleen_chr11$loci<-tmp$V2
same_Spleen_chr11<-intersect(Spleen_chr11$pig_pos,eqtl_Spleen_chr11$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen_chr12$variant_id,split="_"))))
eqtl_Spleen_chr12$loci<-tmp$V2
same_Spleen_chr12<-intersect(Spleen_chr12$pig_pos,eqtl_Spleen_chr12$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen_chr13$variant_id,split="_"))))
eqtl_Spleen_chr13$loci<-tmp$V2
same_Spleen_chr13<-intersect(Spleen_chr13$pig_pos,eqtl_Spleen_chr13$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen_chr14$variant_id,split="_"))))
eqtl_Spleen_chr14$loci<-tmp$V2
same_Spleen_chr14<-intersect(Spleen_chr14$pig_pos,eqtl_Spleen_chr14$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen_chr15$variant_id,split="_"))))
eqtl_Spleen_chr15$loci<-tmp$V2
same_Spleen_chr15<-intersect(Spleen_chr15$pig_pos,eqtl_Spleen_chr15$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen_chr16$variant_id,split="_"))))
eqtl_Spleen_chr16$loci<-tmp$V2
same_Spleen_chr16<-intersect(Spleen_chr16$pig_pos,eqtl_Spleen_chr16$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen_chr17$variant_id,split="_"))))
eqtl_Spleen_chr17$loci<-tmp$V2
same_Spleen_chr17<-intersect(Spleen_chr17$pig_pos,eqtl_Spleen_chr17$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Spleen_chr18$variant_id,split="_"))))
eqtl_Spleen_chr18$loci<-tmp$V2
same_Spleen_chr18<-intersect(Spleen_chr18$pig_pos,eqtl_Spleen_chr18$loci)

chr1_Spleen_genes<-NULL
chr1_Spleen_overloci<-NULL
if(length(same_Spleen_chr1!=0)){
  for(i in 1:length(same_Spleen_chr1)){
    a<-eqtl_Spleen_chr1$phenotype_id[grep(same_Spleen_chr1[i],eqtl_Spleen_chr1$loci)]
    b<-eqtl_Spleen_chr1[grep(same_Spleen_chr1[i],eqtl_Spleen_chr1$loci),]
    chr1_Spleen_genes<-c(chr1_Spleen_genes,a)
    chr1_Spleen_overloci<-rbind(chr1_Spleen_overloci,b)
  }
}

chr2_Spleen_genes<-NULL
chr2_Spleen_overloci<-NULL
if(length(same_Spleen_chr2!=0)){
  for(i in 1:length(same_Spleen_chr2)){
    a<-eqtl_Spleen_chr2$phenotype_id[grep(same_Spleen_chr2[i],eqtl_Spleen_chr2$loci)]
    b<-eqtl_Spleen_chr2[grep(same_Spleen_chr2[i],eqtl_Spleen_chr2$loci),]
    chr2_Spleen_genes<-c(chr2_Spleen_genes,a)
    chr2_Spleen_overloci<-rbind(chr2_Spleen_overloci,b)
  }
}

chr3_Spleen_genes<-NULL
chr3_Spleen_overloci<-NULL
if(length(same_Spleen_chr3!=0)){
  for(i in 1:length(same_Spleen_chr3)){
    a<-eqtl_Spleen_chr3$phenotype_id[grep(same_Spleen_chr3[i],eqtl_Spleen_chr3$loci)]
    b<-eqtl_Spleen_chr3[grep(same_Spleen_chr3[i],eqtl_Spleen_chr3$loci),]
    chr3_Spleen_genes<-c(chr3_Spleen_genes,a)
    chr3_Spleen_overloci<-rbind(chr3_Spleen_overloci,b)
  }
}

chr4_Spleen_genes<-NULL
chr4_Spleen_overloci<-NULL
if(length(same_Spleen_chr4!=0)){
  for(i in 1:length(same_Spleen_chr4)){
    a<-eqtl_Spleen_chr4$phenotype_id[grep(same_Spleen_chr4[i],eqtl_Spleen_chr4$loci)]
    b<-eqtl_Spleen_chr4[grep(same_Spleen_chr4[i],eqtl_Spleen_chr4$loci),]
    chr4_Spleen_genes<-c(chr4_Spleen_genes,a)
    chr4_Spleen_overloci<-rbind(chr4_Spleen_overloci,b)
  }
}

chr5_Spleen_genes<-NULL
chr5_Spleen_overloci<-NULL
if(length(same_Spleen_chr5!=0)){
  for(i in 1:length(same_Spleen_chr5)){
    a<-eqtl_Spleen_chr5$phenotype_id[grep(same_Spleen_chr5[i],eqtl_Spleen_chr5$loci)]
    b<-eqtl_Spleen_chr5[grep(same_Spleen_chr5[i],eqtl_Spleen_chr5$loci),]
    chr5_Spleen_genes<-c(chr5_Spleen_genes,a)
    chr5_Spleen_overloci<-rbind(chr5_Spleen_overloci,b)
  }
}

chr6_Spleen_genes<-NULL
chr6_Spleen_overloci<-NULL
if(length(same_Spleen_chr6!=0)){
  for(i in 1:length(same_Spleen_chr6)){
    a<-eqtl_Spleen_chr6$phenotype_id[grep(same_Spleen_chr6[i],eqtl_Spleen_chr6$loci)]
    b<-eqtl_Spleen_chr6[grep(same_Spleen_chr6[i],eqtl_Spleen_chr6$loci),]
    chr6_Spleen_genes<-c(chr6_Spleen_genes,a)
    chr6_Spleen_overloci<-rbind(chr6_Spleen_overloci,b)
  }
}

chr7_Spleen_genes<-NULL
chr7_Spleen_overloci<-NULL
if(length(same_Spleen_chr7!=0)){
  for(i in 1:length(same_Spleen_chr7)){
    a<-eqtl_Spleen_chr7$phenotype_id[grep(same_Spleen_chr7[i],eqtl_Spleen_chr7$loci)]
    b<-eqtl_Spleen_chr7[grep(same_Spleen_chr7[i],eqtl_Spleen_chr7$loci),]
    chr7_Spleen_genes<-c(chr7_Spleen_genes,a)
    chr7_Spleen_overloci<-rbind(chr7_Spleen_overloci,b)
  }
}

chr8_Spleen_genes<-NULL
chr8_Spleen_overloci<-NULL
if(length(same_Spleen_chr8!=0)){
  for(i in 1:length(same_Spleen_chr8)){
    a<-eqtl_Spleen_chr8$phenotype_id[grep(same_Spleen_chr8[i],eqtl_Spleen_chr8$loci)]
    b<-eqtl_Spleen_chr8[grep(same_Spleen_chr8[i],eqtl_Spleen_chr8$loci),]
    chr8_Spleen_genes<-c(chr8_Spleen_genes,a)
    chr8_Spleen_overloci<-rbind(chr8_Spleen_overloci,b)
  }
}

chr9_Spleen_genes<-NULL
chr9_Spleen_overloci<-NULL
if(length(same_Spleen_chr9!=0)){
  for(i in 1:length(same_Spleen_chr9)){
    a<-eqtl_Spleen_chr9$phenotype_id[grep(same_Spleen_chr9[i],eqtl_Spleen_chr9$loci)]
    b<-eqtl_Spleen_chr9[grep(same_Spleen_chr9[i],eqtl_Spleen_chr9$loci),]
    chr9_Spleen_genes<-c(chr9_Spleen_genes,a)
    chr9_Spleen_overloci<-rbind(chr9_Spleen_overloci,b)
  }
}

chr10_Spleen_genes<-NULL
chr10_Spleen_overloci<-NULL
if(length(same_Spleen_chr10!=0)){
  for(i in 1:length(same_Spleen_chr10)){
    a<-eqtl_Spleen_chr10$phenotype_id[grep(same_Spleen_chr10[i],eqtl_Spleen_chr10$loci)]
    b<-eqtl_Spleen_chr10[grep(same_Spleen_chr10[i],eqtl_Spleen_chr10$loci),]
    chr10_Spleen_genes<-c(chr10_Spleen_genes,a)
    chr10_Spleen_overloci<-rbind(chr10_Spleen_overloci,b)
  }
}

chr11_Spleen_genes<-NULL
chr11_Spleen_overloci<-NULL
if(length(same_Spleen_chr11!=0)){
  for(i in 1:length(same_Spleen_chr11)){
    a<-eqtl_Spleen_chr11$phenotype_id[grep(same_Spleen_chr11[i],eqtl_Spleen_chr11$loci)]
    b<-eqtl_Spleen_chr11[grep(same_Spleen_chr11[i],eqtl_Spleen_chr11$loci),]
    chr11_Spleen_genes<-c(chr11_Spleen_genes,a)
    chr11_Spleen_overloci<-rbind(chr11_Spleen_overloci,b)
  }
}

chr12_Spleen_genes<-NULL
chr12_Spleen_overloci<-NULL
if(length(same_Spleen_chr12!=0)){
  for(i in 1:length(same_Spleen_chr12)){
    a<-eqtl_Spleen_chr12$phenotype_id[grep(same_Spleen_chr12[i],eqtl_Spleen_chr12$loci)]
    b<-eqtl_Spleen_chr12[grep(same_Spleen_chr12[i],eqtl_Spleen_chr12$loci),]
    chr12_Spleen_genes<-c(chr12_Spleen_genes,a)
    chr12_Spleen_overloci<-rbind(chr12_Spleen_overloci,b)
  }
}

chr13_Spleen_genes<-NULL
chr13_Spleen_overloci<-NULL
if(length(same_Spleen_chr13!=0)){
  for(i in 1:length(same_Spleen_chr13)){
    a<-eqtl_Spleen_chr13$phenotype_id[grep(same_Spleen_chr13[i],eqtl_Spleen_chr13$loci)]
    b<-eqtl_Spleen_chr13[grep(same_Spleen_chr13[i],eqtl_Spleen_chr13$loci),]
    chr13_Spleen_genes<-c(chr13_Spleen_genes,a)
    chr13_Spleen_overloci<-rbind(chr13_Spleen_overloci,b)
  }
}

chr14_Spleen_genes<-NULL
chr14_Spleen_overloci<-NULL
if(length(same_Spleen_chr14!=0)){
  for(i in 1:length(same_Spleen_chr14)){
    a<-eqtl_Spleen_chr14$phenotype_id[grep(same_Spleen_chr14[i],eqtl_Spleen_chr14$loci)]
    b<-eqtl_Spleen_chr14[grep(same_Spleen_chr14[i],eqtl_Spleen_chr14$loci),]
    chr14_Spleen_genes<-c(chr14_Spleen_genes,a)
    chr14_Spleen_overloci<-rbind(chr14_Spleen_overloci,b)
  }
}

chr15_Spleen_genes<-NULL
chr15_Spleen_overloci<-NULL
if(length(same_Spleen_chr15!=0)){
  for(i in 1:length(same_Spleen_chr15)){
    a<-eqtl_Spleen_chr15$phenotype_id[grep(same_Spleen_chr15[i],eqtl_Spleen_chr15$loci)]
    b<-eqtl_Spleen_chr15[grep(same_Spleen_chr15[i],eqtl_Spleen_chr15$loci),]
    chr15_Spleen_genes<-c(chr15_Spleen_genes,a)
    chr15_Spleen_overloci<-rbind(chr15_Spleen_overloci,b)
  }
}

chr16_Spleen_genes<-NULL
chr16_Spleen_overloci<-NULL
if(length(same_Spleen_chr16!=0)){
  for(i in 1:length(same_Spleen_chr16)){
    a<-eqtl_Spleen_chr16$phenotype_id[grep(same_Spleen_chr16[i],eqtl_Spleen_chr16$loci)]
    b<-eqtl_Spleen_chr16[grep(same_Spleen_chr16[i],eqtl_Spleen_chr16$loci),]
    chr16_Spleen_genes<-c(chr16_Spleen_genes,a)
    chr16_Spleen_overloci<-rbind(chr16_Spleen_overloci,b)
  }
}

chr17_Spleen_genes<-NULL
chr17_Spleen_overloci<-NULL
if(length(same_Spleen_chr17!=0)){
  for(i in 1:length(same_Spleen_chr17)){
    a<-eqtl_Spleen_chr17$phenotype_id[grep(same_Spleen_chr17[i],eqtl_Spleen_chr17$loci)]
    b<-eqtl_Spleen_chr17[grep(same_Spleen_chr17[i],eqtl_Spleen_chr17$loci),]
    chr17_Spleen_genes<-c(chr17_Spleen_genes,a)
    chr17_Spleen_overloci<-rbind(chr17_Spleen_overloci,b)
  }
}

chr18_Spleen_genes<-NULL
chr18_Spleen_overloci<-NULL
if(length(same_Spleen_chr18!=0)){
  for(i in 1:length(same_Spleen_chr18)){
    a<-eqtl_Spleen_chr18$phenotype_id[grep(same_Spleen_chr18[i],eqtl_Spleen_chr18$loci)]
    b<-eqtl_Spleen_chr18[grep(same_Spleen_chr18[i],eqtl_Spleen_chr18$loci),]
    chr18_Spleen_genes<-c(chr18_Spleen_genes,a)
    chr18_Spleen_overloci<-rbind(chr18_Spleen_overloci,b)
  }
}
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
one2one_pig<-annotation$`Pig gene stable ID`[which(annotation$`Pig homology type`=="ortholog_one2one")]

chr1_Spleen_overloci_one2one<-chr1_Spleen_overloci[match(intersect(chr1_Spleen_genes,one2one_pig),chr1_Spleen_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Spleen_chr1$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Spleen_chr1$hum_chr<-tmp$V1
Spleen_chr1$hum_loci<-tmp$V3
FM_Spleen_hum$index<-paste0(FM_Spleen_hum$chr,"_",FM_Spleen_hum$variant_pos)
Spleen_chr1$index<-paste0("chr",Spleen_chr1$hum_chr,"_",Spleen_chr1$hum_loci)
chr1_Spleen_hum<-NULL
for(i in 1:length(same_Spleen_chr1)){
  a<-Spleen_chr1[match(same_Spleen_chr1[i],Spleen_chr1$pig_pos),]
  b<-FM_Spleen_hum[match(intersect(a$index,FM_Spleen_hum$index),FM_Spleen_hum$index),]
  chr1_Spleen_hum<-rbind(chr1_Spleen_hum,b)
}
chr1_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr1_Spleen_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr1_Spleen_hum_one2one<-chr1_Spleen_hum[match(intersect(chr1_pig2hum_one2one,chr1_Spleen_hum$gene_id),chr1_Spleen_hum$gene_id),]

chr2_Spleen_overloci_one2one<-chr2_Spleen_overloci[match(intersect(chr2_Spleen_genes,one2one_pig),chr2_Spleen_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Spleen_chr2$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Spleen_chr2$hum_chr<-tmp$V1
Spleen_chr2$hum_loci<-tmp$V3
FM_Spleen_hum$index<-paste0(FM_Spleen_hum$chr,"_",FM_Spleen_hum$variant_pos)
Spleen_chr2$index<-paste0("chr",Spleen_chr2$hum_chr,"_",Spleen_chr2$hum_loci)
chr2_Spleen_hum<-NULL
for(i in 1:length(same_Spleen_chr2)){
  a<-Spleen_chr2[match(same_Spleen_chr2[i],Spleen_chr2$pig_pos),]
  b<-FM_Spleen_hum[match(intersect(a$index,FM_Spleen_hum$index),FM_Spleen_hum$index),]
  chr2_Spleen_hum<-rbind(chr2_Spleen_hum,b)
}
chr2_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr2_Spleen_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr2_Spleen_hum_one2one<-chr2_Spleen_hum[match(intersect(chr2_pig2hum_one2one,chr2_Spleen_hum$gene_id),chr2_Spleen_hum$gene_id),]

chr3_Spleen_overloci_one2one<-chr3_Spleen_overloci[match(intersect(chr3_Spleen_genes,one2one_pig),chr3_Spleen_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Spleen_chr3$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Spleen_chr3$hum_chr<-tmp$V1
Spleen_chr3$hum_loci<-tmp$V3
FM_Spleen_hum$index<-paste0(FM_Spleen_hum$chr,"_",FM_Spleen_hum$variant_pos)
Spleen_chr3$index<-paste0("chr",Spleen_chr3$hum_chr,"_",Spleen_chr3$hum_loci)
chr3_Spleen_hum<-NULL
for(i in 1:length(same_Spleen_chr3)){
  a<-Spleen_chr3[match(same_Spleen_chr3[i],Spleen_chr3$pig_pos),]
  b<-FM_Spleen_hum[match(intersect(a$index,FM_Spleen_hum$index),FM_Spleen_hum$index),]
  chr3_Spleen_hum<-rbind(chr3_Spleen_hum,b)
}
chr3_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr3_Spleen_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr3_Spleen_hum_one2one<-chr3_Spleen_hum[match(intersect(chr3_pig2hum_one2one,chr3_Spleen_hum$gene_id),chr3_Spleen_hum$gene_id),]

chr4_Spleen_overloci_one2one<-chr4_Spleen_overloci[match(intersect(chr4_Spleen_genes,one2one_pig),chr4_Spleen_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Spleen_chr4$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Spleen_chr4$hum_chr<-tmp$V1
Spleen_chr4$hum_loci<-tmp$V3
FM_Spleen_hum$index<-paste0(FM_Spleen_hum$chr,"_",FM_Spleen_hum$variant_pos)
Spleen_chr4$index<-paste0("chr",Spleen_chr4$hum_chr,"_",Spleen_chr4$hum_loci)
chr4_Spleen_hum<-NULL
for(i in 1:length(same_Spleen_chr4)){
  a<-Spleen_chr4[match(same_Spleen_chr4[i],Spleen_chr4$pig_pos),]
  b<-FM_Spleen_hum[match(intersect(a$index,FM_Spleen_hum$index),FM_Spleen_hum$index),]
  chr4_Spleen_hum<-rbind(chr4_Spleen_hum,b)
}
chr4_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr4_Spleen_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr4_Spleen_hum_one2one<-chr4_Spleen_hum[match(intersect(chr4_pig2hum_one2one,chr4_Spleen_hum$gene_id),chr4_Spleen_hum$gene_id),]

chr5_Spleen_overloci_one2one<-chr5_Spleen_overloci[match(intersect(chr5_Spleen_genes,one2one_pig),chr5_Spleen_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Spleen_chr5$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Spleen_chr5$hum_chr<-tmp$V1
Spleen_chr5$hum_loci<-tmp$V3
FM_Spleen_hum$index<-paste0(FM_Spleen_hum$chr,"_",FM_Spleen_hum$variant_pos)
Spleen_chr5$index<-paste0("chr",Spleen_chr5$hum_chr,"_",Spleen_chr5$hum_loci)
chr5_Spleen_hum<-NULL
for(i in 1:length(same_Spleen_chr5)){
  a<-Spleen_chr5[match(same_Spleen_chr5[i],Spleen_chr5$pig_pos),]
  b<-FM_Spleen_hum[match(intersect(a$index,FM_Spleen_hum$index),FM_Spleen_hum$index),]
  chr5_Spleen_hum<-rbind(chr5_Spleen_hum,b)
}
chr5_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr5_Spleen_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr5_Spleen_hum_one2one<-chr5_Spleen_hum[match(intersect(chr5_pig2hum_one2one,chr5_Spleen_hum$gene_id),chr5_Spleen_hum$gene_id),]

chr6_Spleen_overloci_one2one<-chr6_Spleen_overloci[match(intersect(chr6_Spleen_genes,one2one_pig),chr6_Spleen_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Spleen_chr6$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Spleen_chr6$hum_chr<-tmp$V1
Spleen_chr6$hum_loci<-tmp$V3
FM_Spleen_hum$index<-paste0(FM_Spleen_hum$chr,"_",FM_Spleen_hum$variant_pos)
Spleen_chr6$index<-paste0("chr",Spleen_chr6$hum_chr,"_",Spleen_chr6$hum_loci)
chr6_Spleen_hum<-NULL
for(i in 1:length(same_Spleen_chr6)){
  a<-Spleen_chr6[match(same_Spleen_chr6[i],Spleen_chr6$pig_pos),]
  b<-FM_Spleen_hum[match(intersect(a$index,FM_Spleen_hum$index),FM_Spleen_hum$index),]
  chr6_Spleen_hum<-rbind(chr6_Spleen_hum,b)
}
chr6_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr6_Spleen_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr6_Spleen_hum_one2one<-chr6_Spleen_hum[match(intersect(chr6_pig2hum_one2one,chr6_Spleen_hum$gene_id),chr6_Spleen_hum$gene_id),]

chr7_Spleen_overloci_one2one<-chr7_Spleen_overloci[match(intersect(chr7_Spleen_genes,one2one_pig),chr7_Spleen_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Spleen_chr7$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Spleen_chr7$hum_chr<-tmp$V1
Spleen_chr7$hum_loci<-tmp$V3
FM_Spleen_hum$index<-paste0(FM_Spleen_hum$chr,"_",FM_Spleen_hum$variant_pos)
Spleen_chr7$index<-paste0("chr",Spleen_chr7$hum_chr,"_",Spleen_chr7$hum_loci)
chr7_Spleen_hum<-NULL
for(i in 1:length(same_Spleen_chr7)){
  a<-Spleen_chr7[match(same_Spleen_chr7[i],Spleen_chr7$pig_pos),]
  b<-FM_Spleen_hum[match(intersect(a$index,FM_Spleen_hum$index),FM_Spleen_hum$index),]
  chr7_Spleen_hum<-rbind(chr7_Spleen_hum,b)
}
chr7_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr7_Spleen_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr7_Spleen_hum_one2one<-chr7_Spleen_hum[match(intersect(chr7_pig2hum_one2one,chr7_Spleen_hum$gene_id),chr7_Spleen_hum$gene_id),]

chr8_Spleen_overloci_one2one<-chr8_Spleen_overloci[match(intersect(chr8_Spleen_genes,one2one_pig),chr8_Spleen_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Spleen_chr8$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Spleen_chr8$hum_chr<-tmp$V1
Spleen_chr8$hum_loci<-tmp$V3
FM_Spleen_hum$index<-paste0(FM_Spleen_hum$chr,"_",FM_Spleen_hum$variant_pos)
Spleen_chr8$index<-paste0("chr",Spleen_chr8$hum_chr,"_",Spleen_chr8$hum_loci)
chr8_Spleen_hum<-NULL
for(i in 1:length(same_Spleen_chr8)){
  a<-Spleen_chr8[match(same_Spleen_chr8[i],Spleen_chr8$pig_pos),]
  b<-FM_Spleen_hum[match(intersect(a$index,FM_Spleen_hum$index),FM_Spleen_hum$index),]
  chr8_Spleen_hum<-rbind(chr8_Spleen_hum,b)
}
chr8_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr8_Spleen_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr8_Spleen_hum_one2one<-chr8_Spleen_hum[match(intersect(chr8_pig2hum_one2one,chr8_Spleen_hum$gene_id),chr8_Spleen_hum$gene_id),]

chr9_Spleen_overloci_one2one<-chr9_Spleen_overloci[match(intersect(chr9_Spleen_genes,one2one_pig),chr9_Spleen_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Spleen_chr9$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Spleen_chr9$hum_chr<-tmp$V1
Spleen_chr9$hum_loci<-tmp$V3
FM_Spleen_hum$index<-paste0(FM_Spleen_hum$chr,"_",FM_Spleen_hum$variant_pos)
Spleen_chr9$index<-paste0("chr",Spleen_chr9$hum_chr,"_",Spleen_chr9$hum_loci)
chr9_Spleen_hum<-NULL
for(i in 1:length(same_Spleen_chr9)){
  a<-Spleen_chr9[match(same_Spleen_chr9[i],Spleen_chr9$pig_pos),]
  b<-FM_Spleen_hum[match(intersect(a$index,FM_Spleen_hum$index),FM_Spleen_hum$index),]
  chr9_Spleen_hum<-rbind(chr9_Spleen_hum,b)
}
chr9_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr9_Spleen_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr9_Spleen_hum_one2one<-chr9_Spleen_hum[match(intersect(chr9_pig2hum_one2one,chr9_Spleen_hum$gene_id),chr9_Spleen_hum$gene_id),]

chr10_Spleen_overloci_one2one<-chr10_Spleen_overloci[match(intersect(chr10_Spleen_genes,one2one_pig),chr10_Spleen_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Spleen_chr10$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Spleen_chr10$hum_chr<-tmp$V1
Spleen_chr10$hum_loci<-tmp$V3
FM_Spleen_hum$index<-paste0(FM_Spleen_hum$chr,"_",FM_Spleen_hum$variant_pos)
Spleen_chr10$index<-paste0("chr",Spleen_chr10$hum_chr,"_",Spleen_chr10$hum_loci)
chr10_Spleen_hum<-NULL
for(i in 1:length(same_Spleen_chr10)){
  a<-Spleen_chr10[match(same_Spleen_chr10[i],Spleen_chr10$pig_pos),]
  b<-FM_Spleen_hum[match(intersect(a$index,FM_Spleen_hum$index),FM_Spleen_hum$index),]
  chr10_Spleen_hum<-rbind(chr10_Spleen_hum,b)
}
chr10_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr10_Spleen_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr10_Spleen_hum_one2one<-chr10_Spleen_hum[match(intersect(chr10_pig2hum_one2one,chr10_Spleen_hum$gene_id),chr10_Spleen_hum$gene_id),]

chr11_Spleen_overloci_one2one<-chr11_Spleen_overloci[match(intersect(chr11_Spleen_genes,one2one_pig),chr11_Spleen_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Spleen_chr11$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Spleen_chr11$hum_chr<-tmp$V1
Spleen_chr11$hum_loci<-tmp$V3
FM_Spleen_hum$index<-paste0(FM_Spleen_hum$chr,"_",FM_Spleen_hum$variant_pos)
Spleen_chr11$index<-paste0("chr",Spleen_chr11$hum_chr,"_",Spleen_chr11$hum_loci)
chr11_Spleen_hum<-NULL
for(i in 1:length(same_Spleen_chr11)){
  a<-Spleen_chr11[match(same_Spleen_chr11[i],Spleen_chr11$pig_pos),]
  b<-FM_Spleen_hum[match(intersect(a$index,FM_Spleen_hum$index),FM_Spleen_hum$index),]
  chr11_Spleen_hum<-rbind(chr11_Spleen_hum,b)
}
chr11_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr11_Spleen_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr11_Spleen_hum_one2one<-chr11_Spleen_hum[match(intersect(chr11_pig2hum_one2one,chr11_Spleen_hum$gene_id),chr11_Spleen_hum$gene_id),]

chr12_Spleen_overloci_one2one<-chr12_Spleen_overloci[match(intersect(chr12_Spleen_genes,one2one_pig),chr12_Spleen_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Spleen_chr12$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Spleen_chr12$hum_chr<-tmp$V1
Spleen_chr12$hum_loci<-tmp$V3
FM_Spleen_hum$index<-paste0(FM_Spleen_hum$chr,"_",FM_Spleen_hum$variant_pos)
Spleen_chr12$index<-paste0("chr",Spleen_chr12$hum_chr,"_",Spleen_chr12$hum_loci)
chr12_Spleen_hum<-NULL
for(i in 1:length(same_Spleen_chr12)){
  a<-Spleen_chr12[match(same_Spleen_chr12[i],Spleen_chr12$pig_pos),]
  b<-FM_Spleen_hum[match(intersect(a$index,FM_Spleen_hum$index),FM_Spleen_hum$index),]
  chr12_Spleen_hum<-rbind(chr12_Spleen_hum,b)
}
chr12_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr12_Spleen_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr12_Spleen_hum_one2one<-chr12_Spleen_hum[match(intersect(chr12_pig2hum_one2one,chr12_Spleen_hum$gene_id),chr12_Spleen_hum$gene_id),]

chr13_Spleen_overloci_one2one<-chr13_Spleen_overloci[match(intersect(chr13_Spleen_genes,one2one_pig),chr13_Spleen_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Spleen_chr13$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Spleen_chr13$hum_chr<-tmp$V1
Spleen_chr13$hum_loci<-tmp$V3
FM_Spleen_hum$index<-paste0(FM_Spleen_hum$chr,"_",FM_Spleen_hum$variant_pos)
Spleen_chr13$index<-paste0("chr",Spleen_chr13$hum_chr,"_",Spleen_chr13$hum_loci)
chr13_Spleen_hum<-NULL
for(i in 1:length(same_Spleen_chr13)){
  a<-Spleen_chr13[match(same_Spleen_chr13[i],Spleen_chr13$pig_pos),]
  b<-FM_Spleen_hum[match(intersect(a$index,FM_Spleen_hum$index),FM_Spleen_hum$index),]
  chr13_Spleen_hum<-rbind(chr13_Spleen_hum,b)
}
chr13_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr13_Spleen_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr13_Spleen_hum_one2one<-chr13_Spleen_hum[match(intersect(chr13_pig2hum_one2one,chr13_Spleen_hum$gene_id),chr13_Spleen_hum$gene_id),]

chr14_Spleen_overloci_one2one<-chr14_Spleen_overloci[match(intersect(chr14_Spleen_genes,one2one_pig),chr14_Spleen_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Spleen_chr14$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Spleen_chr14$hum_chr<-tmp$V1
Spleen_chr14$hum_loci<-tmp$V3
FM_Spleen_hum$index<-paste0(FM_Spleen_hum$chr,"_",FM_Spleen_hum$variant_pos)
Spleen_chr14$index<-paste0("chr",Spleen_chr14$hum_chr,"_",Spleen_chr14$hum_loci)
chr14_Spleen_hum<-NULL
for(i in 1:length(same_Spleen_chr14)){
  a<-Spleen_chr14[match(same_Spleen_chr14[i],Spleen_chr14$pig_pos),]
  b<-FM_Spleen_hum[match(intersect(a$index,FM_Spleen_hum$index),FM_Spleen_hum$index),]
  chr14_Spleen_hum<-rbind(chr14_Spleen_hum,b)
}
chr14_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr14_Spleen_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr14_Spleen_hum_one2one<-chr14_Spleen_hum[match(intersect(chr14_pig2hum_one2one,chr14_Spleen_hum$gene_id),chr14_Spleen_hum$gene_id),]

chr15_Spleen_overloci_one2one<-chr15_Spleen_overloci[match(intersect(chr15_Spleen_genes,one2one_pig),chr15_Spleen_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Spleen_chr15$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Spleen_chr15$hum_chr<-tmp$V1
Spleen_chr15$hum_loci<-tmp$V3
FM_Spleen_hum$index<-paste0(FM_Spleen_hum$chr,"_",FM_Spleen_hum$variant_pos)
Spleen_chr15$index<-paste0("chr",Spleen_chr15$hum_chr,"_",Spleen_chr15$hum_loci)
chr15_Spleen_hum<-NULL
for(i in 1:length(same_Spleen_chr15)){
  a<-Spleen_chr15[match(same_Spleen_chr15[i],Spleen_chr15$pig_pos),]
  b<-FM_Spleen_hum[match(intersect(a$index,FM_Spleen_hum$index),FM_Spleen_hum$index),]
  chr15_Spleen_hum<-rbind(chr15_Spleen_hum,b)
}
chr15_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr15_Spleen_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr15_Spleen_hum_one2one<-chr15_Spleen_hum[match(intersect(chr15_pig2hum_one2one,chr15_Spleen_hum$gene_id),chr15_Spleen_hum$gene_id),]

chr16_Spleen_overloci_one2one<-chr16_Spleen_overloci[match(intersect(chr16_Spleen_genes,one2one_pig),chr16_Spleen_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Spleen_chr16$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Spleen_chr16$hum_chr<-tmp$V1
Spleen_chr16$hum_loci<-tmp$V3
FM_Spleen_hum$index<-paste0(FM_Spleen_hum$chr,"_",FM_Spleen_hum$variant_pos)
Spleen_chr16$index<-paste0("chr",Spleen_chr16$hum_chr,"_",Spleen_chr16$hum_loci)
chr16_Spleen_hum<-NULL
for(i in 1:length(same_Spleen_chr16)){
  a<-Spleen_chr16[match(same_Spleen_chr16[i],Spleen_chr16$pig_pos),]
  b<-FM_Spleen_hum[match(intersect(a$index,FM_Spleen_hum$index),FM_Spleen_hum$index),]
  chr16_Spleen_hum<-rbind(chr16_Spleen_hum,b)
}
chr16_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr16_Spleen_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr16_Spleen_hum_one2one<-chr16_Spleen_hum[match(intersect(chr16_pig2hum_one2one,chr16_Spleen_hum$gene_id),chr16_Spleen_hum$gene_id),]

chr17_Spleen_overloci_one2one<-chr17_Spleen_overloci[match(intersect(chr17_Spleen_genes,one2one_pig),chr17_Spleen_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Spleen_chr17$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Spleen_chr17$hum_chr<-tmp$V1
Spleen_chr17$hum_loci<-tmp$V3
FM_Spleen_hum$index<-paste0(FM_Spleen_hum$chr,"_",FM_Spleen_hum$variant_pos)
Spleen_chr17$index<-paste0("chr",Spleen_chr17$hum_chr,"_",Spleen_chr17$hum_loci)
chr17_Spleen_hum<-NULL
for(i in 1:length(same_Spleen_chr17)){
  a<-Spleen_chr17[match(same_Spleen_chr17[i],Spleen_chr17$pig_pos),]
  b<-FM_Spleen_hum[match(intersect(a$index,FM_Spleen_hum$index),FM_Spleen_hum$index),]
  chr17_Spleen_hum<-rbind(chr17_Spleen_hum,b)
}
chr17_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr17_Spleen_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr17_Spleen_hum_one2one<-chr17_Spleen_hum[match(intersect(chr17_pig2hum_one2one,chr17_Spleen_hum$gene_id),chr17_Spleen_hum$gene_id),]

chr18_Spleen_overloci_one2one<-chr18_Spleen_overloci[match(intersect(chr18_Spleen_genes,one2one_pig),chr18_Spleen_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Spleen_chr18$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Spleen_chr18$hum_chr<-tmp$V1
Spleen_chr18$hum_loci<-tmp$V3
FM_Spleen_hum$index<-paste0(FM_Spleen_hum$chr,"_",FM_Spleen_hum$variant_pos)
Spleen_chr18$index<-paste0("chr",Spleen_chr18$hum_chr,"_",Spleen_chr18$hum_loci)
chr18_Spleen_hum<-NULL
for(i in 1:length(same_Spleen_chr18)){
  a<-Spleen_chr18[match(same_Spleen_chr18[i],Spleen_chr18$pig_pos),]
  b<-FM_Spleen_hum[match(intersect(a$index,FM_Spleen_hum$index),FM_Spleen_hum$index),]
  chr18_Spleen_hum<-rbind(chr18_Spleen_hum,b)
}
chr18_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr18_Spleen_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr18_Spleen_hum_one2one<-chr18_Spleen_hum[match(intersect(chr18_pig2hum_one2one,chr18_Spleen_hum$gene_id),chr18_Spleen_hum$gene_id),]

Spleen_one2one_SNP_hum<-rbind(chr1_Spleen_hum_one2one,chr2_Spleen_hum_one2one,chr3_Spleen_hum_one2one,chr4_Spleen_hum_one2one,chr5_Spleen_hum_one2one,
                              chr6_Spleen_hum_one2one,chr7_Spleen_hum_one2one,chr8_Spleen_hum_one2one,chr9_Spleen_hum_one2one,chr10_Spleen_hum_one2one,
                              chr11_Spleen_hum_one2one,chr12_Spleen_hum_one2one,chr13_Spleen_hum_one2one,chr14_Spleen_hum_one2one,chr15_Spleen_hum_one2one,
                              chr16_Spleen_hum_one2one,chr17_Spleen_hum_one2one,chr18_Spleen_hum_one2one)

chr1_Spleen_pig_one2one<-chr1_Spleen_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr1_Spleen_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr1_Spleen_overloci_one2one$phenotype_id),1:9]
chr2_Spleen_pig_one2one<-chr2_Spleen_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr2_Spleen_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr2_Spleen_overloci_one2one$phenotype_id),1:9]
chr3_Spleen_pig_one2one<-chr3_Spleen_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr3_Spleen_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr3_Spleen_overloci_one2one$phenotype_id),1:9]
chr4_Spleen_pig_one2one<-chr4_Spleen_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr4_Spleen_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr4_Spleen_overloci_one2one$phenotype_id),1:9]
chr5_Spleen_pig_one2one<-chr5_Spleen_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr5_Spleen_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr5_Spleen_overloci_one2one$phenotype_id),1:9]
chr6_Spleen_pig_one2one<-chr6_Spleen_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr6_Spleen_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr6_Spleen_overloci_one2one$phenotype_id),1:9]
chr7_Spleen_pig_one2one<-chr7_Spleen_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr7_Spleen_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr7_Spleen_overloci_one2one$phenotype_id),1:9]
chr8_Spleen_pig_one2one<-chr8_Spleen_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr8_Spleen_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr8_Spleen_overloci_one2one$phenotype_id),1:9]
chr9_Spleen_pig_one2one<-chr9_Spleen_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr9_Spleen_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr9_Spleen_overloci_one2one$phenotype_id),1:9]
chr10_Spleen_pig_one2one<-chr10_Spleen_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr10_Spleen_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr10_Spleen_overloci_one2one$phenotype_id),1:9]
chr11_Spleen_pig_one2one<-chr11_Spleen_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr11_Spleen_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr11_Spleen_overloci_one2one$phenotype_id),1:9]
chr12_Spleen_pig_one2one<-chr12_Spleen_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr12_Spleen_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr12_Spleen_overloci_one2one$phenotype_id),1:9]
chr13_Spleen_pig_one2one<-chr13_Spleen_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr13_Spleen_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr13_Spleen_overloci_one2one$phenotype_id),1:9]
chr14_Spleen_pig_one2one<-chr14_Spleen_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr14_Spleen_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr14_Spleen_overloci_one2one$phenotype_id),1:9]
chr15_Spleen_pig_one2one<-chr15_Spleen_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr15_Spleen_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr15_Spleen_overloci_one2one$phenotype_id),1:9]
chr16_Spleen_pig_one2one<-chr16_Spleen_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr16_Spleen_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr16_Spleen_overloci_one2one$phenotype_id),1:9]
chr17_Spleen_pig_one2one<-chr17_Spleen_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr17_Spleen_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr17_Spleen_overloci_one2one$phenotype_id),1:9]
chr18_Spleen_pig_one2one<-chr18_Spleen_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr18_Spleen_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr18_Spleen_overloci_one2one$phenotype_id),1:9]

Spleen_one2one_SNP_pig<-rbind(chr1_Spleen_pig_one2one,chr2_Spleen_pig_one2one,chr3_Spleen_pig_one2one,chr4_Spleen_pig_one2one,chr5_Spleen_pig_one2one,
                              chr6_Spleen_pig_one2one,chr7_Spleen_pig_one2one,chr8_Spleen_pig_one2one,chr9_Spleen_pig_one2one,chr10_Spleen_pig_one2one,
                              chr11_Spleen_pig_one2one,chr12_Spleen_pig_one2one,chr13_Spleen_pig_one2one,chr14_Spleen_pig_one2one,chr15_Spleen_pig_one2one,
                              chr16_Spleen_pig_one2one,chr17_Spleen_pig_one2one,chr18_Spleen_pig_one2one)

Spleen_SNP_sum<-array(NA,dim=c(nrow(Spleen_one2one_SNP_hum),2))
colnames(Spleen_SNP_sum)<-c("Human","Pig")
Spleen_SNP_sum<-as.data.frame(Spleen_SNP_sum)
Spleen_SNP_sum$Human<-Spleen_one2one_SNP_hum$slope / Spleen_one2one_SNP_hum$slope_se
Spleen_SNP_sum$Pig<-Spleen_one2one_SNP_pig$slope / Spleen_one2one_SNP_pig$slope_se
cor<-cor(abs(Spleen_SNP_sum$Human),abs(Spleen_SNP_sum$Pig))
p_val<-t.test(abs(Spleen_SNP_sum$Human),abs(Spleen_SNP_sum$Pig))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/Spleen_SNP_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(Spleen_SNP_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The number of eGenes expressed in common tissues",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

save(Spleen_one2one_SNP_hum,Spleen_one2one_SNP_pig,file="/Users/baizhonghao/Downloads/human-pig GTEx/Spleen_SNP.Rdata")

#Testis_SNP_overlaploci#
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

eqtl_Testis_chr1<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Testis/Testis.cis_qtl_pairs.1.txt"))
eqtl_Testis_chr2<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Testis/Testis.cis_qtl_pairs.2.txt"))
eqtl_Testis_chr3<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Testis/Testis.cis_qtl_pairs.3.txt"))
eqtl_Testis_chr4<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Testis/Testis.cis_qtl_pairs.4.txt"))
eqtl_Testis_chr5<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Testis/Testis.cis_qtl_pairs.5.txt"))
eqtl_Testis_chr6<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Testis/Testis.cis_qtl_pairs.6.txt"))
eqtl_Testis_chr7<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Testis/Testis.cis_qtl_pairs.7.txt"))
eqtl_Testis_chr8<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Testis/Testis.cis_qtl_pairs.8.txt"))
eqtl_Testis_chr9<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Testis/Testis.cis_qtl_pairs.9.txt"))
eqtl_Testis_chr10<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Testis/Testis.cis_qtl_pairs.10.txt"))
eqtl_Testis_chr11<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Testis/Testis.cis_qtl_pairs.11.txt"))
eqtl_Testis_chr12<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Testis/Testis.cis_qtl_pairs.12.txt"))
eqtl_Testis_chr13<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Testis/Testis.cis_qtl_pairs.13.txt"))
eqtl_Testis_chr14<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Testis/Testis.cis_qtl_pairs.14.txt"))
eqtl_Testis_chr15<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Testis/Testis.cis_qtl_pairs.15.txt"))
eqtl_Testis_chr16<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Testis/Testis.cis_qtl_pairs.16.txt"))
eqtl_Testis_chr17<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Testis/Testis.cis_qtl_pairs.17.txt"))
eqtl_Testis_chr18<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Testis/Testis.cis_qtl_pairs.18.txt"))

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis_chr1$variant_id,split="_"))))
eqtl_Testis_chr1$chr<-tmp$V1
eqtl_Testis_chr1$loci<-tmp$V2

eqtl_Testis_chr1$index<-paste0(eqtl_Testis_chr1$chr,"-",eqtl_Testis_chr1$loci)
Testis_chr1$index<-paste0(Testis_chr1$chr,"-",Testis_chr1$pig_pos)
same_Testis_chr1<-intersect(Testis_chr1$pig_pos,eqtl_Testis_chr1$loci)


tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis_chr2$variant_id,split="_"))))
eqtl_Testis_chr2$chr<-tmp$V1
eqtl_Testis_chr2$loci<-tmp$V2

eqtl_Testis_chr2$index<-paste0(eqtl_Testis_chr2$chr,"-",eqtl_Testis_chr2$loci)
Testis_chr2$index<-paste0(Testis_chr2$chr,"-",Testis_chr2$pig_pos)
same_Testis_chr2<-intersect(Testis_chr2$pig_pos,eqtl_Testis_chr2$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis_chr3$variant_id,split="_"))))
eqtl_Testis_chr3$chr<-tmp$V1
eqtl_Testis_chr3$loci<-tmp$V2

eqtl_Testis_chr3$index<-paste0(eqtl_Testis_chr3$chr,"-",eqtl_Testis_chr3$loci)
Testis_chr3$index<-paste0(Testis_chr3$chr,"-",Testis_chr3$pig_pos)
same_Testis_chr3<-intersect(Testis_chr3$pig_pos,eqtl_Testis_chr3$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis_chr4$variant_id,split="_"))))
eqtl_Testis_chr4$chr<-tmp$V1
eqtl_Testis_chr4$loci<-tmp$V2

eqtl_Testis_chr4$index<-paste0(eqtl_Testis_chr4$chr,"-",eqtl_Testis_chr4$loci)
Testis_chr4$index<-paste0(Testis_chr4$chr,"-",Testis_chr4$pig_pos)
same_Testis_chr4<-intersect(Testis_chr4$pig_pos,eqtl_Testis_chr4$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis_chr5$variant_id,split="_"))))
eqtl_Testis_chr5$loci<-tmp$V2
same_Testis_chr5<-intersect(Testis_chr5$pig_pos,eqtl_Testis_chr5$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis_chr6$variant_id,split="_"))))
eqtl_Testis_chr6$loci<-tmp$V2
same_Testis_chr6<-intersect(Testis_chr6$pig_pos,eqtl_Testis_chr6$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis_chr7$variant_id,split="_"))))
eqtl_Testis_chr7$loci<-tmp$V2
same_Testis_chr7<-intersect(Testis_chr7$pig_pos,eqtl_Testis_chr7$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis_chr8$variant_id,split="_"))))
eqtl_Testis_chr8$loci<-tmp$V2
same_Testis_chr8<-intersect(Testis_chr8$pig_pos,eqtl_Testis_chr8$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis_chr9$variant_id,split="_"))))
eqtl_Testis_chr9$loci<-tmp$V2
same_Testis_chr9<-intersect(Testis_chr9$pig_pos,eqtl_Testis_chr9$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis_chr10$variant_id,split="_"))))
eqtl_Testis_chr10$loci<-tmp$V2
same_Testis_chr10<-intersect(Testis_chr10$pig_pos,eqtl_Testis_chr10$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis_chr11$variant_id,split="_"))))
eqtl_Testis_chr11$loci<-tmp$V2
same_Testis_chr11<-intersect(Testis_chr11$pig_pos,eqtl_Testis_chr11$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis_chr12$variant_id,split="_"))))
eqtl_Testis_chr12$loci<-tmp$V2
same_Testis_chr12<-intersect(Testis_chr12$pig_pos,eqtl_Testis_chr12$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis_chr13$variant_id,split="_"))))
eqtl_Testis_chr13$loci<-tmp$V2
same_Testis_chr13<-intersect(Testis_chr13$pig_pos,eqtl_Testis_chr13$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis_chr14$variant_id,split="_"))))
eqtl_Testis_chr14$loci<-tmp$V2
same_Testis_chr14<-intersect(Testis_chr14$pig_pos,eqtl_Testis_chr14$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis_chr15$variant_id,split="_"))))
eqtl_Testis_chr15$loci<-tmp$V2
same_Testis_chr15<-intersect(Testis_chr15$pig_pos,eqtl_Testis_chr15$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis_chr16$variant_id,split="_"))))
eqtl_Testis_chr16$loci<-tmp$V2
same_Testis_chr16<-intersect(Testis_chr16$pig_pos,eqtl_Testis_chr16$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis_chr17$variant_id,split="_"))))
eqtl_Testis_chr17$loci<-tmp$V2
same_Testis_chr17<-intersect(Testis_chr17$pig_pos,eqtl_Testis_chr17$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Testis_chr18$variant_id,split="_"))))
eqtl_Testis_chr18$loci<-tmp$V2
same_Testis_chr18<-intersect(Testis_chr18$pig_pos,eqtl_Testis_chr18$loci)

chr1_Testis_genes<-NULL
chr1_Testis_overloci<-NULL
if(length(same_Testis_chr1!=0)){
  for(i in 1:length(same_Testis_chr1)){
    a<-eqtl_Testis_chr1$phenotype_id[grep(same_Testis_chr1[i],eqtl_Testis_chr1$loci)]
    b<-eqtl_Testis_chr1[grep(same_Testis_chr1[i],eqtl_Testis_chr1$loci),]
    chr1_Testis_genes<-c(chr1_Testis_genes,a)
    chr1_Testis_overloci<-rbind(chr1_Testis_overloci,b)
  }
}

chr2_Testis_genes<-NULL
chr2_Testis_overloci<-NULL
if(length(same_Testis_chr2!=0)){
  for(i in 1:length(same_Testis_chr2)){
    a<-eqtl_Testis_chr2$phenotype_id[grep(same_Testis_chr2[i],eqtl_Testis_chr2$loci)]
    b<-eqtl_Testis_chr2[grep(same_Testis_chr2[i],eqtl_Testis_chr2$loci),]
    chr2_Testis_genes<-c(chr2_Testis_genes,a)
    chr2_Testis_overloci<-rbind(chr2_Testis_overloci,b)
  }
}

chr3_Testis_genes<-NULL
chr3_Testis_overloci<-NULL
if(length(same_Testis_chr3!=0)){
  for(i in 1:length(same_Testis_chr3)){
    a<-eqtl_Testis_chr3$phenotype_id[grep(same_Testis_chr3[i],eqtl_Testis_chr3$loci)]
    b<-eqtl_Testis_chr3[grep(same_Testis_chr3[i],eqtl_Testis_chr3$loci),]
    chr3_Testis_genes<-c(chr3_Testis_genes,a)
    chr3_Testis_overloci<-rbind(chr3_Testis_overloci,b)
  }
}

chr4_Testis_genes<-NULL
chr4_Testis_overloci<-NULL
if(length(same_Testis_chr4!=0)){
  for(i in 1:length(same_Testis_chr4)){
    a<-eqtl_Testis_chr4$phenotype_id[grep(same_Testis_chr4[i],eqtl_Testis_chr4$loci)]
    b<-eqtl_Testis_chr4[grep(same_Testis_chr4[i],eqtl_Testis_chr4$loci),]
    chr4_Testis_genes<-c(chr4_Testis_genes,a)
    chr4_Testis_overloci<-rbind(chr4_Testis_overloci,b)
  }
}

chr5_Testis_genes<-NULL
chr5_Testis_overloci<-NULL
if(length(same_Testis_chr5!=0)){
  for(i in 1:length(same_Testis_chr5)){
    a<-eqtl_Testis_chr5$phenotype_id[grep(same_Testis_chr5[i],eqtl_Testis_chr5$loci)]
    b<-eqtl_Testis_chr5[grep(same_Testis_chr5[i],eqtl_Testis_chr5$loci),]
    chr5_Testis_genes<-c(chr5_Testis_genes,a)
    chr5_Testis_overloci<-rbind(chr5_Testis_overloci,b)
  }
}

chr6_Testis_genes<-NULL
chr6_Testis_overloci<-NULL
if(length(same_Testis_chr6!=0)){
  for(i in 1:length(same_Testis_chr6)){
    a<-eqtl_Testis_chr6$phenotype_id[grep(same_Testis_chr6[i],eqtl_Testis_chr6$loci)]
    b<-eqtl_Testis_chr6[grep(same_Testis_chr6[i],eqtl_Testis_chr6$loci),]
    chr6_Testis_genes<-c(chr6_Testis_genes,a)
    chr6_Testis_overloci<-rbind(chr6_Testis_overloci,b)
  }
}

chr7_Testis_genes<-NULL
chr7_Testis_overloci<-NULL
if(length(same_Testis_chr7!=0)){
  for(i in 1:length(same_Testis_chr7)){
    a<-eqtl_Testis_chr7$phenotype_id[grep(same_Testis_chr7[i],eqtl_Testis_chr7$loci)]
    b<-eqtl_Testis_chr7[grep(same_Testis_chr7[i],eqtl_Testis_chr7$loci),]
    chr7_Testis_genes<-c(chr7_Testis_genes,a)
    chr7_Testis_overloci<-rbind(chr7_Testis_overloci,b)
  }
}

chr8_Testis_genes<-NULL
chr8_Testis_overloci<-NULL
if(length(same_Testis_chr8!=0)){
  for(i in 1:length(same_Testis_chr8)){
    a<-eqtl_Testis_chr8$phenotype_id[grep(same_Testis_chr8[i],eqtl_Testis_chr8$loci)]
    b<-eqtl_Testis_chr8[grep(same_Testis_chr8[i],eqtl_Testis_chr8$loci),]
    chr8_Testis_genes<-c(chr8_Testis_genes,a)
    chr8_Testis_overloci<-rbind(chr8_Testis_overloci,b)
  }
}

chr9_Testis_genes<-NULL
chr9_Testis_overloci<-NULL
if(length(same_Testis_chr9!=0)){
  for(i in 1:length(same_Testis_chr9)){
    a<-eqtl_Testis_chr9$phenotype_id[grep(same_Testis_chr9[i],eqtl_Testis_chr9$loci)]
    b<-eqtl_Testis_chr9[grep(same_Testis_chr9[i],eqtl_Testis_chr9$loci),]
    chr9_Testis_genes<-c(chr9_Testis_genes,a)
    chr9_Testis_overloci<-rbind(chr9_Testis_overloci,b)
  }
}

chr10_Testis_genes<-NULL
chr10_Testis_overloci<-NULL
if(length(same_Testis_chr10!=0)){
  for(i in 1:length(same_Testis_chr10)){
    a<-eqtl_Testis_chr10$phenotype_id[grep(same_Testis_chr10[i],eqtl_Testis_chr10$loci)]
    b<-eqtl_Testis_chr10[grep(same_Testis_chr10[i],eqtl_Testis_chr10$loci),]
    chr10_Testis_genes<-c(chr10_Testis_genes,a)
    chr10_Testis_overloci<-rbind(chr10_Testis_overloci,b)
  }
}

chr11_Testis_genes<-NULL
chr11_Testis_overloci<-NULL
if(length(same_Testis_chr11!=0)){
  for(i in 1:length(same_Testis_chr11)){
    a<-eqtl_Testis_chr11$phenotype_id[grep(same_Testis_chr11[i],eqtl_Testis_chr11$loci)]
    b<-eqtl_Testis_chr11[grep(same_Testis_chr11[i],eqtl_Testis_chr11$loci),]
    chr11_Testis_genes<-c(chr11_Testis_genes,a)
    chr11_Testis_overloci<-rbind(chr11_Testis_overloci,b)
  }
}

chr12_Testis_genes<-NULL
chr12_Testis_overloci<-NULL
if(length(same_Testis_chr12!=0)){
  for(i in 1:length(same_Testis_chr12)){
    a<-eqtl_Testis_chr12$phenotype_id[grep(same_Testis_chr12[i],eqtl_Testis_chr12$loci)]
    b<-eqtl_Testis_chr12[grep(same_Testis_chr12[i],eqtl_Testis_chr12$loci),]
    chr12_Testis_genes<-c(chr12_Testis_genes,a)
    chr12_Testis_overloci<-rbind(chr12_Testis_overloci,b)
  }
}

chr13_Testis_genes<-NULL
chr13_Testis_overloci<-NULL
if(length(same_Testis_chr13!=0)){
  for(i in 1:length(same_Testis_chr13)){
    a<-eqtl_Testis_chr13$phenotype_id[grep(same_Testis_chr13[i],eqtl_Testis_chr13$loci)]
    b<-eqtl_Testis_chr13[grep(same_Testis_chr13[i],eqtl_Testis_chr13$loci),]
    chr13_Testis_genes<-c(chr13_Testis_genes,a)
    chr13_Testis_overloci<-rbind(chr13_Testis_overloci,b)
  }
}

chr14_Testis_genes<-NULL
chr14_Testis_overloci<-NULL
if(length(same_Testis_chr14!=0)){
  for(i in 1:length(same_Testis_chr14)){
    a<-eqtl_Testis_chr14$phenotype_id[grep(same_Testis_chr14[i],eqtl_Testis_chr14$loci)]
    b<-eqtl_Testis_chr14[grep(same_Testis_chr14[i],eqtl_Testis_chr14$loci),]
    chr14_Testis_genes<-c(chr14_Testis_genes,a)
    chr14_Testis_overloci<-rbind(chr14_Testis_overloci,b)
  }
}

chr15_Testis_genes<-NULL
chr15_Testis_overloci<-NULL
if(length(same_Testis_chr15!=0)){
  for(i in 1:length(same_Testis_chr15)){
    a<-eqtl_Testis_chr15$phenotype_id[grep(same_Testis_chr15[i],eqtl_Testis_chr15$loci)]
    b<-eqtl_Testis_chr15[grep(same_Testis_chr15[i],eqtl_Testis_chr15$loci),]
    chr15_Testis_genes<-c(chr15_Testis_genes,a)
    chr15_Testis_overloci<-rbind(chr15_Testis_overloci,b)
  }
}

chr16_Testis_genes<-NULL
chr16_Testis_overloci<-NULL
if(length(same_Testis_chr16!=0)){
  for(i in 1:length(same_Testis_chr16)){
    a<-eqtl_Testis_chr16$phenotype_id[grep(same_Testis_chr16[i],eqtl_Testis_chr16$loci)]
    b<-eqtl_Testis_chr16[grep(same_Testis_chr16[i],eqtl_Testis_chr16$loci),]
    chr16_Testis_genes<-c(chr16_Testis_genes,a)
    chr16_Testis_overloci<-rbind(chr16_Testis_overloci,b)
  }
}

chr17_Testis_genes<-NULL
chr17_Testis_overloci<-NULL
if(length(same_Testis_chr17!=0)){
  for(i in 1:length(same_Testis_chr17)){
    a<-eqtl_Testis_chr17$phenotype_id[grep(same_Testis_chr17[i],eqtl_Testis_chr17$loci)]
    b<-eqtl_Testis_chr17[grep(same_Testis_chr17[i],eqtl_Testis_chr17$loci),]
    chr17_Testis_genes<-c(chr17_Testis_genes,a)
    chr17_Testis_overloci<-rbind(chr17_Testis_overloci,b)
  }
}

chr18_Testis_genes<-NULL
chr18_Testis_overloci<-NULL
if(length(same_Testis_chr18!=0)){
  for(i in 1:length(same_Testis_chr18)){
    a<-eqtl_Testis_chr18$phenotype_id[grep(same_Testis_chr18[i],eqtl_Testis_chr18$loci)]
    b<-eqtl_Testis_chr18[grep(same_Testis_chr18[i],eqtl_Testis_chr18$loci),]
    chr18_Testis_genes<-c(chr18_Testis_genes,a)
    chr18_Testis_overloci<-rbind(chr18_Testis_overloci,b)
  }
}
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
one2one_pig<-annotation$`Pig gene stable ID`[which(annotation$`Pig homology type`=="ortholog_one2one")]

chr1_Testis_overloci_one2one<-chr1_Testis_overloci[match(intersect(chr1_Testis_genes,one2one_pig),chr1_Testis_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Testis_chr1$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Testis_chr1$hum_chr<-tmp$V1
Testis_chr1$hum_loci<-tmp$V3
FM_Testis_hum$index<-paste0(FM_Testis_hum$chr,"_",FM_Testis_hum$variant_pos)
Testis_chr1$index<-paste0("chr",Testis_chr1$hum_chr,"_",Testis_chr1$hum_loci)
chr1_Testis_hum<-NULL
for(i in 1:length(same_Testis_chr1)){
  a<-Testis_chr1[match(same_Testis_chr1[i],Testis_chr1$pig_pos),]
  b<-FM_Testis_hum[match(intersect(a$index,FM_Testis_hum$index),FM_Testis_hum$index),]
  chr1_Testis_hum<-rbind(chr1_Testis_hum,b)
}
chr1_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr1_Testis_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr1_Testis_hum_one2one<-chr1_Testis_hum[match(intersect(chr1_pig2hum_one2one,chr1_Testis_hum$gene_id),chr1_Testis_hum$gene_id),]

chr2_Testis_overloci_one2one<-chr2_Testis_overloci[match(intersect(chr2_Testis_genes,one2one_pig),chr2_Testis_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Testis_chr2$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Testis_chr2$hum_chr<-tmp$V1
Testis_chr2$hum_loci<-tmp$V3
FM_Testis_hum$index<-paste0(FM_Testis_hum$chr,"_",FM_Testis_hum$variant_pos)
Testis_chr2$index<-paste0("chr",Testis_chr2$hum_chr,"_",Testis_chr2$hum_loci)
chr2_Testis_hum<-NULL
for(i in 1:length(same_Testis_chr2)){
  a<-Testis_chr2[match(same_Testis_chr2[i],Testis_chr2$pig_pos),]
  b<-FM_Testis_hum[match(intersect(a$index,FM_Testis_hum$index),FM_Testis_hum$index),]
  chr2_Testis_hum<-rbind(chr2_Testis_hum,b)
}
chr2_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr2_Testis_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr2_Testis_hum_one2one<-chr2_Testis_hum[match(intersect(chr2_pig2hum_one2one,chr2_Testis_hum$gene_id),chr2_Testis_hum$gene_id),]

chr3_Testis_overloci_one2one<-chr3_Testis_overloci[match(intersect(chr3_Testis_genes,one2one_pig),chr3_Testis_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Testis_chr3$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Testis_chr3$hum_chr<-tmp$V1
Testis_chr3$hum_loci<-tmp$V3
FM_Testis_hum$index<-paste0(FM_Testis_hum$chr,"_",FM_Testis_hum$variant_pos)
Testis_chr3$index<-paste0("chr",Testis_chr3$hum_chr,"_",Testis_chr3$hum_loci)
chr3_Testis_hum<-NULL
for(i in 1:length(same_Testis_chr3)){
  a<-Testis_chr3[match(same_Testis_chr3[i],Testis_chr3$pig_pos),]
  b<-FM_Testis_hum[match(intersect(a$index,FM_Testis_hum$index),FM_Testis_hum$index),]
  chr3_Testis_hum<-rbind(chr3_Testis_hum,b)
}
chr3_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr3_Testis_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr3_Testis_hum_one2one<-chr3_Testis_hum[match(intersect(chr3_pig2hum_one2one,chr3_Testis_hum$gene_id),chr3_Testis_hum$gene_id),]

chr4_Testis_overloci_one2one<-chr4_Testis_overloci[match(intersect(chr4_Testis_genes,one2one_pig),chr4_Testis_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Testis_chr4$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Testis_chr4$hum_chr<-tmp$V1
Testis_chr4$hum_loci<-tmp$V3
FM_Testis_hum$index<-paste0(FM_Testis_hum$chr,"_",FM_Testis_hum$variant_pos)
Testis_chr4$index<-paste0("chr",Testis_chr4$hum_chr,"_",Testis_chr4$hum_loci)
chr4_Testis_hum<-NULL
for(i in 1:length(same_Testis_chr4)){
  a<-Testis_chr4[match(same_Testis_chr4[i],Testis_chr4$pig_pos),]
  b<-FM_Testis_hum[match(intersect(a$index,FM_Testis_hum$index),FM_Testis_hum$index),]
  chr4_Testis_hum<-rbind(chr4_Testis_hum,b)
}
chr4_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr4_Testis_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr4_Testis_hum_one2one<-chr4_Testis_hum[match(intersect(chr4_pig2hum_one2one,chr4_Testis_hum$gene_id),chr4_Testis_hum$gene_id),]

chr5_Testis_overloci_one2one<-chr5_Testis_overloci[match(intersect(chr5_Testis_genes,one2one_pig),chr5_Testis_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Testis_chr5$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Testis_chr5$hum_chr<-tmp$V1
Testis_chr5$hum_loci<-tmp$V3
FM_Testis_hum$index<-paste0(FM_Testis_hum$chr,"_",FM_Testis_hum$variant_pos)
Testis_chr5$index<-paste0("chr",Testis_chr5$hum_chr,"_",Testis_chr5$hum_loci)
chr5_Testis_hum<-NULL
for(i in 1:length(same_Testis_chr5)){
  a<-Testis_chr5[match(same_Testis_chr5[i],Testis_chr5$pig_pos),]
  b<-FM_Testis_hum[match(intersect(a$index,FM_Testis_hum$index),FM_Testis_hum$index),]
  chr5_Testis_hum<-rbind(chr5_Testis_hum,b)
}
chr5_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr5_Testis_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr5_Testis_hum_one2one<-chr5_Testis_hum[match(intersect(chr5_pig2hum_one2one,chr5_Testis_hum$gene_id),chr5_Testis_hum$gene_id),]

chr6_Testis_overloci_one2one<-chr6_Testis_overloci[match(intersect(chr6_Testis_genes,one2one_pig),chr6_Testis_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Testis_chr6$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Testis_chr6$hum_chr<-tmp$V1
Testis_chr6$hum_loci<-tmp$V3
FM_Testis_hum$index<-paste0(FM_Testis_hum$chr,"_",FM_Testis_hum$variant_pos)
Testis_chr6$index<-paste0("chr",Testis_chr6$hum_chr,"_",Testis_chr6$hum_loci)
chr6_Testis_hum<-NULL
for(i in 1:length(same_Testis_chr6)){
  a<-Testis_chr6[match(same_Testis_chr6[i],Testis_chr6$pig_pos),]
  b<-FM_Testis_hum[match(intersect(a$index,FM_Testis_hum$index),FM_Testis_hum$index),]
  chr6_Testis_hum<-rbind(chr6_Testis_hum,b)
}
chr6_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr6_Testis_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr6_Testis_hum_one2one<-chr6_Testis_hum[match(intersect(chr6_pig2hum_one2one,chr6_Testis_hum$gene_id),chr6_Testis_hum$gene_id),]

chr7_Testis_overloci_one2one<-chr7_Testis_overloci[match(intersect(chr7_Testis_genes,one2one_pig),chr7_Testis_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Testis_chr7$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Testis_chr7$hum_chr<-tmp$V1
Testis_chr7$hum_loci<-tmp$V3
FM_Testis_hum$index<-paste0(FM_Testis_hum$chr,"_",FM_Testis_hum$variant_pos)
Testis_chr7$index<-paste0("chr",Testis_chr7$hum_chr,"_",Testis_chr7$hum_loci)
chr7_Testis_hum<-NULL
for(i in 1:length(same_Testis_chr7)){
  a<-Testis_chr7[match(same_Testis_chr7[i],Testis_chr7$pig_pos),]
  b<-FM_Testis_hum[match(intersect(a$index,FM_Testis_hum$index),FM_Testis_hum$index),]
  chr7_Testis_hum<-rbind(chr7_Testis_hum,b)
}
chr7_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr7_Testis_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr7_Testis_hum_one2one<-chr7_Testis_hum[match(intersect(chr7_pig2hum_one2one,chr7_Testis_hum$gene_id),chr7_Testis_hum$gene_id),]

chr8_Testis_overloci_one2one<-chr8_Testis_overloci[match(intersect(chr8_Testis_genes,one2one_pig),chr8_Testis_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Testis_chr8$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Testis_chr8$hum_chr<-tmp$V1
Testis_chr8$hum_loci<-tmp$V3
FM_Testis_hum$index<-paste0(FM_Testis_hum$chr,"_",FM_Testis_hum$variant_pos)
Testis_chr8$index<-paste0("chr",Testis_chr8$hum_chr,"_",Testis_chr8$hum_loci)
chr8_Testis_hum<-NULL
for(i in 1:length(same_Testis_chr8)){
  a<-Testis_chr8[match(same_Testis_chr8[i],Testis_chr8$pig_pos),]
  b<-FM_Testis_hum[match(intersect(a$index,FM_Testis_hum$index),FM_Testis_hum$index),]
  chr8_Testis_hum<-rbind(chr8_Testis_hum,b)
}
chr8_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr8_Testis_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr8_Testis_hum_one2one<-chr8_Testis_hum[match(intersect(chr8_pig2hum_one2one,chr8_Testis_hum$gene_id),chr8_Testis_hum$gene_id),]

chr9_Testis_overloci_one2one<-chr9_Testis_overloci[match(intersect(chr9_Testis_genes,one2one_pig),chr9_Testis_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Testis_chr9$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Testis_chr9$hum_chr<-tmp$V1
Testis_chr9$hum_loci<-tmp$V3
FM_Testis_hum$index<-paste0(FM_Testis_hum$chr,"_",FM_Testis_hum$variant_pos)
Testis_chr9$index<-paste0("chr",Testis_chr9$hum_chr,"_",Testis_chr9$hum_loci)
chr9_Testis_hum<-NULL
for(i in 1:length(same_Testis_chr9)){
  a<-Testis_chr9[match(same_Testis_chr9[i],Testis_chr9$pig_pos),]
  b<-FM_Testis_hum[match(intersect(a$index,FM_Testis_hum$index),FM_Testis_hum$index),]
  chr9_Testis_hum<-rbind(chr9_Testis_hum,b)
}
chr9_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr9_Testis_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr9_Testis_hum_one2one<-chr9_Testis_hum[match(intersect(chr9_pig2hum_one2one,chr9_Testis_hum$gene_id),chr9_Testis_hum$gene_id),]

chr10_Testis_overloci_one2one<-chr10_Testis_overloci[match(intersect(chr10_Testis_genes,one2one_pig),chr10_Testis_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Testis_chr10$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Testis_chr10$hum_chr<-tmp$V1
Testis_chr10$hum_loci<-tmp$V3
FM_Testis_hum$index<-paste0(FM_Testis_hum$chr,"_",FM_Testis_hum$variant_pos)
Testis_chr10$index<-paste0("chr",Testis_chr10$hum_chr,"_",Testis_chr10$hum_loci)
chr10_Testis_hum<-NULL
for(i in 1:length(same_Testis_chr10)){
  a<-Testis_chr10[match(same_Testis_chr10[i],Testis_chr10$pig_pos),]
  b<-FM_Testis_hum[match(intersect(a$index,FM_Testis_hum$index),FM_Testis_hum$index),]
  chr10_Testis_hum<-rbind(chr10_Testis_hum,b)
}
chr10_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr10_Testis_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr10_Testis_hum_one2one<-chr10_Testis_hum[match(intersect(chr10_pig2hum_one2one,chr10_Testis_hum$gene_id),chr10_Testis_hum$gene_id),]

chr11_Testis_overloci_one2one<-chr11_Testis_overloci[match(intersect(chr11_Testis_genes,one2one_pig),chr11_Testis_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Testis_chr11$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Testis_chr11$hum_chr<-tmp$V1
Testis_chr11$hum_loci<-tmp$V3
FM_Testis_hum$index<-paste0(FM_Testis_hum$chr,"_",FM_Testis_hum$variant_pos)
Testis_chr11$index<-paste0("chr",Testis_chr11$hum_chr,"_",Testis_chr11$hum_loci)
chr11_Testis_hum<-NULL
for(i in 1:length(same_Testis_chr11)){
  a<-Testis_chr11[match(same_Testis_chr11[i],Testis_chr11$pig_pos),]
  b<-FM_Testis_hum[match(intersect(a$index,FM_Testis_hum$index),FM_Testis_hum$index),]
  chr11_Testis_hum<-rbind(chr11_Testis_hum,b)
}
chr11_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr11_Testis_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr11_Testis_hum_one2one<-chr11_Testis_hum[match(intersect(chr11_pig2hum_one2one,chr11_Testis_hum$gene_id),chr11_Testis_hum$gene_id),]

chr12_Testis_overloci_one2one<-chr12_Testis_overloci[match(intersect(chr12_Testis_genes,one2one_pig),chr12_Testis_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Testis_chr12$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Testis_chr12$hum_chr<-tmp$V1
Testis_chr12$hum_loci<-tmp$V3
FM_Testis_hum$index<-paste0(FM_Testis_hum$chr,"_",FM_Testis_hum$variant_pos)
Testis_chr12$index<-paste0("chr",Testis_chr12$hum_chr,"_",Testis_chr12$hum_loci)
chr12_Testis_hum<-NULL
for(i in 1:length(same_Testis_chr12)){
  a<-Testis_chr12[match(same_Testis_chr12[i],Testis_chr12$pig_pos),]
  b<-FM_Testis_hum[match(intersect(a$index,FM_Testis_hum$index),FM_Testis_hum$index),]
  chr12_Testis_hum<-rbind(chr12_Testis_hum,b)
}
chr12_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr12_Testis_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr12_Testis_hum_one2one<-chr12_Testis_hum[match(intersect(chr12_pig2hum_one2one,chr12_Testis_hum$gene_id),chr12_Testis_hum$gene_id),]

chr13_Testis_overloci_one2one<-chr13_Testis_overloci[match(intersect(chr13_Testis_genes,one2one_pig),chr13_Testis_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Testis_chr13$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Testis_chr13$hum_chr<-tmp$V1
Testis_chr13$hum_loci<-tmp$V3
FM_Testis_hum$index<-paste0(FM_Testis_hum$chr,"_",FM_Testis_hum$variant_pos)
Testis_chr13$index<-paste0("chr",Testis_chr13$hum_chr,"_",Testis_chr13$hum_loci)
chr13_Testis_hum<-NULL
for(i in 1:length(same_Testis_chr13)){
  a<-Testis_chr13[match(same_Testis_chr13[i],Testis_chr13$pig_pos),]
  b<-FM_Testis_hum[match(intersect(a$index,FM_Testis_hum$index),FM_Testis_hum$index),]
  chr13_Testis_hum<-rbind(chr13_Testis_hum,b)
}
chr13_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr13_Testis_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr13_Testis_hum_one2one<-chr13_Testis_hum[match(intersect(chr13_pig2hum_one2one,chr13_Testis_hum$gene_id),chr13_Testis_hum$gene_id),]

chr14_Testis_overloci_one2one<-chr14_Testis_overloci[match(intersect(chr14_Testis_genes,one2one_pig),chr14_Testis_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Testis_chr14$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Testis_chr14$hum_chr<-tmp$V1
Testis_chr14$hum_loci<-tmp$V3
FM_Testis_hum$index<-paste0(FM_Testis_hum$chr,"_",FM_Testis_hum$variant_pos)
Testis_chr14$index<-paste0("chr",Testis_chr14$hum_chr,"_",Testis_chr14$hum_loci)
chr14_Testis_hum<-NULL
for(i in 1:length(same_Testis_chr14)){
  a<-Testis_chr14[match(same_Testis_chr14[i],Testis_chr14$pig_pos),]
  b<-FM_Testis_hum[match(intersect(a$index,FM_Testis_hum$index),FM_Testis_hum$index),]
  chr14_Testis_hum<-rbind(chr14_Testis_hum,b)
}
chr14_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr14_Testis_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr14_Testis_hum_one2one<-chr14_Testis_hum[match(intersect(chr14_pig2hum_one2one,chr14_Testis_hum$gene_id),chr14_Testis_hum$gene_id),]

chr15_Testis_overloci_one2one<-chr15_Testis_overloci[match(intersect(chr15_Testis_genes,one2one_pig),chr15_Testis_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Testis_chr15$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Testis_chr15$hum_chr<-tmp$V1
Testis_chr15$hum_loci<-tmp$V3
FM_Testis_hum$index<-paste0(FM_Testis_hum$chr,"_",FM_Testis_hum$variant_pos)
Testis_chr15$index<-paste0("chr",Testis_chr15$hum_chr,"_",Testis_chr15$hum_loci)
chr15_Testis_hum<-NULL
for(i in 1:length(same_Testis_chr15)){
  a<-Testis_chr15[match(same_Testis_chr15[i],Testis_chr15$pig_pos),]
  b<-FM_Testis_hum[match(intersect(a$index,FM_Testis_hum$index),FM_Testis_hum$index),]
  chr15_Testis_hum<-rbind(chr15_Testis_hum,b)
}
chr15_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr15_Testis_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr15_Testis_hum_one2one<-chr15_Testis_hum[match(intersect(chr15_pig2hum_one2one,chr15_Testis_hum$gene_id),chr15_Testis_hum$gene_id),]

chr16_Testis_overloci_one2one<-chr16_Testis_overloci[match(intersect(chr16_Testis_genes,one2one_pig),chr16_Testis_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Testis_chr16$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Testis_chr16$hum_chr<-tmp$V1
Testis_chr16$hum_loci<-tmp$V3
FM_Testis_hum$index<-paste0(FM_Testis_hum$chr,"_",FM_Testis_hum$variant_pos)
Testis_chr16$index<-paste0("chr",Testis_chr16$hum_chr,"_",Testis_chr16$hum_loci)
chr16_Testis_hum<-NULL
for(i in 1:length(same_Testis_chr16)){
  a<-Testis_chr16[match(same_Testis_chr16[i],Testis_chr16$pig_pos),]
  b<-FM_Testis_hum[match(intersect(a$index,FM_Testis_hum$index),FM_Testis_hum$index),]
  chr16_Testis_hum<-rbind(chr16_Testis_hum,b)
}
chr16_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr16_Testis_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr16_Testis_hum_one2one<-chr16_Testis_hum[match(intersect(chr16_pig2hum_one2one,chr16_Testis_hum$gene_id),chr16_Testis_hum$gene_id),]

chr17_Testis_overloci_one2one<-chr17_Testis_overloci[match(intersect(chr17_Testis_genes,one2one_pig),chr17_Testis_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Testis_chr17$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Testis_chr17$hum_chr<-tmp$V1
Testis_chr17$hum_loci<-tmp$V3
FM_Testis_hum$index<-paste0(FM_Testis_hum$chr,"_",FM_Testis_hum$variant_pos)
Testis_chr17$index<-paste0("chr",Testis_chr17$hum_chr,"_",Testis_chr17$hum_loci)
chr17_Testis_hum<-NULL
for(i in 1:length(same_Testis_chr17)){
  a<-Testis_chr17[match(same_Testis_chr17[i],Testis_chr17$pig_pos),]
  b<-FM_Testis_hum[match(intersect(a$index,FM_Testis_hum$index),FM_Testis_hum$index),]
  chr17_Testis_hum<-rbind(chr17_Testis_hum,b)
}
chr17_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr17_Testis_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr17_Testis_hum_one2one<-chr17_Testis_hum[match(intersect(chr17_pig2hum_one2one,chr17_Testis_hum$gene_id),chr17_Testis_hum$gene_id),]

chr18_Testis_overloci_one2one<-chr18_Testis_overloci[match(intersect(chr18_Testis_genes,one2one_pig),chr18_Testis_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Testis_chr18$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Testis_chr18$hum_chr<-tmp$V1
Testis_chr18$hum_loci<-tmp$V3
FM_Testis_hum$index<-paste0(FM_Testis_hum$chr,"_",FM_Testis_hum$variant_pos)
Testis_chr18$index<-paste0("chr",Testis_chr18$hum_chr,"_",Testis_chr18$hum_loci)
chr18_Testis_hum<-NULL
for(i in 1:length(same_Testis_chr18)){
  a<-Testis_chr18[match(same_Testis_chr18[i],Testis_chr18$pig_pos),]
  b<-FM_Testis_hum[match(intersect(a$index,FM_Testis_hum$index),FM_Testis_hum$index),]
  chr18_Testis_hum<-rbind(chr18_Testis_hum,b)
}
chr18_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr18_Testis_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr18_Testis_hum_one2one<-chr18_Testis_hum[match(intersect(chr18_pig2hum_one2one,chr18_Testis_hum$gene_id),chr18_Testis_hum$gene_id),]

Testis_one2one_SNP_hum<-rbind(chr1_Testis_hum_one2one,chr2_Testis_hum_one2one,chr3_Testis_hum_one2one,chr4_Testis_hum_one2one,chr5_Testis_hum_one2one,
                              chr6_Testis_hum_one2one,chr7_Testis_hum_one2one,chr8_Testis_hum_one2one,chr9_Testis_hum_one2one,chr10_Testis_hum_one2one,
                              chr11_Testis_hum_one2one,chr12_Testis_hum_one2one,chr13_Testis_hum_one2one,chr14_Testis_hum_one2one,chr15_Testis_hum_one2one,
                              chr16_Testis_hum_one2one,chr17_Testis_hum_one2one,chr18_Testis_hum_one2one)

chr1_Testis_pig_one2one<-chr1_Testis_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr1_Testis_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr1_Testis_overloci_one2one$phenotype_id),1:9]
chr2_Testis_pig_one2one<-chr2_Testis_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr2_Testis_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr2_Testis_overloci_one2one$phenotype_id),1:9]
chr3_Testis_pig_one2one<-chr3_Testis_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr3_Testis_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr3_Testis_overloci_one2one$phenotype_id),1:9]
chr4_Testis_pig_one2one<-chr4_Testis_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr4_Testis_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr4_Testis_overloci_one2one$phenotype_id),1:9]
chr5_Testis_pig_one2one<-chr5_Testis_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr5_Testis_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr5_Testis_overloci_one2one$phenotype_id),1:9]
chr6_Testis_pig_one2one<-chr6_Testis_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr6_Testis_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr6_Testis_overloci_one2one$phenotype_id),1:9]
chr7_Testis_pig_one2one<-chr7_Testis_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr7_Testis_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr7_Testis_overloci_one2one$phenotype_id),1:9]
chr8_Testis_pig_one2one<-chr8_Testis_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr8_Testis_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr8_Testis_overloci_one2one$phenotype_id),1:9]
chr9_Testis_pig_one2one<-chr9_Testis_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr9_Testis_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr9_Testis_overloci_one2one$phenotype_id),1:9]
chr10_Testis_pig_one2one<-chr10_Testis_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr10_Testis_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr10_Testis_overloci_one2one$phenotype_id),1:9]
chr11_Testis_pig_one2one<-chr11_Testis_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr11_Testis_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr11_Testis_overloci_one2one$phenotype_id),1:9]
chr12_Testis_pig_one2one<-chr12_Testis_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr12_Testis_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr12_Testis_overloci_one2one$phenotype_id),1:9]
chr13_Testis_pig_one2one<-chr13_Testis_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr13_Testis_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr13_Testis_overloci_one2one$phenotype_id),1:9]
chr14_Testis_pig_one2one<-chr14_Testis_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr14_Testis_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr14_Testis_overloci_one2one$phenotype_id),1:9]
chr15_Testis_pig_one2one<-chr15_Testis_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr15_Testis_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr15_Testis_overloci_one2one$phenotype_id),1:9]
chr16_Testis_pig_one2one<-chr16_Testis_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr16_Testis_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr16_Testis_overloci_one2one$phenotype_id),1:9]
chr17_Testis_pig_one2one<-chr17_Testis_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr17_Testis_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr17_Testis_overloci_one2one$phenotype_id),1:9]
chr18_Testis_pig_one2one<-chr18_Testis_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr18_Testis_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr18_Testis_overloci_one2one$phenotype_id),1:9]

Testis_one2one_SNP_pig<-rbind(chr1_Testis_pig_one2one,chr2_Testis_pig_one2one,chr3_Testis_pig_one2one,chr4_Testis_pig_one2one,chr5_Testis_pig_one2one,
                              chr6_Testis_pig_one2one,chr7_Testis_pig_one2one,chr8_Testis_pig_one2one,chr9_Testis_pig_one2one,chr10_Testis_pig_one2one,
                              chr11_Testis_pig_one2one,chr12_Testis_pig_one2one,chr13_Testis_pig_one2one,chr14_Testis_pig_one2one,chr15_Testis_pig_one2one,
                              chr16_Testis_pig_one2one,chr17_Testis_pig_one2one,chr18_Testis_pig_one2one)

Testis_SNP_sum<-array(NA,dim=c(nrow(Testis_one2one_SNP_hum),2))
colnames(Testis_SNP_sum)<-c("Human","Pig")
Testis_SNP_sum<-as.data.frame(Testis_SNP_sum)
Testis_SNP_sum$Human<-Testis_one2one_SNP_hum$slope / Testis_one2one_SNP_hum$slope_se
Testis_SNP_sum$Pig<-Testis_one2one_SNP_pig$slope / Testis_one2one_SNP_pig$slope_se
cor<-cor(abs(Testis_SNP_sum$Human),abs(Testis_SNP_sum$Pig))
p_val<-t.test(abs(Testis_SNP_sum$Human),abs(Testis_SNP_sum$Pig))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/Testis_SNP_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(Testis_SNP_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The number of eGenes expressed in common tissues",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

save(Testis_one2one_SNP_hum,Testis_one2one_SNP_pig,file="/Users/baizhonghao/Downloads/human-pig GTEx/Testis_SNP.Rdata")

#Uterus_SNP_overlaploci#
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

eqtl_Uterus_chr1<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Uterus/Uterus.cis_qtl_pairs.1.txt"))
eqtl_Uterus_chr2<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Uterus/Uterus.cis_qtl_pairs.2.txt"))
eqtl_Uterus_chr3<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Uterus/Uterus.cis_qtl_pairs.3.txt"))
eqtl_Uterus_chr4<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Uterus/Uterus.cis_qtl_pairs.4.txt"))
eqtl_Uterus_chr5<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Uterus/Uterus.cis_qtl_pairs.5.txt"))
eqtl_Uterus_chr6<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Uterus/Uterus.cis_qtl_pairs.6.txt"))
eqtl_Uterus_chr7<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Uterus/Uterus.cis_qtl_pairs.7.txt"))
eqtl_Uterus_chr8<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Uterus/Uterus.cis_qtl_pairs.8.txt"))
eqtl_Uterus_chr9<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Uterus/Uterus.cis_qtl_pairs.9.txt"))
eqtl_Uterus_chr10<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Uterus/Uterus.cis_qtl_pairs.10.txt"))
eqtl_Uterus_chr11<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Uterus/Uterus.cis_qtl_pairs.11.txt"))
eqtl_Uterus_chr12<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Uterus/Uterus.cis_qtl_pairs.12.txt"))
eqtl_Uterus_chr13<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Uterus/Uterus.cis_qtl_pairs.13.txt"))
eqtl_Uterus_chr14<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Uterus/Uterus.cis_qtl_pairs.14.txt"))
eqtl_Uterus_chr15<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Uterus/Uterus.cis_qtl_pairs.15.txt"))
eqtl_Uterus_chr16<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Uterus/Uterus.cis_qtl_pairs.16.txt"))
eqtl_Uterus_chr17<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Uterus/Uterus.cis_qtl_pairs.17.txt"))
eqtl_Uterus_chr18<-as.data.frame(fread("/Users/baizhonghao/Downloads/output_pcg/Uterus/Uterus.cis_qtl_pairs.18.txt"))

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus_chr1$variant_id,split="_"))))
eqtl_Uterus_chr1$chr<-tmp$V1
eqtl_Uterus_chr1$loci<-tmp$V2

eqtl_Uterus_chr1$index<-paste0(eqtl_Uterus_chr1$chr,"-",eqtl_Uterus_chr1$loci)
Uterus_chr1$index<-paste0(Uterus_chr1$chr,"-",Uterus_chr1$pig_pos)
same_Uterus_chr1<-intersect(Uterus_chr1$pig_pos,eqtl_Uterus_chr1$loci)


tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus_chr2$variant_id,split="_"))))
eqtl_Uterus_chr2$chr<-tmp$V1
eqtl_Uterus_chr2$loci<-tmp$V2

eqtl_Uterus_chr2$index<-paste0(eqtl_Uterus_chr2$chr,"-",eqtl_Uterus_chr2$loci)
Uterus_chr2$index<-paste0(Uterus_chr2$chr,"-",Uterus_chr2$pig_pos)
same_Uterus_chr2<-intersect(Uterus_chr2$pig_pos,eqtl_Uterus_chr2$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus_chr3$variant_id,split="_"))))
eqtl_Uterus_chr3$chr<-tmp$V1
eqtl_Uterus_chr3$loci<-tmp$V2

eqtl_Uterus_chr3$index<-paste0(eqtl_Uterus_chr3$chr,"-",eqtl_Uterus_chr3$loci)
Uterus_chr3$index<-paste0(Uterus_chr3$chr,"-",Uterus_chr3$pig_pos)
same_Uterus_chr3<-intersect(Uterus_chr3$pig_pos,eqtl_Uterus_chr3$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus_chr4$variant_id,split="_"))))
eqtl_Uterus_chr4$chr<-tmp$V1
eqtl_Uterus_chr4$loci<-tmp$V2

eqtl_Uterus_chr4$index<-paste0(eqtl_Uterus_chr4$chr,"-",eqtl_Uterus_chr4$loci)
Uterus_chr4$index<-paste0(Uterus_chr4$chr,"-",Uterus_chr4$pig_pos)
same_Uterus_chr4<-intersect(Uterus_chr4$pig_pos,eqtl_Uterus_chr4$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus_chr5$variant_id,split="_"))))
eqtl_Uterus_chr5$loci<-tmp$V2
same_Uterus_chr5<-intersect(Uterus_chr5$pig_pos,eqtl_Uterus_chr5$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus_chr6$variant_id,split="_"))))
eqtl_Uterus_chr6$loci<-tmp$V2
same_Uterus_chr6<-intersect(Uterus_chr6$pig_pos,eqtl_Uterus_chr6$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus_chr7$variant_id,split="_"))))
eqtl_Uterus_chr7$loci<-tmp$V2
same_Uterus_chr7<-intersect(Uterus_chr7$pig_pos,eqtl_Uterus_chr7$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus_chr8$variant_id,split="_"))))
eqtl_Uterus_chr8$loci<-tmp$V2
same_Uterus_chr8<-intersect(Uterus_chr8$pig_pos,eqtl_Uterus_chr8$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus_chr9$variant_id,split="_"))))
eqtl_Uterus_chr9$loci<-tmp$V2
same_Uterus_chr9<-intersect(Uterus_chr9$pig_pos,eqtl_Uterus_chr9$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus_chr10$variant_id,split="_"))))
eqtl_Uterus_chr10$loci<-tmp$V2
same_Uterus_chr10<-intersect(Uterus_chr10$pig_pos,eqtl_Uterus_chr10$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus_chr11$variant_id,split="_"))))
eqtl_Uterus_chr11$loci<-tmp$V2
same_Uterus_chr11<-intersect(Uterus_chr11$pig_pos,eqtl_Uterus_chr11$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus_chr12$variant_id,split="_"))))
eqtl_Uterus_chr12$loci<-tmp$V2
same_Uterus_chr12<-intersect(Uterus_chr12$pig_pos,eqtl_Uterus_chr12$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus_chr13$variant_id,split="_"))))
eqtl_Uterus_chr13$loci<-tmp$V2
same_Uterus_chr13<-intersect(Uterus_chr13$pig_pos,eqtl_Uterus_chr13$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus_chr14$variant_id,split="_"))))
eqtl_Uterus_chr14$loci<-tmp$V2
same_Uterus_chr14<-intersect(Uterus_chr14$pig_pos,eqtl_Uterus_chr14$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus_chr15$variant_id,split="_"))))
eqtl_Uterus_chr15$loci<-tmp$V2
same_Uterus_chr15<-intersect(Uterus_chr15$pig_pos,eqtl_Uterus_chr15$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus_chr16$variant_id,split="_"))))
eqtl_Uterus_chr16$loci<-tmp$V2
same_Uterus_chr16<-intersect(Uterus_chr16$pig_pos,eqtl_Uterus_chr16$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus_chr17$variant_id,split="_"))))
eqtl_Uterus_chr17$loci<-tmp$V2
same_Uterus_chr17<-intersect(Uterus_chr17$pig_pos,eqtl_Uterus_chr17$loci)

tmp<-as.data.frame(t(as.data.frame(strsplit(eqtl_Uterus_chr18$variant_id,split="_"))))
eqtl_Uterus_chr18$loci<-tmp$V2
same_Uterus_chr18<-intersect(Uterus_chr18$pig_pos,eqtl_Uterus_chr18$loci)

chr1_Uterus_genes<-NULL
chr1_Uterus_overloci<-NULL
if(length(same_Uterus_chr1!=0)){
  for(i in 1:length(same_Uterus_chr1)){
    a<-eqtl_Uterus_chr1$phenotype_id[grep(same_Uterus_chr1[i],eqtl_Uterus_chr1$loci)]
    b<-eqtl_Uterus_chr1[grep(same_Uterus_chr1[i],eqtl_Uterus_chr1$loci),]
    chr1_Uterus_genes<-c(chr1_Uterus_genes,a)
    chr1_Uterus_overloci<-rbind(chr1_Uterus_overloci,b)
  }
}

chr2_Uterus_genes<-NULL
chr2_Uterus_overloci<-NULL
if(length(same_Uterus_chr2!=0)){
  for(i in 1:length(same_Uterus_chr2)){
    a<-eqtl_Uterus_chr2$phenotype_id[grep(same_Uterus_chr2[i],eqtl_Uterus_chr2$loci)]
    b<-eqtl_Uterus_chr2[grep(same_Uterus_chr2[i],eqtl_Uterus_chr2$loci),]
    chr2_Uterus_genes<-c(chr2_Uterus_genes,a)
    chr2_Uterus_overloci<-rbind(chr2_Uterus_overloci,b)
  }
}

chr3_Uterus_genes<-NULL
chr3_Uterus_overloci<-NULL
if(length(same_Uterus_chr3!=0)){
  for(i in 1:length(same_Uterus_chr3)){
    a<-eqtl_Uterus_chr3$phenotype_id[grep(same_Uterus_chr3[i],eqtl_Uterus_chr3$loci)]
    b<-eqtl_Uterus_chr3[grep(same_Uterus_chr3[i],eqtl_Uterus_chr3$loci),]
    chr3_Uterus_genes<-c(chr3_Uterus_genes,a)
    chr3_Uterus_overloci<-rbind(chr3_Uterus_overloci,b)
  }
}

chr4_Uterus_genes<-NULL
chr4_Uterus_overloci<-NULL
if(length(same_Uterus_chr4!=0)){
  for(i in 1:length(same_Uterus_chr4)){
    a<-eqtl_Uterus_chr4$phenotype_id[grep(same_Uterus_chr4[i],eqtl_Uterus_chr4$loci)]
    b<-eqtl_Uterus_chr4[grep(same_Uterus_chr4[i],eqtl_Uterus_chr4$loci),]
    chr4_Uterus_genes<-c(chr4_Uterus_genes,a)
    chr4_Uterus_overloci<-rbind(chr4_Uterus_overloci,b)
  }
}

chr5_Uterus_genes<-NULL
chr5_Uterus_overloci<-NULL
if(length(same_Uterus_chr5!=0)){
  for(i in 1:length(same_Uterus_chr5)){
    a<-eqtl_Uterus_chr5$phenotype_id[grep(same_Uterus_chr5[i],eqtl_Uterus_chr5$loci)]
    b<-eqtl_Uterus_chr5[grep(same_Uterus_chr5[i],eqtl_Uterus_chr5$loci),]
    chr5_Uterus_genes<-c(chr5_Uterus_genes,a)
    chr5_Uterus_overloci<-rbind(chr5_Uterus_overloci,b)
  }
}

chr6_Uterus_genes<-NULL
chr6_Uterus_overloci<-NULL
if(length(same_Uterus_chr6!=0)){
  for(i in 1:length(same_Uterus_chr6)){
    a<-eqtl_Uterus_chr6$phenotype_id[grep(same_Uterus_chr6[i],eqtl_Uterus_chr6$loci)]
    b<-eqtl_Uterus_chr6[grep(same_Uterus_chr6[i],eqtl_Uterus_chr6$loci),]
    chr6_Uterus_genes<-c(chr6_Uterus_genes,a)
    chr6_Uterus_overloci<-rbind(chr6_Uterus_overloci,b)
  }
}

chr7_Uterus_genes<-NULL
chr7_Uterus_overloci<-NULL
if(length(same_Uterus_chr7!=0)){
  for(i in 1:length(same_Uterus_chr7)){
    a<-eqtl_Uterus_chr7$phenotype_id[grep(same_Uterus_chr7[i],eqtl_Uterus_chr7$loci)]
    b<-eqtl_Uterus_chr7[grep(same_Uterus_chr7[i],eqtl_Uterus_chr7$loci),]
    chr7_Uterus_genes<-c(chr7_Uterus_genes,a)
    chr7_Uterus_overloci<-rbind(chr7_Uterus_overloci,b)
  }
}

chr8_Uterus_genes<-NULL
chr8_Uterus_overloci<-NULL
if(length(same_Uterus_chr8!=0)){
  for(i in 1:length(same_Uterus_chr8)){
    a<-eqtl_Uterus_chr8$phenotype_id[grep(same_Uterus_chr8[i],eqtl_Uterus_chr8$loci)]
    b<-eqtl_Uterus_chr8[grep(same_Uterus_chr8[i],eqtl_Uterus_chr8$loci),]
    chr8_Uterus_genes<-c(chr8_Uterus_genes,a)
    chr8_Uterus_overloci<-rbind(chr8_Uterus_overloci,b)
  }
}

chr9_Uterus_genes<-NULL
chr9_Uterus_overloci<-NULL
if(length(same_Uterus_chr9!=0)){
  for(i in 1:length(same_Uterus_chr9)){
    a<-eqtl_Uterus_chr9$phenotype_id[grep(same_Uterus_chr9[i],eqtl_Uterus_chr9$loci)]
    b<-eqtl_Uterus_chr9[grep(same_Uterus_chr9[i],eqtl_Uterus_chr9$loci),]
    chr9_Uterus_genes<-c(chr9_Uterus_genes,a)
    chr9_Uterus_overloci<-rbind(chr9_Uterus_overloci,b)
  }
}

chr10_Uterus_genes<-NULL
chr10_Uterus_overloci<-NULL
if(length(same_Uterus_chr10!=0)){
  for(i in 1:length(same_Uterus_chr10)){
    a<-eqtl_Uterus_chr10$phenotype_id[grep(same_Uterus_chr10[i],eqtl_Uterus_chr10$loci)]
    b<-eqtl_Uterus_chr10[grep(same_Uterus_chr10[i],eqtl_Uterus_chr10$loci),]
    chr10_Uterus_genes<-c(chr10_Uterus_genes,a)
    chr10_Uterus_overloci<-rbind(chr10_Uterus_overloci,b)
  }
}

chr11_Uterus_genes<-NULL
chr11_Uterus_overloci<-NULL
if(length(same_Uterus_chr11!=0)){
  for(i in 1:length(same_Uterus_chr11)){
    a<-eqtl_Uterus_chr11$phenotype_id[grep(same_Uterus_chr11[i],eqtl_Uterus_chr11$loci)]
    b<-eqtl_Uterus_chr11[grep(same_Uterus_chr11[i],eqtl_Uterus_chr11$loci),]
    chr11_Uterus_genes<-c(chr11_Uterus_genes,a)
    chr11_Uterus_overloci<-rbind(chr11_Uterus_overloci,b)
  }
}

chr12_Uterus_genes<-NULL
chr12_Uterus_overloci<-NULL
if(length(same_Uterus_chr12!=0)){
  for(i in 1:length(same_Uterus_chr12)){
    a<-eqtl_Uterus_chr12$phenotype_id[grep(same_Uterus_chr12[i],eqtl_Uterus_chr12$loci)]
    b<-eqtl_Uterus_chr12[grep(same_Uterus_chr12[i],eqtl_Uterus_chr12$loci),]
    chr12_Uterus_genes<-c(chr12_Uterus_genes,a)
    chr12_Uterus_overloci<-rbind(chr12_Uterus_overloci,b)
  }
}

chr13_Uterus_genes<-NULL
chr13_Uterus_overloci<-NULL
if(length(same_Uterus_chr13!=0)){
  for(i in 1:length(same_Uterus_chr13)){
    a<-eqtl_Uterus_chr13$phenotype_id[grep(same_Uterus_chr13[i],eqtl_Uterus_chr13$loci)]
    b<-eqtl_Uterus_chr13[grep(same_Uterus_chr13[i],eqtl_Uterus_chr13$loci),]
    chr13_Uterus_genes<-c(chr13_Uterus_genes,a)
    chr13_Uterus_overloci<-rbind(chr13_Uterus_overloci,b)
  }
}

chr14_Uterus_genes<-NULL
chr14_Uterus_overloci<-NULL
if(length(same_Uterus_chr14!=0)){
  for(i in 1:length(same_Uterus_chr14)){
    a<-eqtl_Uterus_chr14$phenotype_id[grep(same_Uterus_chr14[i],eqtl_Uterus_chr14$loci)]
    b<-eqtl_Uterus_chr14[grep(same_Uterus_chr14[i],eqtl_Uterus_chr14$loci),]
    chr14_Uterus_genes<-c(chr14_Uterus_genes,a)
    chr14_Uterus_overloci<-rbind(chr14_Uterus_overloci,b)
  }
}

chr15_Uterus_genes<-NULL
chr15_Uterus_overloci<-NULL
if(length(same_Uterus_chr15!=0)){
  for(i in 1:length(same_Uterus_chr15)){
    a<-eqtl_Uterus_chr15$phenotype_id[grep(same_Uterus_chr15[i],eqtl_Uterus_chr15$loci)]
    b<-eqtl_Uterus_chr15[grep(same_Uterus_chr15[i],eqtl_Uterus_chr15$loci),]
    chr15_Uterus_genes<-c(chr15_Uterus_genes,a)
    chr15_Uterus_overloci<-rbind(chr15_Uterus_overloci,b)
  }
}

chr16_Uterus_genes<-NULL
chr16_Uterus_overloci<-NULL
if(length(same_Uterus_chr16!=0)){
  for(i in 1:length(same_Uterus_chr16)){
    a<-eqtl_Uterus_chr16$phenotype_id[grep(same_Uterus_chr16[i],eqtl_Uterus_chr16$loci)]
    b<-eqtl_Uterus_chr16[grep(same_Uterus_chr16[i],eqtl_Uterus_chr16$loci),]
    chr16_Uterus_genes<-c(chr16_Uterus_genes,a)
    chr16_Uterus_overloci<-rbind(chr16_Uterus_overloci,b)
  }
}

chr17_Uterus_genes<-NULL
chr17_Uterus_overloci<-NULL
if(length(same_Uterus_chr17!=0)){
  for(i in 1:length(same_Uterus_chr17)){
    a<-eqtl_Uterus_chr17$phenotype_id[grep(same_Uterus_chr17[i],eqtl_Uterus_chr17$loci)]
    b<-eqtl_Uterus_chr17[grep(same_Uterus_chr17[i],eqtl_Uterus_chr17$loci),]
    chr17_Uterus_genes<-c(chr17_Uterus_genes,a)
    chr17_Uterus_overloci<-rbind(chr17_Uterus_overloci,b)
  }
}

chr18_Uterus_genes<-NULL
chr18_Uterus_overloci<-NULL
if(length(same_Uterus_chr18!=0)){
  for(i in 1:length(same_Uterus_chr18)){
    a<-eqtl_Uterus_chr18$phenotype_id[grep(same_Uterus_chr18[i],eqtl_Uterus_chr18$loci)]
    b<-eqtl_Uterus_chr18[grep(same_Uterus_chr18[i],eqtl_Uterus_chr18$loci),]
    chr18_Uterus_genes<-c(chr18_Uterus_genes,a)
    chr18_Uterus_overloci<-rbind(chr18_Uterus_overloci,b)
  }
}
annotation<-as.data.frame(fread("/Users/baizhonghao/Downloads/human-pig GTEx/annotation.txt"))
one2one_pig<-annotation$`Pig gene stable ID`[which(annotation$`Pig homology type`=="ortholog_one2one")]

chr1_Uterus_overloci_one2one<-chr1_Uterus_overloci[match(intersect(chr1_Uterus_genes,one2one_pig),chr1_Uterus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Uterus_chr1$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Uterus_chr1$hum_chr<-tmp$V1
Uterus_chr1$hum_loci<-tmp$V3
FM_Uterus_hum$index<-paste0(FM_Uterus_hum$chr,"_",FM_Uterus_hum$variant_pos)
Uterus_chr1$index<-paste0("chr",Uterus_chr1$hum_chr,"_",Uterus_chr1$hum_loci)
chr1_Uterus_hum<-NULL
for(i in 1:length(same_Uterus_chr1)){
  a<-Uterus_chr1[match(same_Uterus_chr1[i],Uterus_chr1$pig_pos),]
  b<-FM_Uterus_hum[match(intersect(a$index,FM_Uterus_hum$index),FM_Uterus_hum$index),]
  chr1_Uterus_hum<-rbind(chr1_Uterus_hum,b)
}
chr1_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr1_Uterus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr1_Uterus_hum_one2one<-chr1_Uterus_hum[match(intersect(chr1_pig2hum_one2one,chr1_Uterus_hum$gene_id),chr1_Uterus_hum$gene_id),]

chr2_Uterus_overloci_one2one<-chr2_Uterus_overloci[match(intersect(chr2_Uterus_genes,one2one_pig),chr2_Uterus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Uterus_chr2$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Uterus_chr2$hum_chr<-tmp$V1
Uterus_chr2$hum_loci<-tmp$V3
FM_Uterus_hum$index<-paste0(FM_Uterus_hum$chr,"_",FM_Uterus_hum$variant_pos)
Uterus_chr2$index<-paste0("chr",Uterus_chr2$hum_chr,"_",Uterus_chr2$hum_loci)
chr2_Uterus_hum<-NULL
for(i in 1:length(same_Uterus_chr2)){
  a<-Uterus_chr2[match(same_Uterus_chr2[i],Uterus_chr2$pig_pos),]
  b<-FM_Uterus_hum[match(intersect(a$index,FM_Uterus_hum$index),FM_Uterus_hum$index),]
  chr2_Uterus_hum<-rbind(chr2_Uterus_hum,b)
}
chr2_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr2_Uterus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr2_Uterus_hum_one2one<-chr2_Uterus_hum[match(intersect(chr2_pig2hum_one2one,chr2_Uterus_hum$gene_id),chr2_Uterus_hum$gene_id),]

chr3_Uterus_overloci_one2one<-chr3_Uterus_overloci[match(intersect(chr3_Uterus_genes,one2one_pig),chr3_Uterus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Uterus_chr3$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Uterus_chr3$hum_chr<-tmp$V1
Uterus_chr3$hum_loci<-tmp$V3
FM_Uterus_hum$index<-paste0(FM_Uterus_hum$chr,"_",FM_Uterus_hum$variant_pos)
Uterus_chr3$index<-paste0("chr",Uterus_chr3$hum_chr,"_",Uterus_chr3$hum_loci)
chr3_Uterus_hum<-NULL
for(i in 1:length(same_Uterus_chr3)){
  a<-Uterus_chr3[match(same_Uterus_chr3[i],Uterus_chr3$pig_pos),]
  b<-FM_Uterus_hum[match(intersect(a$index,FM_Uterus_hum$index),FM_Uterus_hum$index),]
  chr3_Uterus_hum<-rbind(chr3_Uterus_hum,b)
}
chr3_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr3_Uterus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr3_Uterus_hum_one2one<-chr3_Uterus_hum[match(intersect(chr3_pig2hum_one2one,chr3_Uterus_hum$gene_id),chr3_Uterus_hum$gene_id),]

chr4_Uterus_overloci_one2one<-chr4_Uterus_overloci[match(intersect(chr4_Uterus_genes,one2one_pig),chr4_Uterus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Uterus_chr4$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Uterus_chr4$hum_chr<-tmp$V1
Uterus_chr4$hum_loci<-tmp$V3
FM_Uterus_hum$index<-paste0(FM_Uterus_hum$chr,"_",FM_Uterus_hum$variant_pos)
Uterus_chr4$index<-paste0("chr",Uterus_chr4$hum_chr,"_",Uterus_chr4$hum_loci)
chr4_Uterus_hum<-NULL
for(i in 1:length(same_Uterus_chr4)){
  a<-Uterus_chr4[match(same_Uterus_chr4[i],Uterus_chr4$pig_pos),]
  b<-FM_Uterus_hum[match(intersect(a$index,FM_Uterus_hum$index),FM_Uterus_hum$index),]
  chr4_Uterus_hum<-rbind(chr4_Uterus_hum,b)
}
chr4_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr4_Uterus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr4_Uterus_hum_one2one<-chr4_Uterus_hum[match(intersect(chr4_pig2hum_one2one,chr4_Uterus_hum$gene_id),chr4_Uterus_hum$gene_id),]

chr5_Uterus_overloci_one2one<-chr5_Uterus_overloci[match(intersect(chr5_Uterus_genes,one2one_pig),chr5_Uterus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Uterus_chr5$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Uterus_chr5$hum_chr<-tmp$V1
Uterus_chr5$hum_loci<-tmp$V3
FM_Uterus_hum$index<-paste0(FM_Uterus_hum$chr,"_",FM_Uterus_hum$variant_pos)
Uterus_chr5$index<-paste0("chr",Uterus_chr5$hum_chr,"_",Uterus_chr5$hum_loci)
chr5_Uterus_hum<-NULL
for(i in 1:length(same_Uterus_chr5)){
  a<-Uterus_chr5[match(same_Uterus_chr5[i],Uterus_chr5$pig_pos),]
  b<-FM_Uterus_hum[match(intersect(a$index,FM_Uterus_hum$index),FM_Uterus_hum$index),]
  chr5_Uterus_hum<-rbind(chr5_Uterus_hum,b)
}
chr5_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr5_Uterus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr5_Uterus_hum_one2one<-chr5_Uterus_hum[match(intersect(chr5_pig2hum_one2one,chr5_Uterus_hum$gene_id),chr5_Uterus_hum$gene_id),]

chr6_Uterus_overloci_one2one<-chr6_Uterus_overloci[match(intersect(chr6_Uterus_genes,one2one_pig),chr6_Uterus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Uterus_chr6$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Uterus_chr6$hum_chr<-tmp$V1
Uterus_chr6$hum_loci<-tmp$V3
FM_Uterus_hum$index<-paste0(FM_Uterus_hum$chr,"_",FM_Uterus_hum$variant_pos)
Uterus_chr6$index<-paste0("chr",Uterus_chr6$hum_chr,"_",Uterus_chr6$hum_loci)
chr6_Uterus_hum<-NULL
for(i in 1:length(same_Uterus_chr6)){
  a<-Uterus_chr6[match(same_Uterus_chr6[i],Uterus_chr6$pig_pos),]
  b<-FM_Uterus_hum[match(intersect(a$index,FM_Uterus_hum$index),FM_Uterus_hum$index),]
  chr6_Uterus_hum<-rbind(chr6_Uterus_hum,b)
}
chr6_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr6_Uterus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr6_Uterus_hum_one2one<-chr6_Uterus_hum[match(intersect(chr6_pig2hum_one2one,chr6_Uterus_hum$gene_id),chr6_Uterus_hum$gene_id),]

chr7_Uterus_overloci_one2one<-chr7_Uterus_overloci[match(intersect(chr7_Uterus_genes,one2one_pig),chr7_Uterus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Uterus_chr7$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Uterus_chr7$hum_chr<-tmp$V1
Uterus_chr7$hum_loci<-tmp$V3
FM_Uterus_hum$index<-paste0(FM_Uterus_hum$chr,"_",FM_Uterus_hum$variant_pos)
Uterus_chr7$index<-paste0("chr",Uterus_chr7$hum_chr,"_",Uterus_chr7$hum_loci)
chr7_Uterus_hum<-NULL
for(i in 1:length(same_Uterus_chr7)){
  a<-Uterus_chr7[match(same_Uterus_chr7[i],Uterus_chr7$pig_pos),]
  b<-FM_Uterus_hum[match(intersect(a$index,FM_Uterus_hum$index),FM_Uterus_hum$index),]
  chr7_Uterus_hum<-rbind(chr7_Uterus_hum,b)
}
chr7_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr7_Uterus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr7_Uterus_hum_one2one<-chr7_Uterus_hum[match(intersect(chr7_pig2hum_one2one,chr7_Uterus_hum$gene_id),chr7_Uterus_hum$gene_id),]

chr8_Uterus_overloci_one2one<-chr8_Uterus_overloci[match(intersect(chr8_Uterus_genes,one2one_pig),chr8_Uterus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Uterus_chr8$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Uterus_chr8$hum_chr<-tmp$V1
Uterus_chr8$hum_loci<-tmp$V3
FM_Uterus_hum$index<-paste0(FM_Uterus_hum$chr,"_",FM_Uterus_hum$variant_pos)
Uterus_chr8$index<-paste0("chr",Uterus_chr8$hum_chr,"_",Uterus_chr8$hum_loci)
chr8_Uterus_hum<-NULL
for(i in 1:length(same_Uterus_chr8)){
  a<-Uterus_chr8[match(same_Uterus_chr8[i],Uterus_chr8$pig_pos),]
  b<-FM_Uterus_hum[match(intersect(a$index,FM_Uterus_hum$index),FM_Uterus_hum$index),]
  chr8_Uterus_hum<-rbind(chr8_Uterus_hum,b)
}
chr8_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr8_Uterus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr8_Uterus_hum_one2one<-chr8_Uterus_hum[match(intersect(chr8_pig2hum_one2one,chr8_Uterus_hum$gene_id),chr8_Uterus_hum$gene_id),]

chr9_Uterus_overloci_one2one<-chr9_Uterus_overloci[match(intersect(chr9_Uterus_genes,one2one_pig),chr9_Uterus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Uterus_chr9$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Uterus_chr9$hum_chr<-tmp$V1
Uterus_chr9$hum_loci<-tmp$V3
FM_Uterus_hum$index<-paste0(FM_Uterus_hum$chr,"_",FM_Uterus_hum$variant_pos)
Uterus_chr9$index<-paste0("chr",Uterus_chr9$hum_chr,"_",Uterus_chr9$hum_loci)
chr9_Uterus_hum<-NULL
for(i in 1:length(same_Uterus_chr9)){
  a<-Uterus_chr9[match(same_Uterus_chr9[i],Uterus_chr9$pig_pos),]
  b<-FM_Uterus_hum[match(intersect(a$index,FM_Uterus_hum$index),FM_Uterus_hum$index),]
  chr9_Uterus_hum<-rbind(chr9_Uterus_hum,b)
}
chr9_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr9_Uterus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr9_Uterus_hum_one2one<-chr9_Uterus_hum[match(intersect(chr9_pig2hum_one2one,chr9_Uterus_hum$gene_id),chr9_Uterus_hum$gene_id),]

chr10_Uterus_overloci_one2one<-chr10_Uterus_overloci[match(intersect(chr10_Uterus_genes,one2one_pig),chr10_Uterus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Uterus_chr10$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Uterus_chr10$hum_chr<-tmp$V1
Uterus_chr10$hum_loci<-tmp$V3
FM_Uterus_hum$index<-paste0(FM_Uterus_hum$chr,"_",FM_Uterus_hum$variant_pos)
Uterus_chr10$index<-paste0("chr",Uterus_chr10$hum_chr,"_",Uterus_chr10$hum_loci)
chr10_Uterus_hum<-NULL
for(i in 1:length(same_Uterus_chr10)){
  a<-Uterus_chr10[match(same_Uterus_chr10[i],Uterus_chr10$pig_pos),]
  b<-FM_Uterus_hum[match(intersect(a$index,FM_Uterus_hum$index),FM_Uterus_hum$index),]
  chr10_Uterus_hum<-rbind(chr10_Uterus_hum,b)
}
chr10_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr10_Uterus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr10_Uterus_hum_one2one<-chr10_Uterus_hum[match(intersect(chr10_pig2hum_one2one,chr10_Uterus_hum$gene_id),chr10_Uterus_hum$gene_id),]

chr11_Uterus_overloci_one2one<-chr11_Uterus_overloci[match(intersect(chr11_Uterus_genes,one2one_pig),chr11_Uterus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Uterus_chr11$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Uterus_chr11$hum_chr<-tmp$V1
Uterus_chr11$hum_loci<-tmp$V3
FM_Uterus_hum$index<-paste0(FM_Uterus_hum$chr,"_",FM_Uterus_hum$variant_pos)
Uterus_chr11$index<-paste0("chr",Uterus_chr11$hum_chr,"_",Uterus_chr11$hum_loci)
chr11_Uterus_hum<-NULL
for(i in 1:length(same_Uterus_chr11)){
  a<-Uterus_chr11[match(same_Uterus_chr11[i],Uterus_chr11$pig_pos),]
  b<-FM_Uterus_hum[match(intersect(a$index,FM_Uterus_hum$index),FM_Uterus_hum$index),]
  chr11_Uterus_hum<-rbind(chr11_Uterus_hum,b)
}
chr11_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr11_Uterus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr11_Uterus_hum_one2one<-chr11_Uterus_hum[match(intersect(chr11_pig2hum_one2one,chr11_Uterus_hum$gene_id),chr11_Uterus_hum$gene_id),]

chr12_Uterus_overloci_one2one<-chr12_Uterus_overloci[match(intersect(chr12_Uterus_genes,one2one_pig),chr12_Uterus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Uterus_chr12$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Uterus_chr12$hum_chr<-tmp$V1
Uterus_chr12$hum_loci<-tmp$V3
FM_Uterus_hum$index<-paste0(FM_Uterus_hum$chr,"_",FM_Uterus_hum$variant_pos)
Uterus_chr12$index<-paste0("chr",Uterus_chr12$hum_chr,"_",Uterus_chr12$hum_loci)
chr12_Uterus_hum<-NULL
for(i in 1:length(same_Uterus_chr12)){
  a<-Uterus_chr12[match(same_Uterus_chr12[i],Uterus_chr12$pig_pos),]
  b<-FM_Uterus_hum[match(intersect(a$index,FM_Uterus_hum$index),FM_Uterus_hum$index),]
  chr12_Uterus_hum<-rbind(chr12_Uterus_hum,b)
}
chr12_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr12_Uterus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr12_Uterus_hum_one2one<-chr12_Uterus_hum[match(intersect(chr12_pig2hum_one2one,chr12_Uterus_hum$gene_id),chr12_Uterus_hum$gene_id),]

chr13_Uterus_overloci_one2one<-chr13_Uterus_overloci[match(intersect(chr13_Uterus_genes,one2one_pig),chr13_Uterus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Uterus_chr13$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Uterus_chr13$hum_chr<-tmp$V1
Uterus_chr13$hum_loci<-tmp$V3
FM_Uterus_hum$index<-paste0(FM_Uterus_hum$chr,"_",FM_Uterus_hum$variant_pos)
Uterus_chr13$index<-paste0("chr",Uterus_chr13$hum_chr,"_",Uterus_chr13$hum_loci)
chr13_Uterus_hum<-NULL
for(i in 1:length(same_Uterus_chr13)){
  a<-Uterus_chr13[match(same_Uterus_chr13[i],Uterus_chr13$pig_pos),]
  b<-FM_Uterus_hum[match(intersect(a$index,FM_Uterus_hum$index),FM_Uterus_hum$index),]
  chr13_Uterus_hum<-rbind(chr13_Uterus_hum,b)
}
chr13_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr13_Uterus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr13_Uterus_hum_one2one<-chr13_Uterus_hum[match(intersect(chr13_pig2hum_one2one,chr13_Uterus_hum$gene_id),chr13_Uterus_hum$gene_id),]

chr14_Uterus_overloci_one2one<-chr14_Uterus_overloci[match(intersect(chr14_Uterus_genes,one2one_pig),chr14_Uterus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Uterus_chr14$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Uterus_chr14$hum_chr<-tmp$V1
Uterus_chr14$hum_loci<-tmp$V3
FM_Uterus_hum$index<-paste0(FM_Uterus_hum$chr,"_",FM_Uterus_hum$variant_pos)
Uterus_chr14$index<-paste0("chr",Uterus_chr14$hum_chr,"_",Uterus_chr14$hum_loci)
chr14_Uterus_hum<-NULL
for(i in 1:length(same_Uterus_chr14)){
  a<-Uterus_chr14[match(same_Uterus_chr14[i],Uterus_chr14$pig_pos),]
  b<-FM_Uterus_hum[match(intersect(a$index,FM_Uterus_hum$index),FM_Uterus_hum$index),]
  chr14_Uterus_hum<-rbind(chr14_Uterus_hum,b)
}
chr14_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr14_Uterus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr14_Uterus_hum_one2one<-chr14_Uterus_hum[match(intersect(chr14_pig2hum_one2one,chr14_Uterus_hum$gene_id),chr14_Uterus_hum$gene_id),]

chr15_Uterus_overloci_one2one<-chr15_Uterus_overloci[match(intersect(chr15_Uterus_genes,one2one_pig),chr15_Uterus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Uterus_chr15$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Uterus_chr15$hum_chr<-tmp$V1
Uterus_chr15$hum_loci<-tmp$V3
FM_Uterus_hum$index<-paste0(FM_Uterus_hum$chr,"_",FM_Uterus_hum$variant_pos)
Uterus_chr15$index<-paste0("chr",Uterus_chr15$hum_chr,"_",Uterus_chr15$hum_loci)
chr15_Uterus_hum<-NULL
for(i in 1:length(same_Uterus_chr15)){
  a<-Uterus_chr15[match(same_Uterus_chr15[i],Uterus_chr15$pig_pos),]
  b<-FM_Uterus_hum[match(intersect(a$index,FM_Uterus_hum$index),FM_Uterus_hum$index),]
  chr15_Uterus_hum<-rbind(chr15_Uterus_hum,b)
}
chr15_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr15_Uterus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr15_Uterus_hum_one2one<-chr15_Uterus_hum[match(intersect(chr15_pig2hum_one2one,chr15_Uterus_hum$gene_id),chr15_Uterus_hum$gene_id),]

chr16_Uterus_overloci_one2one<-chr16_Uterus_overloci[match(intersect(chr16_Uterus_genes,one2one_pig),chr16_Uterus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Uterus_chr16$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Uterus_chr16$hum_chr<-tmp$V1
Uterus_chr16$hum_loci<-tmp$V3
FM_Uterus_hum$index<-paste0(FM_Uterus_hum$chr,"_",FM_Uterus_hum$variant_pos)
Uterus_chr16$index<-paste0("chr",Uterus_chr16$hum_chr,"_",Uterus_chr16$hum_loci)
chr16_Uterus_hum<-NULL
for(i in 1:length(same_Uterus_chr16)){
  a<-Uterus_chr16[match(same_Uterus_chr16[i],Uterus_chr16$pig_pos),]
  b<-FM_Uterus_hum[match(intersect(a$index,FM_Uterus_hum$index),FM_Uterus_hum$index),]
  chr16_Uterus_hum<-rbind(chr16_Uterus_hum,b)
}
chr16_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr16_Uterus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr16_Uterus_hum_one2one<-chr16_Uterus_hum[match(intersect(chr16_pig2hum_one2one,chr16_Uterus_hum$gene_id),chr16_Uterus_hum$gene_id),]

chr17_Uterus_overloci_one2one<-chr17_Uterus_overloci[match(intersect(chr17_Uterus_genes,one2one_pig),chr17_Uterus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Uterus_chr17$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Uterus_chr17$hum_chr<-tmp$V1
Uterus_chr17$hum_loci<-tmp$V3
FM_Uterus_hum$index<-paste0(FM_Uterus_hum$chr,"_",FM_Uterus_hum$variant_pos)
Uterus_chr17$index<-paste0("chr",Uterus_chr17$hum_chr,"_",Uterus_chr17$hum_loci)
chr17_Uterus_hum<-NULL
for(i in 1:length(same_Uterus_chr17)){
  a<-Uterus_chr17[match(same_Uterus_chr17[i],Uterus_chr17$pig_pos),]
  b<-FM_Uterus_hum[match(intersect(a$index,FM_Uterus_hum$index),FM_Uterus_hum$index),]
  chr17_Uterus_hum<-rbind(chr17_Uterus_hum,b)
}
chr17_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr17_Uterus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr17_Uterus_hum_one2one<-chr17_Uterus_hum[match(intersect(chr17_pig2hum_one2one,chr17_Uterus_hum$gene_id),chr17_Uterus_hum$gene_id),]

chr18_Uterus_overloci_one2one<-chr18_Uterus_overloci[match(intersect(chr18_Uterus_genes,one2one_pig),chr18_Uterus_overloci$phenotype_id),]
tmp<-as.data.frame(strsplit(Uterus_chr18$hum_pos,split = "_"))
tmp<-as.data.frame(t(tmp))
Uterus_chr18$hum_chr<-tmp$V1
Uterus_chr18$hum_loci<-tmp$V3
FM_Uterus_hum$index<-paste0(FM_Uterus_hum$chr,"_",FM_Uterus_hum$variant_pos)
Uterus_chr18$index<-paste0("chr",Uterus_chr18$hum_chr,"_",Uterus_chr18$hum_loci)
chr18_Uterus_hum<-NULL
for(i in 1:length(same_Uterus_chr18)){
  a<-Uterus_chr18[match(same_Uterus_chr18[i],Uterus_chr18$pig_pos),]
  b<-FM_Uterus_hum[match(intersect(a$index,FM_Uterus_hum$index),FM_Uterus_hum$index),]
  chr18_Uterus_hum<-rbind(chr18_Uterus_hum,b)
}
chr18_pig2hum_one2one<-annotation$`Gene stable ID`[match(chr18_Uterus_overloci_one2one$phenotype_id,annotation$`Pig gene stable ID`)]
chr18_Uterus_hum_one2one<-chr18_Uterus_hum[match(intersect(chr18_pig2hum_one2one,chr18_Uterus_hum$gene_id),chr18_Uterus_hum$gene_id),]

Uterus_one2one_SNP_hum<-rbind(chr1_Uterus_hum_one2one,chr2_Uterus_hum_one2one,chr3_Uterus_hum_one2one,chr4_Uterus_hum_one2one,chr5_Uterus_hum_one2one,
                              chr6_Uterus_hum_one2one,chr7_Uterus_hum_one2one,chr8_Uterus_hum_one2one,chr9_Uterus_hum_one2one,chr10_Uterus_hum_one2one,
                              chr11_Uterus_hum_one2one,chr12_Uterus_hum_one2one,chr13_Uterus_hum_one2one,chr14_Uterus_hum_one2one,chr15_Uterus_hum_one2one,
                              chr16_Uterus_hum_one2one,chr17_Uterus_hum_one2one,chr18_Uterus_hum_one2one)

chr1_Uterus_pig_one2one<-chr1_Uterus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr1_Uterus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr1_Uterus_overloci_one2one$phenotype_id),1:9]
chr2_Uterus_pig_one2one<-chr2_Uterus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr2_Uterus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr2_Uterus_overloci_one2one$phenotype_id),1:9]
chr3_Uterus_pig_one2one<-chr3_Uterus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr3_Uterus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr3_Uterus_overloci_one2one$phenotype_id),1:9]
chr4_Uterus_pig_one2one<-chr4_Uterus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr4_Uterus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr4_Uterus_overloci_one2one$phenotype_id),1:9]
chr5_Uterus_pig_one2one<-chr5_Uterus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr5_Uterus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr5_Uterus_overloci_one2one$phenotype_id),1:9]
chr6_Uterus_pig_one2one<-chr6_Uterus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr6_Uterus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr6_Uterus_overloci_one2one$phenotype_id),1:9]
chr7_Uterus_pig_one2one<-chr7_Uterus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr7_Uterus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr7_Uterus_overloci_one2one$phenotype_id),1:9]
chr8_Uterus_pig_one2one<-chr8_Uterus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr8_Uterus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr8_Uterus_overloci_one2one$phenotype_id),1:9]
chr9_Uterus_pig_one2one<-chr9_Uterus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr9_Uterus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr9_Uterus_overloci_one2one$phenotype_id),1:9]
chr10_Uterus_pig_one2one<-chr10_Uterus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr10_Uterus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr10_Uterus_overloci_one2one$phenotype_id),1:9]
chr11_Uterus_pig_one2one<-chr11_Uterus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr11_Uterus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr11_Uterus_overloci_one2one$phenotype_id),1:9]
chr12_Uterus_pig_one2one<-chr12_Uterus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr12_Uterus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr12_Uterus_overloci_one2one$phenotype_id),1:9]
chr13_Uterus_pig_one2one<-chr13_Uterus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr13_Uterus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr13_Uterus_overloci_one2one$phenotype_id),1:9]
chr14_Uterus_pig_one2one<-chr14_Uterus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr14_Uterus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr14_Uterus_overloci_one2one$phenotype_id),1:9]
chr15_Uterus_pig_one2one<-chr15_Uterus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr15_Uterus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr15_Uterus_overloci_one2one$phenotype_id),1:9]
chr16_Uterus_pig_one2one<-chr16_Uterus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr16_Uterus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr16_Uterus_overloci_one2one$phenotype_id),1:9]
chr17_Uterus_pig_one2one<-chr17_Uterus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr17_Uterus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr17_Uterus_overloci_one2one$phenotype_id),1:9]
chr18_Uterus_pig_one2one<-chr18_Uterus_overloci_one2one[match(annotation$`Pig gene stable ID`[match(chr18_Uterus_hum_one2one$gene_id,annotation$`Gene stable ID`)],chr18_Uterus_overloci_one2one$phenotype_id),1:9]

Uterus_one2one_SNP_pig<-rbind(chr1_Uterus_pig_one2one,chr2_Uterus_pig_one2one,chr3_Uterus_pig_one2one,chr4_Uterus_pig_one2one,chr5_Uterus_pig_one2one,
                              chr6_Uterus_pig_one2one,chr7_Uterus_pig_one2one,chr8_Uterus_pig_one2one,chr9_Uterus_pig_one2one,chr10_Uterus_pig_one2one,
                              chr11_Uterus_pig_one2one,chr12_Uterus_pig_one2one,chr13_Uterus_pig_one2one,chr14_Uterus_pig_one2one,chr15_Uterus_pig_one2one,
                              chr16_Uterus_pig_one2one,chr17_Uterus_pig_one2one,chr18_Uterus_pig_one2one)

Uterus_SNP_sum<-array(NA,dim=c(nrow(Uterus_one2one_SNP_hum),2))
colnames(Uterus_SNP_sum)<-c("Human","Pig")
Uterus_SNP_sum<-as.data.frame(Uterus_SNP_sum)
Uterus_SNP_sum$Human<-Uterus_one2one_SNP_hum$slope / Uterus_one2one_SNP_hum$slope_se
Uterus_SNP_sum$Pig<-Uterus_one2one_SNP_pig$slope / Uterus_one2one_SNP_pig$slope_se
cor<-cor(abs(Uterus_SNP_sum$Human),abs(Uterus_SNP_sum$Pig))
p_val<-t.test(abs(Uterus_SNP_sum$Human),abs(Uterus_SNP_sum$Pig))
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/Uterus_SNP_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(Uterus_SNP_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The number of eGenes expressed in common tissues",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

save(Uterus_one2one_SNP_hum,Uterus_one2one_SNP_pig,file="/Users/baizhonghao/Downloads/human-pig GTEx/Uterus_SNP.Rdata")


SNP_sum<- rbind(Adipose_SNP_sum, Artery_SNP_sum, Blood_SNP_sum, Colon_SNP_sum, Frontal_cortex_SNP_sum, Heart_SNP_sum, Hypothalamus_SNP_sum,
                Ileum_SNP_sum, Kidney_SNP_sum, Liver_SNP_sum, Lung_SNP_sum, Muscle_SNP_sum, Ovary_SNP_sum, Pituitary_SNP_sum, Spleen_SNP_sum,
                Testis_SNP_sum, Uterus_SNP_sum)
SNP_sum$Human<-abs(SNP_sum$Human)
SNP_sum$Pig<-abs(SNP_sum$Pig)
cor<-cor(SNP_sum$Human,SNP_sum$Pig)
p_val<-cor.test(SNP_sum$Human,SNP_sum$Pig,method = "spearman")
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/all_SNP_plot_se.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(SNP_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("Same loci SNPs in all common tissues",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()


a<-as.data.frame(t(as.data.frame(strsplit(Adipose_one2one_SNP_pig$variant_id,split = "_"))))
Adipose_one2one_SNP_pig$ref<-a$V3
Adipose_one2one_SNP_pig$alt<-a$V4
Adipose_one2one_SNP_pig$var_index<-paste0(Adipose_one2one_SNP_pig$ref,"_",Adipose_one2one_SNP_pig$alt)
Adipose_one2one_SNP_pig$var_index1<-paste0(Adipose_one2one_SNP_pig$alt,"_",Adipose_one2one_SNP_pig$ref)
Adipose_one2one_SNP_hum$var_index<-paste0(Adipose_one2one_SNP_hum$ref,"_",Adipose_one2one_SNP_hum$alt)
Adipose_one2one_SNP_hum$var_index1<-paste0(Adipose_one2one_SNP_hum$alt,"_",Adipose_one2one_SNP_hum$ref)
Adipose_samevariant_hum<-NULL
Adipose_samevariant_pig<-NULL
for(i in 1:nrow(Adipose_one2one_SNP_hum)){
  if(Adipose_one2one_SNP_hum$var_index[i] == Adipose_one2one_SNP_pig$var_index[i] | Adipose_one2one_SNP_hum$var_index1[i] == Adipose_one2one_SNP_pig$var_index[i]){
    Adipose_samevariant_hum<-rbind(Adipose_samevariant_hum,Adipose_one2one_SNP_hum[i,])
    Adipose_samevariant_pig<-rbind(Adipose_samevariant_pig,Adipose_one2one_SNP_pig[i,])
  }
}

a<-as.data.frame(t(as.data.frame(strsplit(Artery_one2one_SNP_pig$variant_id,split = "_"))))
Artery_one2one_SNP_pig$ref<-a$V3
Artery_one2one_SNP_pig$alt<-a$V4
Artery_one2one_SNP_pig$var_index<-paste0(Artery_one2one_SNP_pig$ref,"_",Artery_one2one_SNP_pig$alt)
Artery_one2one_SNP_pig$var_index1<-paste0(Artery_one2one_SNP_pig$alt,"_",Artery_one2one_SNP_pig$ref)
Artery_one2one_SNP_hum$var_index<-paste0(Artery_one2one_SNP_hum$ref,"_",Artery_one2one_SNP_hum$alt)
Artery_one2one_SNP_hum$var_index1<-paste0(Artery_one2one_SNP_hum$alt,"_",Artery_one2one_SNP_hum$ref)
Artery_samevariant_hum<-NULL
Artery_samevariant_pig<-NULL
for(i in 1:nrow(Artery_one2one_SNP_hum)){
  if(Artery_one2one_SNP_hum$var_index[i] == Artery_one2one_SNP_pig$var_index[i] | Artery_one2one_SNP_hum$var_index1[i] == Artery_one2one_SNP_pig$var_index[i]){
    Artery_samevariant_hum<-rbind(Artery_samevariant_hum,Artery_one2one_SNP_hum[i,])
    Artery_samevariant_pig<-rbind(Artery_samevariant_pig,Artery_one2one_SNP_pig[i,])
  }
}

a<-as.data.frame(t(as.data.frame(strsplit(Blood_one2one_SNP_pig$variant_id,split = "_"))))
Blood_one2one_SNP_pig$ref<-a$V3
Blood_one2one_SNP_pig$alt<-a$V4
Blood_one2one_SNP_pig$var_index<-paste0(Blood_one2one_SNP_pig$ref,"_",Blood_one2one_SNP_pig$alt)
Blood_one2one_SNP_pig$var_index1<-paste0(Blood_one2one_SNP_pig$alt,"_",Blood_one2one_SNP_pig$ref)
Blood_one2one_SNP_hum$var_index<-paste0(Blood_one2one_SNP_hum$ref,"_",Blood_one2one_SNP_hum$alt)
Blood_one2one_SNP_hum$var_index1<-paste0(Blood_one2one_SNP_hum$alt,"_",Blood_one2one_SNP_hum$ref)
Blood_samevariant_hum<-NULL
Blood_samevariant_pig<-NULL
for(i in 1:nrow(Blood_one2one_SNP_hum)){
  if(Blood_one2one_SNP_hum$var_index[i] == Blood_one2one_SNP_pig$var_index[i] | Blood_one2one_SNP_hum$var_index1[i] == Blood_one2one_SNP_pig$var_index[i]){
    Blood_samevariant_hum<-rbind(Blood_samevariant_hum,Blood_one2one_SNP_hum[i,])
    Blood_samevariant_pig<-rbind(Blood_samevariant_pig,Blood_one2one_SNP_pig[i,])
  }
}

a<-as.data.frame(t(as.data.frame(strsplit(Colon_one2one_SNP_pig$variant_id,split = "_"))))
Colon_one2one_SNP_pig$ref<-a$V3
Colon_one2one_SNP_pig$alt<-a$V4
Colon_one2one_SNP_pig$var_index<-paste0(Colon_one2one_SNP_pig$ref,"_",Colon_one2one_SNP_pig$alt)
Colon_one2one_SNP_pig$var_index1<-paste0(Colon_one2one_SNP_pig$alt,"_",Colon_one2one_SNP_pig$ref)
Colon_one2one_SNP_hum$var_index<-paste0(Colon_one2one_SNP_hum$ref,"_",Colon_one2one_SNP_hum$alt)
Colon_one2one_SNP_hum$var_index1<-paste0(Colon_one2one_SNP_hum$alt,"_",Colon_one2one_SNP_hum$ref)
Colon_samevariant_hum<-NULL
Colon_samevariant_pig<-NULL
for(i in 1:nrow(Colon_one2one_SNP_hum)){
  if(Colon_one2one_SNP_hum$var_index[i] == Colon_one2one_SNP_pig$var_index[i] | Colon_one2one_SNP_hum$var_index1[i] == Colon_one2one_SNP_pig$var_index[i]){
    Colon_samevariant_hum<-rbind(Colon_samevariant_hum,Colon_one2one_SNP_hum[i,])
    Colon_samevariant_pig<-rbind(Colon_samevariant_pig,Colon_one2one_SNP_pig[i,])
  }
}

a<-as.data.frame(t(as.data.frame(strsplit(Frontal_cortex_one2one_SNP_pig$variant_id,split = "_"))))
Frontal_cortex_one2one_SNP_pig$ref<-a$V3
Frontal_cortex_one2one_SNP_pig$alt<-a$V4
Frontal_cortex_one2one_SNP_pig$var_index<-paste0(Frontal_cortex_one2one_SNP_pig$ref,"_",Frontal_cortex_one2one_SNP_pig$alt)
Frontal_cortex_one2one_SNP_pig$var_index1<-paste0(Frontal_cortex_one2one_SNP_pig$alt,"_",Frontal_cortex_one2one_SNP_pig$ref)
Frontal_cortex_one2one_SNP_hum$var_index<-paste0(Frontal_cortex_one2one_SNP_hum$ref,"_",Frontal_cortex_one2one_SNP_hum$alt)
Frontal_cortex_one2one_SNP_hum$var_index1<-paste0(Frontal_cortex_one2one_SNP_hum$alt,"_",Frontal_cortex_one2one_SNP_hum$ref)
Frontal_cortex_samevariant_hum<-NULL
Frontal_cortex_samevariant_pig<-NULL
for(i in 1:nrow(Frontal_cortex_one2one_SNP_hum)){
  if(Frontal_cortex_one2one_SNP_hum$var_index[i] == Frontal_cortex_one2one_SNP_pig$var_index[i] | Frontal_cortex_one2one_SNP_hum$var_index1[i] == Frontal_cortex_one2one_SNP_pig$var_index[i]){
    Frontal_cortex_samevariant_hum<-rbind(Frontal_cortex_samevariant_hum,Frontal_cortex_one2one_SNP_hum[i,])
    Frontal_cortex_samevariant_pig<-rbind(Frontal_cortex_samevariant_pig,Frontal_cortex_one2one_SNP_pig[i,])
  }
}

a<-as.data.frame(t(as.data.frame(strsplit(Heart_one2one_SNP_pig$variant_id,split = "_"))))
Heart_one2one_SNP_pig$ref<-a$V3
Heart_one2one_SNP_pig$alt<-a$V4
Heart_one2one_SNP_pig$var_index<-paste0(Heart_one2one_SNP_pig$ref,"_",Heart_one2one_SNP_pig$alt)
Heart_one2one_SNP_pig$var_index1<-paste0(Heart_one2one_SNP_pig$alt,"_",Heart_one2one_SNP_pig$ref)
Heart_one2one_SNP_hum$var_index<-paste0(Heart_one2one_SNP_hum$ref,"_",Heart_one2one_SNP_hum$alt)
Heart_one2one_SNP_hum$var_index1<-paste0(Heart_one2one_SNP_hum$alt,"_",Heart_one2one_SNP_hum$ref)
Heart_samevariant_hum<-NULL
Heart_samevariant_pig<-NULL
for(i in 1:nrow(Heart_one2one_SNP_hum)){
  if(Heart_one2one_SNP_hum$var_index[i] == Heart_one2one_SNP_pig$var_index[i] | Heart_one2one_SNP_hum$var_index1[i] == Heart_one2one_SNP_pig$var_index[i]){
    Heart_samevariant_hum<-rbind(Heart_samevariant_hum,Heart_one2one_SNP_hum[i,])
    Heart_samevariant_pig<-rbind(Heart_samevariant_pig,Heart_one2one_SNP_pig[i,])
  }
}

a<-as.data.frame(t(as.data.frame(strsplit(Hypothalamus_one2one_SNP_pig$variant_id,split = "_"))))
Hypothalamus_one2one_SNP_pig$ref<-a$V3
Hypothalamus_one2one_SNP_pig$alt<-a$V4
Hypothalamus_one2one_SNP_pig$var_index<-paste0(Hypothalamus_one2one_SNP_pig$ref,"_",Hypothalamus_one2one_SNP_pig$alt)
Hypothalamus_one2one_SNP_pig$var_index1<-paste0(Hypothalamus_one2one_SNP_pig$alt,"_",Hypothalamus_one2one_SNP_pig$ref)
Hypothalamus_one2one_SNP_hum$var_index<-paste0(Hypothalamus_one2one_SNP_hum$ref,"_",Hypothalamus_one2one_SNP_hum$alt)
Hypothalamus_one2one_SNP_hum$var_index1<-paste0(Hypothalamus_one2one_SNP_hum$alt,"_",Hypothalamus_one2one_SNP_hum$ref)
Hypothalamus_samevariant_hum<-NULL
Hypothalamus_samevariant_pig<-NULL
for(i in 1:nrow(Hypothalamus_one2one_SNP_hum)){
  if(Hypothalamus_one2one_SNP_hum$var_index[i] == Hypothalamus_one2one_SNP_pig$var_index[i] | Hypothalamus_one2one_SNP_hum$var_index1[i] == Hypothalamus_one2one_SNP_pig$var_index[i]){
    Hypothalamus_samevariant_hum<-rbind(Hypothalamus_samevariant_hum,Hypothalamus_one2one_SNP_hum[i,])
    Hypothalamus_samevariant_pig<-rbind(Hypothalamus_samevariant_pig,Hypothalamus_one2one_SNP_pig[i,])
  }
}

a<-as.data.frame(t(as.data.frame(strsplit(Ileum_one2one_SNP_pig$variant_id,split = "_"))))
Ileum_one2one_SNP_pig$ref<-a$V3
Ileum_one2one_SNP_pig$alt<-a$V4
Ileum_one2one_SNP_pig$var_index<-paste0(Ileum_one2one_SNP_pig$ref,"_",Ileum_one2one_SNP_pig$alt)
Ileum_one2one_SNP_pig$var_index1<-paste0(Ileum_one2one_SNP_pig$alt,"_",Ileum_one2one_SNP_pig$ref)
Ileum_one2one_SNP_hum$var_index<-paste0(Ileum_one2one_SNP_hum$ref,"_",Ileum_one2one_SNP_hum$alt)
Ileum_one2one_SNP_hum$var_index1<-paste0(Ileum_one2one_SNP_hum$alt,"_",Ileum_one2one_SNP_hum$ref)
Ileum_samevariant_hum<-NULL
Ileum_samevariant_pig<-NULL
for(i in 1:nrow(Ileum_one2one_SNP_hum)){
  if(Ileum_one2one_SNP_hum$var_index[i] == Ileum_one2one_SNP_pig$var_index[i] | Ileum_one2one_SNP_hum$var_index1[i] == Ileum_one2one_SNP_pig$var_index[i]){
    Ileum_samevariant_hum<-rbind(Ileum_samevariant_hum,Ileum_one2one_SNP_hum[i,])
    Ileum_samevariant_pig<-rbind(Ileum_samevariant_pig,Ileum_one2one_SNP_pig[i,])
  }
}

a<-as.data.frame(t(as.data.frame(strsplit(Kidney_one2one_SNP_pig$variant_id,split = "_"))))
Kidney_one2one_SNP_pig$ref<-a$V3
Kidney_one2one_SNP_pig$alt<-a$V4
Kidney_one2one_SNP_pig$var_index<-paste0(Kidney_one2one_SNP_pig$ref,"_",Kidney_one2one_SNP_pig$alt)
Kidney_one2one_SNP_pig$var_index1<-paste0(Kidney_one2one_SNP_pig$alt,"_",Kidney_one2one_SNP_pig$ref)
Kidney_one2one_SNP_hum$var_index<-paste0(Kidney_one2one_SNP_hum$ref,"_",Kidney_one2one_SNP_hum$alt)
Kidney_one2one_SNP_hum$var_index1<-paste0(Kidney_one2one_SNP_hum$alt,"_",Kidney_one2one_SNP_hum$ref)
Kidney_samevariant_hum<-NULL
Kidney_samevariant_pig<-NULL
for(i in 1:nrow(Kidney_one2one_SNP_hum)){
  if(Kidney_one2one_SNP_hum$var_index[i] == Kidney_one2one_SNP_pig$var_index[i] | Kidney_one2one_SNP_hum$var_index1[i] == Kidney_one2one_SNP_pig$var_index[i]){
    Kidney_samevariant_hum<-rbind(Kidney_samevariant_hum,Kidney_one2one_SNP_hum[i,])
    Kidney_samevariant_pig<-rbind(Kidney_samevariant_pig,Kidney_one2one_SNP_pig[i,])
  }
}

a<-as.data.frame(t(as.data.frame(strsplit(Liver_one2one_SNP_pig$variant_id,split = "_"))))
Liver_one2one_SNP_pig$ref<-a$V3
Liver_one2one_SNP_pig$alt<-a$V4
Liver_one2one_SNP_pig$var_index<-paste0(Liver_one2one_SNP_pig$ref,"_",Liver_one2one_SNP_pig$alt)
Liver_one2one_SNP_pig$var_index1<-paste0(Liver_one2one_SNP_pig$alt,"_",Liver_one2one_SNP_pig$ref)
Liver_one2one_SNP_hum$var_index<-paste0(Liver_one2one_SNP_hum$ref,"_",Liver_one2one_SNP_hum$alt)
Liver_one2one_SNP_hum$var_index1<-paste0(Liver_one2one_SNP_hum$alt,"_",Liver_one2one_SNP_hum$ref)
Liver_samevariant_hum<-NULL
Liver_samevariant_pig<-NULL
for(i in 1:nrow(Liver_one2one_SNP_hum)){
  if(Liver_one2one_SNP_hum$var_index[i] == Liver_one2one_SNP_pig$var_index[i] | Liver_one2one_SNP_hum$var_index1[i] == Liver_one2one_SNP_pig$var_index[i]){
    Liver_samevariant_hum<-rbind(Liver_samevariant_hum,Liver_one2one_SNP_hum[i,])
    Liver_samevariant_pig<-rbind(Liver_samevariant_pig,Liver_one2one_SNP_pig[i,])
  }
}

a<-as.data.frame(t(as.data.frame(strsplit(Lung_one2one_SNP_pig$variant_id,split = "_"))))
Lung_one2one_SNP_pig$ref<-a$V3
Lung_one2one_SNP_pig$alt<-a$V4
Lung_one2one_SNP_pig$var_index<-paste0(Lung_one2one_SNP_pig$ref,"_",Lung_one2one_SNP_pig$alt)
Lung_one2one_SNP_pig$var_index1<-paste0(Lung_one2one_SNP_pig$alt,"_",Lung_one2one_SNP_pig$ref)
Lung_one2one_SNP_hum$var_index<-paste0(Lung_one2one_SNP_hum$ref,"_",Lung_one2one_SNP_hum$alt)
Lung_one2one_SNP_hum$var_index1<-paste0(Lung_one2one_SNP_hum$alt,"_",Lung_one2one_SNP_hum$ref)
Lung_samevariant_hum<-NULL
Lung_samevariant_pig<-NULL
for(i in 1:nrow(Lung_one2one_SNP_hum)){
  if(Lung_one2one_SNP_hum$var_index[i] == Lung_one2one_SNP_pig$var_index[i] | Lung_one2one_SNP_hum$var_index1[i] == Lung_one2one_SNP_pig$var_index[i]){
    Lung_samevariant_hum<-rbind(Lung_samevariant_hum,Lung_one2one_SNP_hum[i,])
    Lung_samevariant_pig<-rbind(Lung_samevariant_pig,Lung_one2one_SNP_pig[i,])
  }
}

a<-as.data.frame(t(as.data.frame(strsplit(Muscle_one2one_SNP_pig$variant_id,split = "_"))))
Muscle_one2one_SNP_pig$ref<-a$V3
Muscle_one2one_SNP_pig$alt<-a$V4
Muscle_one2one_SNP_pig$var_index<-paste0(Muscle_one2one_SNP_pig$ref,"_",Muscle_one2one_SNP_pig$alt)
Muscle_one2one_SNP_pig$var_index1<-paste0(Muscle_one2one_SNP_pig$alt,"_",Muscle_one2one_SNP_pig$ref)
Muscle_one2one_SNP_hum$var_index<-paste0(Muscle_one2one_SNP_hum$ref,"_",Muscle_one2one_SNP_hum$alt)
Muscle_one2one_SNP_hum$var_index1<-paste0(Muscle_one2one_SNP_hum$alt,"_",Muscle_one2one_SNP_hum$ref)
Muscle_samevariant_hum<-NULL
Muscle_samevariant_pig<-NULL
for(i in 1:nrow(Muscle_one2one_SNP_hum)){
  if(Muscle_one2one_SNP_hum$var_index[i] == Muscle_one2one_SNP_pig$var_index[i] | Muscle_one2one_SNP_hum$var_index1[i] == Muscle_one2one_SNP_pig$var_index[i]){
    Muscle_samevariant_hum<-rbind(Muscle_samevariant_hum,Muscle_one2one_SNP_hum[i,])
    Muscle_samevariant_pig<-rbind(Muscle_samevariant_pig,Muscle_one2one_SNP_pig[i,])
  }
}

a<-as.data.frame(t(as.data.frame(strsplit(Ovary_one2one_SNP_pig$variant_id,split = "_"))))
Ovary_one2one_SNP_pig$ref<-a$V3
Ovary_one2one_SNP_pig$alt<-a$V4
Ovary_one2one_SNP_pig$var_index<-paste0(Ovary_one2one_SNP_pig$ref,"_",Ovary_one2one_SNP_pig$alt)
Ovary_one2one_SNP_pig$var_index1<-paste0(Ovary_one2one_SNP_pig$alt,"_",Ovary_one2one_SNP_pig$ref)
Ovary_one2one_SNP_hum$var_index<-paste0(Ovary_one2one_SNP_hum$ref,"_",Ovary_one2one_SNP_hum$alt)
Ovary_one2one_SNP_hum$var_index1<-paste0(Ovary_one2one_SNP_hum$alt,"_",Ovary_one2one_SNP_hum$ref)
Ovary_samevariant_hum<-NULL
Ovary_samevariant_pig<-NULL
for(i in 1:nrow(Ovary_one2one_SNP_hum)){
  if(Ovary_one2one_SNP_hum$var_index[i] == Ovary_one2one_SNP_pig$var_index[i] | Ovary_one2one_SNP_hum$var_index1[i] == Ovary_one2one_SNP_pig$var_index[i]){
    Ovary_samevariant_hum<-rbind(Ovary_samevariant_hum,Ovary_one2one_SNP_hum[i,])
    Ovary_samevariant_pig<-rbind(Ovary_samevariant_pig,Ovary_one2one_SNP_pig[i,])
  }
}

a<-as.data.frame(t(as.data.frame(strsplit(Pituitary_one2one_SNP_pig$variant_id,split = "_"))))
Pituitary_one2one_SNP_pig$ref<-a$V3
Pituitary_one2one_SNP_pig$alt<-a$V4
Pituitary_one2one_SNP_pig$var_index<-paste0(Pituitary_one2one_SNP_pig$ref,"_",Pituitary_one2one_SNP_pig$alt)
Pituitary_one2one_SNP_pig$var_index1<-paste0(Pituitary_one2one_SNP_pig$alt,"_",Pituitary_one2one_SNP_pig$ref)
Pituitary_one2one_SNP_hum$var_index<-paste0(Pituitary_one2one_SNP_hum$ref,"_",Pituitary_one2one_SNP_hum$alt)
Pituitary_one2one_SNP_hum$var_index1<-paste0(Pituitary_one2one_SNP_hum$alt,"_",Pituitary_one2one_SNP_hum$ref)
Pituitary_samevariant_hum<-NULL
Pituitary_samevariant_pig<-NULL
for(i in 1:nrow(Pituitary_one2one_SNP_hum)){
  if(Pituitary_one2one_SNP_hum$var_index[i] == Pituitary_one2one_SNP_pig$var_index[i] | Pituitary_one2one_SNP_hum$var_index1[i] == Pituitary_one2one_SNP_pig$var_index[i]){
    Pituitary_samevariant_hum<-rbind(Pituitary_samevariant_hum,Pituitary_one2one_SNP_hum[i,])
    Pituitary_samevariant_pig<-rbind(Pituitary_samevariant_pig,Pituitary_one2one_SNP_pig[i,])
  }
}

a<-as.data.frame(t(as.data.frame(strsplit(Spleen_one2one_SNP_pig$variant_id,split = "_"))))
Spleen_one2one_SNP_pig$ref<-a$V3
Spleen_one2one_SNP_pig$alt<-a$V4
Spleen_one2one_SNP_pig$var_index<-paste0(Spleen_one2one_SNP_pig$ref,"_",Spleen_one2one_SNP_pig$alt)
Spleen_one2one_SNP_pig$var_index1<-paste0(Spleen_one2one_SNP_pig$alt,"_",Spleen_one2one_SNP_pig$ref)
Spleen_one2one_SNP_hum$var_index<-paste0(Spleen_one2one_SNP_hum$ref,"_",Spleen_one2one_SNP_hum$alt)
Spleen_one2one_SNP_hum$var_index1<-paste0(Spleen_one2one_SNP_hum$alt,"_",Spleen_one2one_SNP_hum$ref)
Spleen_samevariant_hum<-NULL
Spleen_samevariant_pig<-NULL
for(i in 1:nrow(Spleen_one2one_SNP_hum)){
  if(Spleen_one2one_SNP_hum$var_index[i] == Spleen_one2one_SNP_pig$var_index[i] | Spleen_one2one_SNP_hum$var_index1[i] == Spleen_one2one_SNP_pig$var_index[i]){
    Spleen_samevariant_hum<-rbind(Spleen_samevariant_hum,Spleen_one2one_SNP_hum[i,])
    Spleen_samevariant_pig<-rbind(Spleen_samevariant_pig,Spleen_one2one_SNP_pig[i,])
  }
}

a<-as.data.frame(t(as.data.frame(strsplit(Testis_one2one_SNP_pig$variant_id,split = "_"))))
Testis_one2one_SNP_pig$ref<-a$V3
Testis_one2one_SNP_pig$alt<-a$V4
Testis_one2one_SNP_pig$var_index<-paste0(Testis_one2one_SNP_pig$ref,"_",Testis_one2one_SNP_pig$alt)
Testis_one2one_SNP_pig$var_index1<-paste0(Testis_one2one_SNP_pig$alt,"_",Testis_one2one_SNP_pig$ref)
Testis_one2one_SNP_hum$var_index<-paste0(Testis_one2one_SNP_hum$ref,"_",Testis_one2one_SNP_hum$alt)
Testis_one2one_SNP_hum$var_index1<-paste0(Testis_one2one_SNP_hum$alt,"_",Testis_one2one_SNP_hum$ref)
Testis_samevariant_hum<-NULL
Testis_samevariant_pig<-NULL
for(i in 1:nrow(Testis_one2one_SNP_hum)){
  if(Testis_one2one_SNP_hum$var_index[i] == Testis_one2one_SNP_pig$var_index[i] | Testis_one2one_SNP_hum$var_index1[i] == Testis_one2one_SNP_pig$var_index[i]){
    Testis_samevariant_hum<-rbind(Testis_samevariant_hum,Testis_one2one_SNP_hum[i,])
    Testis_samevariant_pig<-rbind(Testis_samevariant_pig,Testis_one2one_SNP_pig[i,])
  }
}

a<-as.data.frame(t(as.data.frame(strsplit(Uterus_one2one_SNP_pig$variant_id,split = "_"))))
Uterus_one2one_SNP_pig$ref<-a$V3
Uterus_one2one_SNP_pig$alt<-a$V4
Uterus_one2one_SNP_pig$var_index<-paste0(Uterus_one2one_SNP_pig$ref,"_",Uterus_one2one_SNP_pig$alt)
Uterus_one2one_SNP_pig$var_index1<-paste0(Uterus_one2one_SNP_pig$alt,"_",Uterus_one2one_SNP_pig$ref)
Uterus_one2one_SNP_hum$var_index<-paste0(Uterus_one2one_SNP_hum$ref,"_",Uterus_one2one_SNP_hum$alt)
Uterus_one2one_SNP_hum$var_index1<-paste0(Uterus_one2one_SNP_hum$alt,"_",Uterus_one2one_SNP_hum$ref)
Uterus_samevariant_hum<-NULL
Uterus_samevariant_pig<-NULL
for(i in 1:nrow(Uterus_one2one_SNP_hum)){
  if(Uterus_one2one_SNP_hum$var_index[i] == Uterus_one2one_SNP_pig$var_index[i] | Uterus_one2one_SNP_hum$var_index1[i] == Uterus_one2one_SNP_pig$var_index[i]){
    Uterus_samevariant_hum<-rbind(Uterus_samevariant_hum,Uterus_one2one_SNP_hum[i,])
    Uterus_samevariant_pig<-rbind(Uterus_samevariant_pig,Uterus_one2one_SNP_pig[i,])
  }
}

samevariant_sum_hum<-rbind(Adipose_samevariant_hum, Artery_samevariant_hum, Blood_samevariant_hum, Colon_samevariant_hum, Frontal_cortex_samevariant_hum, Heart_samevariant_hum, Hypothalamus_samevariant_hum,
                           Ileum_samevariant_hum, Kidney_samevariant_hum, Liver_samevariant_hum, Lung_samevariant_hum, Muscle_samevariant_hum, Ovary_samevariant_hum, Pituitary_samevariant_hum, Spleen_samevariant_hum,
                           Testis_samevariant_hum, Uterus_samevariant_hum)
samevariant_sum_pig<-rbind(Adipose_samevariant_pig[,1:9], Artery_samevariant_pig[,1:9], Blood_samevariant_pig[,1:9], Colon_samevariant_pig[,1:9], Frontal_cortex_samevariant_pig[,1:9], Heart_samevariant_pig[,1:9], Hypothalamus_samevariant_pig[,1:9],
                           Ileum_samevariant_pig[,1:9], Kidney_samevariant_pig[,1:9], Liver_samevariant_pig[,1:9], Lung_samevariant_pig[,1:9], Muscle_samevariant_pig[,1:9], Ovary_samevariant_pig[,1:9], Pituitary_samevariant_pig[,1:9], Spleen_samevariant_pig[,1:9],
                           Testis_samevariant_pig[,1:9], Uterus_samevariant_pig[,1:9])

sameloci_sum_pig<-rbind(Adipose_one2one_SNP_pig[,1:9], Artery_one2one_SNP_pig[,1:9], Blood_one2one_SNP_pig[,1:9], Colon_one2one_SNP_pig[,1:9], Frontal_cortex_one2one_SNP_pig[,1:9], Heart_one2one_SNP_pig[,1:9], Hypothalamus_one2one_SNP_pig[,1:9],
                           Ileum_one2one_SNP_pig[,1:9], Kidney_one2one_SNP_pig[,1:9], Liver_one2one_SNP_pig[,1:9], Lung_one2one_SNP_pig[,1:9], Muscle_one2one_SNP_pig[,1:9], Ovary_one2one_SNP_pig[,1:9], Pituitary_one2one_SNP_pig[,1:9], Spleen_one2one_SNP_pig[,1:9],
                           Testis_one2one_SNP_pig[,1:9], Uterus_one2one_SNP_pig[,1:9])

sameloci_sum_hum<-rbind(Adipose_one2one_SNP_hum, Artery_one2one_SNP_hum, Blood_one2one_SNP_hum, Colon_one2one_SNP_hum, Frontal_cortex_one2one_SNP_hum, Heart_one2one_SNP_hum, Hypothalamus_one2one_SNP_hum,
                        Ileum_one2one_SNP_hum, Kidney_one2one_SNP_hum, Liver_one2one_SNP_hum, Lung_one2one_SNP_hum, Muscle_one2one_SNP_hum, Ovary_one2one_SNP_hum, Pituitary_one2one_SNP_hum, Spleen_one2one_SNP_hum,
                        Testis_one2one_SNP_hum, Uterus_one2one_SNP_hum)

samevariant_sum<-array(NA,dim=c(nrow(samevariant_sum_hum),2))
colnames(samevariant_sum)<-c("Human","Pig")
samevariant_sum<-as.data.frame(samevariant_sum)
samevariant_sum$Human<-samevariant_sum_hum$slope / samevariant_sum_hum$slope_se
samevariant_sum$Pig<-samevariant_sum_pig$slope / samevariant_sum_pig$slope_se
samevariant_sum$Human<-abs(samevariant_sum$Human)
samevariant_sum$Pig<-abs(samevariant_sum$Pig)
cor<-cor(samevariant_sum$Human,samevariant_sum$Pig)
p_val<-cor.test(samevariant_sum$Human,samevariant_sum$Pig,method = "spearman")
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/samevariant_SNP_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(samevariant_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The effect size(slope/se) of samevariant SNPs",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

samevariant_file<-dir("/Users/baizhonghao/Downloads/output_aFC_one2one_SNP/one2one",full.names = T)
samevariant_file
samevariant<-substr(samevariant_file,nchar("/Users/baizhonghao/Downloads/output_aFC_one2one_SNP/samevariant/")+1,nchar(samevariant_file))
samevariant

samevariant_sumfc<-NULL
for(i in 1:length(samevariant)){
  a<-as.data.frame(fread(samevariant_file[i]))
  samevariant_sumfc<-rbind(samevariant_sumfc,a)
}

afc_sum<-array(NA,dim=c(nrow(sameloci_sum_hum),2))
colnames(afc_sum)<-c("Human","Pig")
afc_sum<-as.data.frame(afc_sum)
afc_sum$Human<-abs(sameloci_sum_hum$log2_aFC)
afc_sum$Pig<-abs(samevariant_sumfc$log2_aFC)
cor<-cor(afc_sum$Human,afc_sum$Pig)
p_val<-cor.test(afc_sum$Human,afc_sum$Pig,method="spearman")
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/sameloci_SNP_FCplot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(afc_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The effect size(abs(log2_aFC)) of sameloci SNPs",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

factor_Adipose<-as.data.frame(fread("/Users/baizhonghao/Downloads/ref_factor/Adipose.ref_factor.txt"))
factor_Artery<-as.data.frame(fread("/Users/baizhonghao/Downloads/ref_factor/Artery.ref_factor.txt"))
factor_Blood<-as.data.frame(fread("/Users/baizhonghao/Downloads/ref_factor/Blood.ref_factor.txt"))
factor_Colon<-as.data.frame(fread("/Users/baizhonghao/Downloads/ref_factor/Colon.ref_factor.txt"))
factor_Frontal_cortex<-as.data.frame(fread("/Users/baizhonghao/Downloads/ref_factor/Frontal_cortex.ref_factor.txt"))
factor_Heart<-as.data.frame(fread("/Users/baizhonghao/Downloads/ref_factor/Heart.ref_factor.txt"))
factor_Hypothalamus<-as.data.frame(fread("/Users/baizhonghao/Downloads/ref_factor/Hypothalamus.ref_factor.txt"))
factor_Ileum<-as.data.frame(fread("/Users/baizhonghao/Downloads/ref_factor/Ileum.ref_factor.txt"))
factor_Kidney<-as.data.frame(fread("/Users/baizhonghao/Downloads/ref_factor/Kidney.ref_factor.txt"))
factor_Liver<-as.data.frame(fread("/Users/baizhonghao/Downloads/ref_factor/Liver.ref_factor.txt"))
factor_Lung<-as.data.frame(fread("/Users/baizhonghao/Downloads/ref_factor/Lung.ref_factor.txt"))
factor_Muscle<-as.data.frame(fread("/Users/baizhonghao/Downloads/ref_factor/Muscle.ref_factor.txt"))
factor_Ovary<-as.data.frame(fread("/Users/baizhonghao/Downloads/ref_factor/Ovary.ref_factor.txt"))
factor_Pituitary<-as.data.frame(fread("/Users/baizhonghao/Downloads/ref_factor/Pituitary.ref_factor.txt"))
factor_Spleen<-as.data.frame(fread("/Users/baizhonghao/Downloads/ref_factor/Spleen.ref_factor.txt"))
factor_Testis<-as.data.frame(fread("/Users/baizhonghao/Downloads/ref_factor/Testis.ref_factor.txt"))
factor_Uterus<-as.data.frame(fread("/Users/baizhonghao/Downloads/ref_factor/Uterus.ref_factor.txt"))
Adipose_samevariant_pig$factor<-factor_Adipose$V2[match(Adipose_samevariant_pig$variant_id,factor_Adipose$V1)]
Artery_samevariant_pig$factor<-factor_Artery$V2[match(Artery_samevariant_pig$variant_id,factor_Artery$V1)]
Blood_samevariant_pig$factor<-factor_Blood$V2[match(Blood_samevariant_pig$variant_id,factor_Blood$V1)]
Colon_samevariant_pig$factor<-factor_Colon$V2[match(Colon_samevariant_pig$variant_id,factor_Colon$V1)]
Frontal_cortex_samevariant_pig$factor<-factor_Frontal_cortex$V2[match(Frontal_cortex_samevariant_pig$variant_id,factor_Frontal_cortex$V1)]
Heart_samevariant_pig$factor<-factor_Heart$V2[match(Heart_samevariant_pig$variant_id,factor_Heart$V1)]
Hypothalamus_samevariant_pig$factor<-factor_Hypothalamus$V2[match(Hypothalamus_samevariant_pig$variant_id,factor_Hypothalamus$V1)]
Ileum_samevariant_pig$factor<-factor_Ileum$V2[match(Ileum_samevariant_pig$variant_id,factor_Ileum$V1)]
Kidney_samevariant_pig$factor<-factor_Kidney$V2[match(Kidney_samevariant_pig$variant_id,factor_Kidney$V1)]
Liver_samevariant_pig$factor<-factor_Liver$V2[match(Liver_samevariant_pig$variant_id,factor_Liver$V1)]
Lung_samevariant_pig$factor<-factor_Lung$V2[match(Lung_samevariant_pig$variant_id,factor_Lung$V1)]
Muscle_samevariant_pig$factor<-factor_Muscle$V2[match(Muscle_samevariant_pig$variant_id,factor_Muscle$V1)]
Ovary_samevariant_pig$factor<-factor_Ovary$V2[match(Ovary_samevariant_pig$variant_id,factor_Ovary$V1)]
Pituitary_samevariant_pig$factor<-factor_Pituitary$V2[match(Pituitary_samevariant_pig$variant_id,factor_Pituitary$V1)]
Spleen_samevariant_pig$factor<-factor_Spleen$V2[match(Spleen_samevariant_pig$variant_id,factor_Spleen$V1)]
Testis_samevariant_pig$factor<-factor_Testis$V2[match(Testis_samevariant_pig$variant_id,factor_Testis$V1)]
Uterus_samevariant_pig$factor<-factor_Uterus$V2[match(Uterus_samevariant_pig$variant_id,factor_Uterus$V1)]

Adipose_one2one_SNP_pig$factor<-factor_Adipose$V2[match(Adipose_one2one_SNP_pig$variant_id,factor_Adipose$V1)]
Artery_one2one_SNP_pig$factor<-factor_Artery$V2[match(Artery_one2one_SNP_pig$variant_id,factor_Artery$V1)]
Blood_one2one_SNP_pig$factor<-factor_Blood$V2[match(Blood_one2one_SNP_pig$variant_id,factor_Blood$V1)]
Colon_one2one_SNP_pig$factor<-factor_Colon$V2[match(Colon_one2one_SNP_pig$variant_id,factor_Colon$V1)]
Frontal_cortex_one2one_SNP_pig$factor<-factor_Frontal_cortex$V2[match(Frontal_cortex_one2one_SNP_pig$variant_id,factor_Frontal_cortex$V1)]
Heart_one2one_SNP_pig$factor<-factor_Heart$V2[match(Heart_one2one_SNP_pig$variant_id,factor_Heart$V1)]
Hypothalamus_one2one_SNP_pig$factor<-factor_Hypothalamus$V2[match(Hypothalamus_one2one_SNP_pig$variant_id,factor_Hypothalamus$V1)]
Ileum_one2one_SNP_pig$factor<-factor_Ileum$V2[match(Ileum_one2one_SNP_pig$variant_id,factor_Ileum$V1)]
Kidney_one2one_SNP_pig$factor<-factor_Kidney$V2[match(Kidney_one2one_SNP_pig$variant_id,factor_Kidney$V1)]
Liver_one2one_SNP_pig$factor<-factor_Liver$V2[match(Liver_one2one_SNP_pig$variant_id,factor_Liver$V1)]
Lung_one2one_SNP_pig$factor<-factor_Lung$V2[match(Lung_one2one_SNP_pig$variant_id,factor_Lung$V1)]
Muscle_one2one_SNP_pig$factor<-factor_Muscle$V2[match(Muscle_one2one_SNP_pig$variant_id,factor_Muscle$V1)]
Ovary_one2one_SNP_pig$factor<-factor_Ovary$V2[match(Ovary_one2one_SNP_pig$variant_id,factor_Ovary$V1)]
Pituitary_one2one_SNP_pig$factor<-factor_Pituitary$V2[match(Pituitary_one2one_SNP_pig$variant_id,factor_Pituitary$V1)]
Spleen_one2one_SNP_pig$factor<-factor_Spleen$V2[match(Spleen_one2one_SNP_pig$variant_id,factor_Spleen$V1)]
Testis_one2one_SNP_pig$factor<-factor_Testis$V2[match(Testis_one2one_SNP_pig$variant_id,factor_Testis$V1)]
Uterus_one2one_SNP_pig$factor<-factor_Uterus$V2[match(Uterus_one2one_SNP_pig$variant_id,factor_Uterus$V1)]

library(xlsx)
for(i in 1:length(tissues)){
  write.xlsx(get(paste0(tissues[i],"_one2one_SNP_pig")),file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/",tissues[i],"_one2one_SNP_pig.xlsx"))
  write.xlsx(get(paste0(tissues[i],"_one2one_SNP_hum")),file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/",tissues[i],"_one2one_SNP_hum.xlsx"))
  write.xlsx(get(paste0(tissues[i],"_samevariant_pig")),file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/",tissues[i],"_samevariant_pig.xlsx"))
  write.xlsx(get(paste0(tissues[i],"_samevariant_hum")),file=paste0("/Users/baizhonghao/Downloads/human-pig GTEx/",tissues[i],"_samevariant_hum.xlsx"))
}



library(xlsx)
summary_var<-read.xlsx(file = '/Users/baizhonghao/Downloads/samevariant_summary.xlsx',sheetName = "Sheet1")
var_sum<-array(NA,dim=c(nrow(summary_var),3))
colnames(var_sum)<-c("Human","Pig",'Tissues')
var_sum<-as.data.frame(var_sum)
var_sum$Pig<-summary_var$slope_pig / summary_var$slope_se_pig * summary_var$ref_value
var_sum$Human<-summary_var$slope_hum / summary_var$slope_se_hum
var_sum$Tissues<-summary_var$tissues
cor<-cor(var_sum$Human,var_sum$Pig)
p_val<-cor.test(var_sum$Human,var_sum$Pig,method="spearman")

#tiff(file = "/Users/baizhonghao/Downloads/samevariant_slope_plot-edited.tiff",##reqiured to change
#     res = 300, width = 1300, height = 1300,compression = "lzw")
p1<-ggplot(var_sum,aes(x=Human,y=Pig))+
  geom_point(size=1, shape=15, color='black')+
  stat_smooth(method="lm")+
  theme_classic() +
  labs(x="the effect size of SNPs with\n same variants in human", 
       y="the effect size of SNPs with\n same variants in pig") + 
  theme(axis.title =element_text(size = 12),axis.text =element_text(size = 6, color = 'black'))+
  theme(axis.text.x = element_text(hjust = 1,size=5.2))+
  theme(legend.position = "none")+
  theme(axis.title.y.left = element_text(vjust = 2))+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 18))+
  theme(axis.text=element_text(color="black"))+
  theme(axis.title=element_text(face = "bold", color="black"))+
  scale_color_manual(values=colors)+
  annotate("text",x=-12,y=10,label=paste0("R=",round(cor,4),", p=",round(p_val$p.value,4),size=6))
#annotate("text",x=12.3,y=15,label=paste0(italic_R,"=",corr,", ",italic_p,"=",p.val$p.value),size=6)
#dev.off()
ggsave(p1,file='/Users/baizhonghao/Downloads/SNPef.pdf',width=5,height=5,dpi=300)

tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/samevariant_slope_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(var_sum,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The effect size(slope/se) of samevariant SNPs",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

summary_var<-read.xlsx(file = '/Users/baizhonghao/Downloads/samevariant_summary.xlsx',sheetName = "Sheet1")
var_sum_afc<-array(NA,dim=c(nrow(summary_var),2))
colnames(var_sum_afc)<-c("Human","Pig")
var_sum_afc<-as.data.frame(var_sum_afc)
var_sum_afc$Pig<-summary_var$log2_aFC_pig * summary_var$ref_value
var_sum_afc$Human<-summary_var$log2_aFC_hum
cor<-cor(var_sum_afc$Human,var_sum_afc$Pig)
p_val<-cor.test(var_sum_afc$Human,var_sum_afc$Pig,method="spearman")
tiff(file = "/Users/baizhonghao/Downloads/human-pig GTEx/samevariant_afc_plot.tiff",
     res = 300, width = 2000, height = 2000,compression = "lzw")
ggplot(var_sum_afc,aes(x=Human,y=Pig))+ geom_point(size=1,shape=15)+geom_smooth(method=lm)+ggtitle("The effect size(aFC) of samevariant SNPs",subtitle =paste0("cor=",cor," pval=",p_val$p.value))
dev.off()

