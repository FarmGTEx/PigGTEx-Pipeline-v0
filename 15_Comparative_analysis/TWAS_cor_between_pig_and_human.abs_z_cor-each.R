library(data.table)
library(doParallel)
setwd("/BIGDATA2/scau_hzhang_1/USER/tengjy/Pig_GTEx_Comparative")
"%&%" = function(a,b) paste0(a,b)
dir_data = "/BIGDATA2/scau_hzhang_1/USER/tengjy/Pig_GTEx_Comparative/data/"

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
    print("No arguments supplied.")
    ##supply default values
    par_i = 1:131
}else{
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
    print(par_i)
}


human_genes = fread(dir_data %&% "15944ortho.txt",header=F)$V1
pig_genes = fread(dir_data %&% "15944pig.txt",header=F)$V1

# PCG_TWAS_SMultiXcan = fread("/BIGDATA2/scau_hzhang_1/USER/tengjy/Pig_GTEx_GWAS_and_eQTLs/TWAS/PCG_TWAS_SMultiXcan.csv.gz")
PCG_TWAS_SPrediXcan = fread("/BIGDATA2/scau_hzhang_1/USER/tengjy/Pig_GTEx_GWAS_and_eQTLs/TWAS/PCG_TWAS_SPrediXcan.csv.gz",nThread=20)
pig_traits = unique(PCG_TWAS_SPrediXcan$Trait)
pig_tissues = unique(PCG_TWAS_SPrediXcan$Tissue)

HM_TWAS_results_UKB = fread("/BIGDATA2/scau_hzhang_1/USER/tengjy/Pig_GTEx_Comparative/data/HM_TWAS_results_UKB.csv.gz",nThread=20)
HM_TWAS_results_first_round = fread("/BIGDATA2/scau_hzhang_1/USER/tengjy/Pig_GTEx_Comparative/data/HM_TWAS_results_first_round.csv.gz",nThread=20)

HM_TWAS_results = rbind(HM_TWAS_results_first_round,HM_TWAS_results_UKB)
HM_TWAS_results$gene = gsub("[.][0-9]*","",HM_TWAS_results$gene)
human_traits = unique(HM_TWAS_results$Trait)
human_tissues = unique(HM_TWAS_results$Tissue)



# # correlation between pig and human
# registerDoParallel(20)
# cnames = unlist(lapply(pig_traits[par_i],function(x)x%&%"."%&%pig_tissues))
# hnames = unlist(lapply(human_traits,function(x)x%&%"."%&%human_tissues))
# cor_res = matrix(NA,length(cnames),length(hnames),dimnames=list(cnames,hnames))
# r_mat = cor_res
# # r_mat[1,1] = 0.0
# p_mat = cor_res
# # p_mat[1,1] = 0.0
# n_mat = cor_res
# # n_mat[1,1] = 0.0
sampleDF = fread("/BIGDATA2/scau_hzhang_1/USER/tengjy/Pig_GTEx_Comparative/example/sampleDF.csv.gz")
ht = unique(sampleDF$Human_Trait)[par_i]
sampleDF = sampleDF[Human_Trait == ht]

sampleDF$gene_overlap_human_fdr0.05_pig_0.05 = ""
sampleDF$gene_overlap_pig_fdr0.05_human_0.05 = ""
sampleDF$gene_overlap_pig_fdr0.05_human_fdr0.05 = ""
sampleDF$pig_zscore = ""
sampleDF$human_zscore = ""
sampleDF$pig_pvalue = ""
sampleDF$human_pvalue = ""
sampleDF$overlap_ortho_gene = ""
sampleDF$abs_z_corr = NA
sampleDF$abs_z_corp = NA

for(sample_i in 1:nrow(sampleDF)){
    tmpDF = sampleDF[sample_i,]
    # if(tmpDF$rank > 5){next;}
    message(sample_i)
    ttt1 = strsplit(tmpDF$Pig_TWAS,".",fixed=T)[[1]][1]
    tis1 = strsplit(tmpDF$Pig_TWAS,".",fixed=T)[[1]][2]
    ttt2 = strsplit(tmpDF$Human_TWAS,".",fixed=T)[[1]][1]
    tis2 = strsplit(tmpDF$Human_TWAS,".",fixed=T)[[1]][2]
    
    # if(!tmpDF$Human_Trait %in% c("BMI","Weight","Height")){next;}
    
    # ttt1 = pig_traits[i1]
    message(ttt1)
    pig_twas_ttt = PCG_TWAS_SPrediXcan[Trait == ttt1]
    
    # ttt2 = human_traits[i2]
    human_twas_ttt = HM_TWAS_results[Trait == ttt2]
    message(ttt2)
    
    message(tis1)
    pig_twas = pig_twas_ttt[Tissue == tis1]
    pig_twas$pvalue_BH = p.adjust(pig_twas$pvalue,method = "BH")
    pig_twas2 = pig_twas[match(pig_genes,pig_twas$gene),]
    
    message(tis2)
    human_twas = human_twas_ttt[Tissue == tis2]
    human_twas$pvalue_BH = p.adjust(human_twas$pvalue,method = "BH")
    human_twas2 = human_twas[match(human_genes,human_twas$gene),]

    idx = !is.na(pig_twas2$zscore) & !is.na(human_twas2$zscore)
    n = sum(idx)
    if(n<3){next;}
    # z_cor = cor.test(pig_twas2$zscore,human_twas2$zscore)
    abs_z_cor = cor.test(abs(pig_twas2$zscore),abs(human_twas2$zscore))
    
    sampleDF$abs_z_corr[sample_i] = abs_z_cor$estimate
    sampleDF$abs_z_corp[sample_i] = abs_z_cor$p.value
    
    g1 = human_twas2$gene_name[which(pig_twas2$pvalue < 0.05 & human_twas2$pvalue_BH < 0.05)]
    g1 = g1[g1!=""]
    g2 = human_twas2$gene_name[which(human_twas2$pvalue < 0.05 & pig_twas2$pvalue_BH < 0.05)]
    g2 = g2[g2!=""]
    g3 = human_twas2$gene_name[which(human_twas2$pvalue_BH < 0.05 & pig_twas2$pvalue_BH < 0.05)]
    g3 = g3[g3!=""]
    
    sampleDF$gene_overlap_human_fdr0.05_pig_0.05[sample_i] = paste0(g1,collapse =";")
    sampleDF$gene_overlap_pig_fdr0.05_human_0.05[sample_i] = paste0(g2,collapse =";")
    sampleDF$gene_overlap_pig_fdr0.05_human_fdr0.05[sample_i] = paste0(g3,collapse =";")
    
    if(length(g1)==0 && length(g2)==0 && length(g3)==0) {next;}
    
    sampleDF$overlap_ortho_gene[sample_i] = paste0(human_twas2$gene_name[idx],collapse=";")
    sampleDF$pig_zscore[sample_i] = paste0(pig_twas2$zscore[idx],collapse=";")
    sampleDF$human_zscore[sample_i] = paste0(human_twas2$zscore[idx],collapse=";")
    sampleDF$pig_pvalue[sample_i] = paste0(pig_twas2$pvalue[idx],collapse=";")
    sampleDF$human_pvalue[sample_i] = paste0(human_twas2$pvalue[idx],collapse=";")
    # print(z_cor$estimate)
    # print(z_cor$p.value)
    # cor.test(abs(human_twas2$zscore),abs(pig_twas2$zscore))
    # cor.test(-log10(human_twas2$pvalue),-log10(pig_twas2$pvalue))
    # cor_res = rbind(cor_res,data.frame(study1 = ttt1%&%"."%&%tis1, study2 = ttt2%&%"."%&%tis2, r = z_cor$estimate, pvalue = z_cor$p.value, n = n))
    # data.frame(study1 = ttt1%&%"."%&%tis1, study2 = ttt2%&%"."%&%tis2, r = z_cor$estimate, pvalue = z_cor$p.value, n = n)
    # idx1 = which(cnames == ttt1%&%"."%&%tis1)
    # idx2 = which(cnames == ttt2%&%"."%&%tis2)
    # r_mat[ttt1%&%"."%&%tis1,ttt2%&%"."%&%tis2] = z_cor$estimate
    # p_mat[ttt1%&%"."%&%tis1,ttt2%&%"."%&%tis2] = z_cor$p.value
    # n_mat[ttt1%&%"."%&%tis1,ttt2%&%"."%&%tis2] = n
    # n_mat[idx1,idx2] = n
}
fwrite(sampleDF,"/BIGDATA2/scau_hzhang_1/USER/tengjy/Pig_GTEx_Comparative/TWAS_cor/Pig_and_Human.matched_tissue/TWAS_cor.matched_tissue.abs_z_cor."%&%par_i%&%".csv.gz")

### combine all results
if(F){
    sampleDF_all = NULL
    for(par_i in 1:131){
        message(par_i)
        sampleDF = fread("/BIGDATA2/scau_hzhang_1/USER/tengjy/Pig_GTEx_Comparative/TWAS_cor/Pig_and_Human.matched_tissue/TWAS_cor.matched_tissue.abs_z_cor."%&%par_i%&%".csv.gz",data.table = F)
        sampleDF$rank_abs_z = NA
        tmp = lapply(unique(sampleDF$Pig_Tissue), function(x) {
            tmprank = rank(sampleDF[sampleDF$Pig_Tissue == x,]$abs_z_corp)
            sampleDF[sampleDF$Pig_Tissue == x,]$rank_abs_z <<- tmprank
        })
        sampleDF[sampleDF$rank_abs_z > 5,c("gene_overlap_human_fdr0.05_pig_0.05","gene_overlap_pig_fdr0.05_human_0.05","gene_overlap_pig_fdr0.05_human_fdr0.05")] = ""
        sampleDF_all = rbind(sampleDF_all,sampleDF)
    }
    sampleDF_all$pig_zscore[sampleDF_all$rank_abs_z > 5] = ""
    sampleDF_all$human_zscore[sampleDF_all$rank_abs_z > 5] = ""
    sampleDF_all$pig_pvalue[sampleDF_all$rank_abs_z > 5] = ""
    sampleDF_all$human_pvalue[sampleDF_all$rank_abs_z > 5] = ""
    sampleDF_all$overlap_ortho_gene[sampleDF_all$rank_abs_z > 5] = ""
    fwrite(sampleDF_all,"/BIGDATA2/scau_hzhang_1/USER/tengjy/Pig_GTEx_Comparative/TWAS_cor/Pig_and_Human.matched_tissue/TWAS_cor.matched_tissue.abs_z_cor.all.csv.gz")
}


