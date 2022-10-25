library(data.table)
"%&%" = function(a, b) paste0(a, b)

ARGS <- commandArgs(trailingOnly = TRUE)
wkdir = ARGS[1]
ieqtl_dir = ARGS[2]
fdr_cutoff = ARGS[3] # 0.05, 0.01

setwd(wkdir)

interQTLs_ieGene_celltype = list()
tiss = unique(gsub("[.](.)*", "", list.files(ieqtl_dir)))
for (tis in tiss) {
    #tis = tiss[3]
    message(tis)

    profile = list.files(ieqtl_dir %&% tis, pattern = "gz")
    pronames = gsub("(.)*interactions_pro", "pro", gsub(".cis_qtl_top_assoc.txt.gz", "", profile))

    ieGenes = list()
    for (pro in pronames) {
        #pro = pronames[1]
        message(pro)
        file_pro = ieqtl_dir %&% "/" %&% tis %&% "/" %&% tis %&% ".interactions_" %&% pro %&% ".cis_qtl_top_assoc.txt.gz"
        top_assoc_pro = fread(file_pro)
        n1 = sum(top_assoc_pro$pval_adj_bh < 0.01)
        top_assoc_pro$is_ieGene = top_assoc_pro$pval_adj_bh < fdr_cutoff &
            p.adjust(top_assoc_pro$pval_g, method = "BH") < fdr_cutoff &
            p.adjust(top_assoc_pro$pval_i, method = "BH") < fdr_cutoff
        ieGenes[[pro]] = top_assoc_pro
        print(sum(top_assoc_pro$is_ieGene))
    }
    interQTLs_ieGene_celltype[[tis]] = ieGenes
}

saveRDS(interQTLs_ieGene_celltype, "ieGene.fdr" %&% fdr_cutoff %&% ".RDS")

# results summary
sumDF_fdr = data.frame()
for (tis in names(interQTLs_ieGene_celltype)) {
    tis_tmp = interQTLs_ieGene_celltype[[tis]]
    for (pro in names(tis_tmp)) {
        sumDF_fdr = rbind(sumDF_fdr, c(tis, pro, sum(tis_tmp[[pro]]$is_ieGene)))
    }
}
colnames(sumDF_fdr) = c("Tissue", "Factor", "ieGenes")
fwrite(sumDF_fdr, "sumDF_fdr" %&% fdr_cutoff %&% ".csv")






