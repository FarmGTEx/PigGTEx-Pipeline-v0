library(data.table)
"%&%" = function(a,b) paste0(a,b)

ARGS <- commandArgs(trailingOnly = TRUE)
file_tss = ARGS[1] # Sus_scrofa.Sscrofa11.1.100.tss.gz
file_gtf = ARGS[2] # ./PreDATA/Sus_scrofa.Sscrofa11.1.100.gene.gtf
dir_frq = ARGS[3] # ~/snp_frq/frq/
dir_BESD_file = ARGS[4] # ./PreDATA/BESD_file_pcg/

annot = fread(file_tss)
annot_gtf = fread(file_gtf)

tis_names = c("Adipose", "Large_intestine", "Small_intestine", "Blood", "Fetal_thymus", "Lymph_node", "Milk", "Spleen", "Cartilage", "Synovial_membrane", "Brain", "Embryo", "Heart", "Kidney", "Liver", "Lung", "Testis", "Muscle", "Oocyte", "Artery", "Ovary", "Uterus", "Placenta", "Colon", "Duodenum", "Ileum", "Jejunum", "Macrophage", "Frontal_cortex", "Hypothalamus", "Pituitary", "Blastocyst", "Blastomere", "Morula")

for (i in 1:34){
    tis = tis_names[i]
    message(tis)
    frq = fread(dir_frq %&% tis %&% ".geno.frq.gz")
    esi = fread(dir_BESD_file %&% tis %&% ".eqtl_allpairs.esi.old")
    epi = fread(dir_BESD_file %&% tis %&% ".eqtl_allpairs.epi.old")

    matIdx = match(esi$V2,frq$SNP)
    esi$V1 = frq$CHR[matIdx]
    esi$V5 = frq$A1[matIdx]
    esi$V6 = frq$A2[matIdx]
    esi$V7 = frq$MAF[matIdx]
    esi$V4 = gsub("^[0-9]+_|_(A|T|C|G)","",esi$V2)
    fwrite(esi,dir_BESD_file %&% tis %&% ".eqtl_allpairs.esi", col.names=F, sep = "\t")

    annotIdx = match(epi$V2,annot$gene_id)
    epi$V1 = annot$`#chr`[annotIdx]
    epi$V4 = annot$end[annotIdx]
    epi$V5 = epi$V2
    epi$V6 = annot_gtf$strand[match(epi$V2,annot_gtf$gene_id)]
    fwrite(epi,dir_BESD_file %&% tis %&% ".eqtl_allpairs.epi", col.names=F, sep = "\t")
}
