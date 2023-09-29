#!
argv <- commandArgs(trailingOnly = TRUE)
workdir <- argv[1]
outdir <- argv[2]
setwd(workdir)
#############========================================================================#############
### Extract relevant annotation information of eQTL
#############========================================================================#############
### Read the corresponding file
gene_annotation = data.table::fread("Sus_scrofa.Sscrofa11.1.100.chr.gtf")
gene_annotation = data.frame(gene_annotation)
### Split the last column and extract the corresponding data column
result = data.frame(chr = gene_annotation[, 1], gene_id = 1, gene_name = 1, start = gene_annotation[, 4], end = gene_annotation[, 5], gene_type = 1)
result[, 2] = apply(gene_annotation, 1, function(x) {
  d = unlist(strsplit(x[9], ";"))
  if (length(grep("gene_id", d)) > 0) {
    ttt = unlist(strsplit(d[grep("gene_id", d)], " "))[2]
  } else {
    ttt = NA	#Missing information is replaced by Na
  }
  return(ttt)
})
result[, 3] = apply(gene_annotation, 1, function(x) {
  d = unlist(strsplit(x[9], ";"))
  if (length(grep("gene_name", d)) > 0) {
    ttt = unlist(strsplit(d[grep("gene_name", d)], " "))[3]
  } else {
    ttt = NA
  }
  return(ttt)
})
result[, 6] = apply(gene_annotation, 1, function(x) {
  d = unlist(strsplit(x[9], ";"))
  if (length(grep("gene_biotype", d)) > 0) {
    ttt = unlist(strsplit(d[grep("gene_biotype", d)], " "))[3]
  } else {
    ttt = NA
  }
  return(ttt)
})
### Extract the information with 'gene' in the third column of annotation information
new_gene_annotation = result[which(gene_annotation[, 3] == "gene"), ]
for (i in 1:nrow(new_gene_annotation)) {
  new_gene_annotation[i, 2] = gsub("\"", "", new_gene_annotation[i, 2])# Replace " with empty
  new_gene_annotation[i, 3] = gsub("\"", "", new_gene_annotation[i, 3])
  new_gene_annotation[i, 6] = gsub("\"", "", new_gene_annotation[i, 6])
}
write.table(new_gene_annotation, paste0(outdir, "eQTL_annot.txt"), row.names = F, quote = F,sep = "\t")



#############========================================================================#############
### Generate annotation files of exon, splicing and lncrna
#############========================================================================#############
### exon
exon = read.table("Sus_scrofa.Sscrofa11.1.100.exon.glist")
chr = exon$V1
start = exon$V2
end = exon$V3
gene_id = exon$V4
gene_id = gsub(":", "_", gene_id)
gene_name = rep(NA, length(gene_id))
gene_type = rep("exon", length(gene_id))
data = cbind(chr, gene_id, gene_name, start, end, gene_type)
write.table(data, paste0(outdir, "exon_annot.txt"), row.names = F, quote = F)
### splicing
splicing = read.table("Sus_scrofa.Sscrofa11.1.100.splicing.glist")
chr = splicing$V1
start = splicing$V2
end = splicing$V3
gene_id = splicing$V4
gene_id = paste("Splicing_", gene_id, sep = "")
gene_id = gsub(":", "_", gene_id)
gene_name = rep(NA, length(gene_id))
gene_type = rep("splicing", length(gene_id))
data = cbind(chr, gene_id, gene_name, start, end, gene_type)
write.table(data, paste0(outdir, "splicing_annot.txt"), row.names = F, quote = F)
### lncRNA
lncRNA = read.table("Sus_scrofa.Sscrofa11.1.100.lncRNA.glist")
chr = lncRNA$V1
start = lncRNA$V2
end = lncRNA$V3
gene_id = lncRNA$V4
gene_id = gsub("&", "_", gene_id)
# gene_id = paste("Splicing_", gene_id, sep = "")
# gene_id = gsub(":", "_", gene_id)
gene_name = rep(NA, length(gene_id))
gene_type = rep("lncRNA", length(gene_id))
data = cbind(chr, gene_id, gene_name, start, end, gene_type)
write.table(data, paste0(outdir, "lncRNA_annot.txt"), row.names = F, quote = F)
