library(parallel)
suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
suppressMessages((library(reshape2)))
suppressMessages(library(methods))

"%&%" <- function(a,b) paste(a,b, sep='')

argv <- commandArgs(trailingOnly = TRUE)
chrom <- argv[1]
tissue <- argv[2]
type <- argv[3]
work_path <- argv[4] 
software_path <- argv[5]
setwd(work_path)
gene_annot_file <- "./type/" %&% type %&% "/tissue_annotation/" %&% type %&% "_annot.txt"
expression_file <- "./type/" %&% type %&% "/tissue_expression/run_expression/" %&% tissue %&% "." %&% type %&% ".transformed_expression.txt"
covariates_file <- "./type/" %&% type %&% "/tissue_covar/" %&% tissue %&% ".covariates.txt"
snp_annot_file <- "./tissue_snp_annotation/" %&% tissue %&% "/" %&% tissue %&% ".snp_annot.chr" %&% chrom %&% ".txt"
genotype_file <- "./tissue_gene_matrix/" %&% tissue %&% "/" %&% tissue %&% ".genotype.chr" %&% chrom %&% ".txt"
prefix <- "Model_training"

get_gene_annotation <- function(gene_annot_file_name, chrom, gene_types=c('protein_coding', 'pseudogene', 'lncRNA', 'splicing', 'exon')){
  gene_df <- read.table(gene_annot_file_name, header = TRUE, stringsAsFactors = FALSE) %>%
    filter((chr == chrom) & gene_type %in% gene_types)
  gene_df
}
get_gene_expression <- function(gene_expression_file_name, gene_annot) {
  expr_df <- as.data.frame(t(read.table(gene_expression_file_name, header = T, stringsAsFactors = F, row.names = 1)))
  expr_df <- expr_df %>% t() %>% as.data.frame()
  expr_df <- expr_df %>% select(one_of(intersect(gene_annot$gene_id, colnames(expr_df))))
  expr_df
}
####
gene_annot <- get_gene_annotation(gene_annot_file, chrom)
expr_df <- get_gene_expression(expression_file, gene_annot)
n_gene <- length(expr_df)


parallel_main = function(i) {
	source(paste0(software_path, "/code/parallel_nested_cv_elnet.R")
	main(snp_annot_file, gene_annot, genotype_file, expr_df, covariates_file, as.numeric(chrom), prefix, tissue, type, as.numeric(i),null_testing=FALSE)
}


cl = makeCluster(11, type = "FORK")
clusterEvalQ(cl, c("dplyr", "glmnet", "reshape2", "methods"))
clusterExport(cl, c("snp_annot_file", "gene_annot", "genotype_file", "expr_df", "covariates_file", "prefix", "tissue", "type", "chrom"))
parLapply(cl, 1:n_gene, parallel_main)
stopCluster(cl)

