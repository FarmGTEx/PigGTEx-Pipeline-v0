### by Jinyan Teng (kingyan312@live.cn) on September 29, 2020
options(stringsAsFactors = FALSE)
library(edgeR)
library(peer)
library(preprocessCore)
library(RNOmni)
library(data.table)
library(R.utils)
library(SNPRelate)
#----------------------------------------------------------------------------
### functions
"%&%" = function(a, b) { paste0(a, b) }
# Transform rows to a standard normal distribution
inverse_normal_transform = function(x) {
    qnorm(rank(x) / (length(x)+1))
}
#----------------------------------------------------------------------------
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
### data input <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
ARGS <- commandArgs(trailingOnly = TRUE)
file_counts = ARGS[1] # Counts file. Row is gene, column is sample; rowname is gene id, colname is sample id
file_tpm = ARGS[2] # TPM file. Row is gene, column is sample; rowname is gene id, colname is sample id
tss_annot_file = ARGS[3] # TSS annotation file
vcf.fn = ARGS[4] # Input data for genotype PCA. genotype data from imputation (VCF format)
tis = ARGS[5] # Prefix of output file, like tissue name

if (!file.exists(file_counts)) { stop("Can not find the file_counts") }
if (!file.exists(file_tpm)) { stop("Can not find the file_tpm") }
if (!file.exists(tss_annot_file)) { stop("Can not find the tss_annot_file") }
if (!file.exists(vcf.fn)) { stop("Can not find the vcf.fn") }

setwd("./")
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
### main program
#----------------------------------------------------------------------------
# Input data for TMM calculation
# Read counts matrix. Row is gene, column is sample; rowname is gene id, colname is sample id
Counts = fread(file_counts,data.table=F)
# TPM matrix. Row is gene, column is sample; rowname is gene id, colname is sample id
TPM = fread(file_tpm,data.table=F)

## 1. prepare TMM
samids = colnames(Counts) # sample id
expr_counts = Counts
expr = DGEList(counts=expr_counts) # counts
nsamples = length(samids) # sample number
ngenes = nrow(expr_counts) # gene number

# calculate TMM
y = calcNormFactors(expr, method="TMM")
TMM = cpm(y,normalized.lib.sizes=T)

# expression thresholds
count_threshold = 6
tpm_threshold = 0.1
sample_frac_threshold = 0.2
sample_count_threshold = 10

#keep the genes with >=0.1 tpm and >=6 read counts in >=20% samples.
expr_tpm = TPM[rownames(expr_counts),samids]
tpm_th = rowSums(expr_tpm >= tpm_threshold)
count_th = rowSums(expr_counts >= count_threshold)
ctrl1 = tpm_th >= (sample_frac_threshold * nsamples)
ctrl2 = count_th >= (sample_frac_threshold * nsamples)
mask = ctrl1 & ctrl2
TMM_pass = TMM[mask,] ##row is gene; column is sample

###expression values (TMM) were inverse normal transformed across samples.
TMM_inv = t(apply(TMM_pass, MARGIN = 1, FUN = inverse_normal_transform)) #apply to each row, each row represents one gene, observed values for all the samples. scale across samples.
#----------------------------------------------------------------------------
### 2. prepare bed file
dir.create("bed",showWarnings=F)

region_annot = fread(tss_annot_file) # load gtf file
geneid = region_annot$gene_id

expr_matrix = TMM_inv[rownames(TMM_inv) %in% geneid,] # expr_matrix TMM_inv

# prepare bed file for tensorQTL
bed_annot = region_annot[region_annot$gene_id %in% rownames(expr_matrix),]
bed = data.frame(bed_annot,expr_matrix[bed_annot$gene_id,])
bed = bed[bed[,1] %in% as.character(1:18),]
bed[,1] = as.numeric(bed[,1])
bed = bed[order(bed[,1],bed[,2]),]
colnames(bed)[1] = "#Chr"

# output bed file
fwrite(bed,file = "./bed/" %&% tis %&% ".expr_tmm_inv.bed", sep = "\t")
system("bgzip ./bed/" %&% tis %&% ".expr_tmm_inv.bed")
system("tabix -p bed ./bed/" %&% tis %&% ".expr_tmm_inv.bed.gz")
#----------------------------------------------------------------------------
### 3. Peer factor estimation
dir.create("pca_peer",showWarnings=F)
rownames(bed) = bed$gene_id
bed$gene_id = NULL
bed$`#Chr` = NULL
bed$start = NULL
bed$end = NULL
#
expr_peer = bed
#
k = 60
model = PEER()
PEER_setPhenoMean(model, as.matrix(t(expr_peer))) #NULL response means no err# #N rows (samples); G columns (Genes) #have been sorted based on genotype samples
dim(PEER_getPhenoMean(model))
PEER_setNk(model, k)
PEER_setNmax_iterations(model, 1000)
PEER_getNk(model)
PEER_update(model)
factors = PEER_getX(model)

rownames(factors) = colnames(expr_peer)
colnames(factors) = paste0("peer",c(1:k))

Alpha = PEER_getAlpha(model)

alpha0 = data.frame(Peer=c(1:k),Alpha=Alpha,Relevance = 1.0 / Alpha)
residuals0 = t(PEER_getResiduals(model))
residuals0 = data.frame(GENE_ID=rownames(expr_peer),residuals0)
colnames(residuals0) = colnames(expr_peer)
covariates0 <- data.frame(SampleID = colnames(expr_peer), factors)

# output
cat("PEER_Nk" %&% k %&% ": writing results ... ")
write.table(covariates0, "./pca_peer/" %&% tis %&% ".PEER_covariates.Nk" %&% k %&% ".txt", sep = "\t", row.names = F, quote = FALSE)
write.table(alpha0, "./pca_peer/" %&% tis %&% ".PEER_alpha.Nk" %&% k %&% ".txt", sep = "\t", row.names = F, quote = FALSE)
write.table(residuals0, "./pca_peer/" %&% tis %&% ".PEER_residuals.Nk" %&% k %&% ".txt", sep = "\t", row.names = F, quote = FALSE)
#----------------------------------------------------------------------------
### 4. Genotype PCA
snpgdsVCF2GDS(vcf.fn, "./pca_peer/" %&% tis %&% ".ccm.gds", method = "biallelic.only")
genofile <- openfn.gds("./pca_peer/" %&% tis %&% ".ccm.gds")
ccm_pca <- snpgdsPCA(genofile,num.thread=23)

pca_genotype <- ccm_pca$eigenvect[, 1:30]
colnames(pca_genotype) <- paste0("pc", 1:30)
rownames(pca_genotype) <- ccm_pca$sample.id
pca_genotype0 = data.frame(SampleID=ccm_pca$sample.id,pca_genotype)
pca_var0 = data.frame(pc=1:30,eigenval=ccm_pca$eigenval[1:30],varprop=ccm_pca$varprop[1:30])

# output
fwrite(pca_genotype0, "./pca_peer/" %&% tis %&% ".PCA_eigenvect.txt", sep = "\t", row.names = F, quote = FALSE)
fwrite(pca_var0, "./pca_peer/" %&% tis %&% ".PCA_var.txt", sep = "\t", row.names = F, quote = FALSE)
#----------------------------------------------------------------------------
### 5. merge all covariates
dir.create("covar",showWarnings=F)
pca_vect = pca_genotype0
peer_vect = covariates0
n_samples = nrow(pca_vect)

if(n_samples<200){
    cov_pc = pca_vect[,2:6]
} else if(n_samples>=200){
    cov_pc = pca_vect[,2:11]
}
cov_peer = peer_vect[,2:11]
cov0 = as.matrix(t(cbind(pca_vect[,1],cov_pc,cov_peer)))
colnames(cov0) = cov0[1,]
cov0 = cov0[-1,]
cov_output = data.frame(rownames(cov0),cov0)
colnames(cov_output)[1] = ""

# output
fwrite(cov_output,"./covar/" %&% tis %&% ".covariates.txt",sep="\t",quote=F)
#----------------------------------------------------------------------------
cat("done.\n")
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
