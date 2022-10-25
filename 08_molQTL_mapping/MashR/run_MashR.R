#!/BIGDATA2/scau_hzhang_1/software/microsoft-r-open/microsoft-r-open-MRO4.0.2-X11/bin/Rscript
options(stringsAsFactors = FALSE)
library(data.table)
library(R.utils)
library(mashr)
library(doParallel)
library(foreach)

"%&%" = function(a, b) { paste0(a, b) }
readFile = function(f){
    if(grepl(".RDS", f, ignore.case=T)){
        DF = readRDS(f)
    } else{
        DF = fread(f)
    }
    return(DF)
}
#----------------------------------------------------------------------------
ARGS <- commandArgs(trailingOnly = TRUE)
file_strong = ARGS[1]
file_random = ARGS[2]
dropna = as.logical(as.numeric(ARGS[3]))
dir_output = ARGS[4]
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
if (length(dir_output) == 0){
    dir_output = "./"
} else {
    dir_output = dir_output %&% "/"
}
dir.create(dir_output)
set.seed(9823)
############ strong pairs
zval_data.strong = as.data.frame(readFile(file_strong))

rownames(zval_data.strong) = zval_data.strong$pair_id
zval_data.strong$pair_id = NULL
if (dropna){
    zval_data.strong = zval_data.strong[rowSums(is.na(zval_data.strong))==0,]
} else {
    zval_data.strong[is.na(zval_data.strong)] = 0
}
zval_data.strong = as.matrix(zval_data.strong)

############ all pairs of a random subset of 1000000 tests
zval_data.random.subset = as.data.frame(readFile(file_random))
rownames(zval_data.random.subset) = zval_data.random.subset$pair_id
zval_data.random.subset$pair_id = NULL
zval_data.random.subset = as.matrix(zval_data.random.subset)

if (sum(colnames(zval_data.random.subset) != colnames(zval_data.strong))>0){
    stop("Different colnames between random subset and strong subset.")
}

## Correlation structure
data.temp = mash_set_data(Bhat=zval_data.random.subset,alpha=1,zero_Bhat_Shat_reset=1e6)
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)

data.random = mash_set_data(Bhat=zval_data.random.subset,alpha=1,V=Vhat,zero_Bhat_Shat_reset=1e6)
data.strong = mash_set_data(Bhat=zval_data.strong,alpha=1,V=Vhat,zero_Bhat_Shat_reset=1e6)

## Data driven covariances
if (file.exists(dir_output %&% "U.ed.RDS")){
    U.ed = readRDS(dir_output %&% "U.ed.RDS")
} else {
    U.pca = cov_pca(data.strong,5)
    U.ed = cov_ed(data.strong, U.pca)
    saveRDS(U.ed, dir_output %&% "U.ed.RDS")
}

## Fit mash model (estimate mixture proportions)
U.c = cov_canonical(data.random)
m.r = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)
saveRDS(m.r, dir_output %&% "m.r.RDS")

## Compute posterior summaries
m.s = mash(data.strong, g=get_fitted_g(m.r), fixg=TRUE)

saveRDS(m.s, dir_output %&% "m.s_zval.RDS")
saveRDS(get_lfsr(m.s), dir_output %&% "lfsr_m.s.RDS")
saveRDS(get_pm(m.s), dir_output %&% "pm_m.s.RDS")
saveRDS(get_psd(m.s), dir_output %&% "psd_m.s.RDS")
saveRDS(get_significant_results(m.s), dir_output %&% "significant_results_m.s.RDS")
saveRDS(get_pairwise_sharing(m.s), dir_output %&% "pairwise_sharing_m.s.RDS")
saveRDS(get_pairwise_sharing(m.s, factor=0), dir_output %&% "pairwise_sharing_factor0_m.s.RDS")
saveRDS(get_loglik(m.s), dir_output %&% "loglik_m.s.RDS")
saveRDS(get_estimated_pi(m.s), dir_output %&% "estimated_pi_m.s.RDS")

