library(data.table)
#----------------------------------------------------------------------------
"%&%" = function(a, b) { paste0(a, b) }
###Transform rows to a standard normal distribution
inverse_normal_transform = function(x) {
    # qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))) # cattle-gtex
    qnorm(rank(x) / (length(x)+1)) # Science 2020
}
#----------------------------------------------------------------------------
ARGS <- commandArgs(trailingOnly = TRUE)
dir_deconvolution_results = ARGS[1]
gsub_pattern=ARGS[2] # to get tissue names
dir_covariates=ARGS[3]
medianPro_cutoff=ARGS[4] # default: 0.1
dir_output = ARGS[5]

fs = list.files(dir_deconvolution_results)
tis_names = gsub(gsub_pattern,"",fs)

dir.create(dir_output)
for (i in 1:length(tis_names)) {
    # i=11
    tis = tis_names[i]
    covar = as.data.frame(fread(dir_covariates %&% "/" %&% tis %&% ".covariates.txt"))
    n_samples = ncol(covar) - 1
    samids = colnames(covar)[2:ncol(covar)]

    proCelltypes = fread(dir_deconvolution_results %&% "/" %&% fs[grep(tis,fs)],data.table=F)
    if (sum(is.na(match(samids,proCelltypes$Mixture)))>0){
        stop(tis)
        next;
    }
    proCelltypes = proCelltypes[match(samids,proCelltypes$Mixture),]
    colnames(proCelltypes) = gsub(" ","_",colnames(proCelltypes))
    colnames(proCelltypes) = gsub("/",".",colnames(proCelltypes))
    celltypes = colnames(proCelltypes)[2:(ncol(proCelltypes)-3)]
    cell_types_pass = celltypes[apply(proCelltypes[, celltypes], 2, median) > medianPro_cutoff]

    for (ct in cell_types_pass) {
        pro_cell = inverse_normal_transform(proCelltypes[[ct]])
        inter_cell = data.frame(samids,pro_cell)
        fwrite(inter_cell,"./"%&%dir_output%&%"/" %&% tis %&% ".interactions_pro_"%&%ct%&%".txt",sep="\t",quote=F,col.names=F)
    }
}

cat(unique(gsub("[.](.)*","",list.files("./"%&%dir_output))))


