#! /usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
suppressMessages((library(reshape2)))
suppressMessages(library(methods))
"%&%" <- function(a,b) paste(a,b, sep = "")


get_filtered_snp_annot <- function(snp_annot_file_name) {
  snp_annot <- read.table(snp_annot_file_name, header = T, stringsAsFactors = F) %>%
    filter(!((ref_vcf == 'A' & alt_vcf == 'T') |
               (ref_vcf == 'T' & alt_vcf == 'A') |
               (ref_vcf == 'C' & alt_vcf == 'G') |
               (ref_vcf == 'G' & alt_vcf == 'C')) &
             !(is.na(rsid))) %>%
    distinct(varID, .keep_all = TRUE)
  snp_annot
}


get_maf_filtered_genotype <- function(genotype_file_name,  maf, samples) {
  gt_df <- read.table(genotype_file_name, header = T, stringsAsFactors = F, row.names = 1)
  gt_df <- gt_df[,(colnames(gt_df) %in% samples )] %>% t() %>% as.data.frame()
  effect_allele_freqs <- colMeans(gt_df) / 2
  gt_df <- gt_df[,which((effect_allele_freqs >= maf) & (effect_allele_freqs <= 1 - maf))]
  gt_df
}

get_gene_annotation <- function(gene_annot_file_name, chrom, gene_types=c('protein_coding', 'pseudogene', 'lncRNA', 'splicing', 'exon')){
  gene_df <- read.table(gene_annot_file_name, header = TRUE, stringsAsFactors = FALSE) %>%
    filter((chr == chrom) & gene_type %in% gene_types)
  gene_df
}

get_gene_type <- function(gene_annot, gene) {
  filter(gene_annot, gene_id == gene)$gene_type
}

get_gene_expression <- function(gene_expression_file_name, gene_annot) {
  expr_df <- as.data.frame(t(read.table(gene_expression_file_name, header = T, stringsAsFactors = F, row.names = 1)))
  expr_df <- expr_df %>% t() %>% as.data.frame()
  expr_df <- expr_df %>% select(one_of(intersect(gene_annot$gene_id, colnames(expr_df))))
  expr_df
}

get_gene_coords <- function(gene_annot, gene) {
  row <- gene_annot[which(gene_annot$gene_id == gene),]
  c(row$start, row$end)
}

get_cis_genotype <- function(gt_df, snp_annot, coords, cis_window) {
  snp_info <- snp_annot %>% filter((pos >= (coords[1] - cis_window) & !is.na(rsid)) & (pos <= (coords[2] + cis_window)))
  if (nrow(snp_info) == 0)
    return(NA)
  #Check of the varID exist in the data
  if (TRUE %in% (snp_info$varID %in% names(gt_df))) {
  cis_gt <- gt_df %>% select(one_of(intersect(snp_info$varID, colnames(gt_df))))
  } else {
    return(NA) # the varID doesn't exist in the gt_df dataset
  }
  column_labels <- colnames(cis_gt)
  row_labels <- rownames(cis_gt)
  # Convert cis_gt to a matrix for glmnet
  cis_gt <- matrix(as.matrix(cis_gt), ncol=ncol(cis_gt)) # R is such a bad language.
  colnames(cis_gt) <- column_labels
  rownames(cis_gt) <- row_labels
  cis_gt
}

get_covariates <- function(covariate_file_name, samples) {
  cov_df <- read.table(covariate_file_name, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  cov_df <- cov_df[,(names(cov_df) %in% samples )] %>% t() %>% as.data.frame()
  cov_df
}

generate_fold_ids <- function(n_samples, n_folds=10) {
  n <- ceiling(n_samples / n_folds)
  fold_ids <- rep(1:n_folds, n)
  sample(fold_ids[1:n_samples])
}

adjust_for_covariates <- function(expression_vec, cov_df) {
  combined_df <- cbind(expression_vec, cov_df)
  expr_resid <- summary(lm(expression_vec ~ ., data=combined_df))$residuals
  expr_resid <- scale(expr_resid, center = TRUE, scale = TRUE)
  expr_resid
}

calc_R2 <- function(y, y_pred) {
  tss <- sum(y**2)
  rss <- sum((y - y_pred)**2)
  1 - rss/tss
}

calc_corr <- function(y, y_pred) {
  sum(y*y_pred) / (sqrt(sum(y**2)) * sqrt(sum(y_pred**2)))
}

nested_cv_elastic_net_perf <- function(x, y, n_samples, n_train_test_folds, n_k_folds, alpha, samples) {
  # Gets performance estimates for k-fold cross-validated elastic-net models.
  # Splits data into n_train_test_folds disjoint folds, roughly equal in size,
  # and for each fold, calculates a n_k_folds cross-validated elastic net model. Lambda parameter is
  # cross validated. Then get performance measures for how the model predicts on the hold-out
  # fold. Get the coefficient of determination, R^2, and a p-value, where the null hypothesis
  # is there is no correlation between prediction and observed.
  #
  # The mean and standard deviation of R^2 over all folds is then reported, and the p-values
  # are combined using Fisher's method.
  R2_folds <- rep(0, n_train_test_folds)
  corr_folds <- rep(0, n_train_test_folds)
  zscore_folds <- rep(0, n_train_test_folds)
  pval_folds <- rep(0, n_train_test_folds)
  # Outer-loop split into training and test set.
  train_test_fold_ids <- generate_fold_ids(n_samples, n_folds=n_train_test_folds)
  for (test_fold in 1:n_train_test_folds) {
    train_idxs <- which(train_test_fold_ids != test_fold)
    test_idxs <- which(train_test_fold_ids == test_fold)
    x_train <- x[(rownames(x) %in% samples[train_idxs]), ]
    y_train <- y[(rownames(y) %in% rownames(x_train))]
    x_test <- x[(rownames(x) %in% samples[test_idxs]), ]
    y_test <- y[(rownames(y) %in% rownames(x_test))]
    # Inner-loop - split up training set for cross-validation to choose lambda.
    cv_fold_ids <- generate_fold_ids(length(y_train), n_k_folds)
    y_pred <- tryCatch({
      # Fit model with training data.
      fit <- cv.glmnet(x_train, y_train, nfolds = n_k_folds, alpha = alpha, type.measure='mse', foldid = cv_fold_ids)
      # Predict test data using model that had minimal mean-squared error in cross validation.
      predict(fit, x_test, s = 'lambda.min')},
      # if the elastic-net model did not converge, predict the mean of the y_train (same as all non-intercept coef=0)
      error = function(cond) rep(mean(y_train), length(y_test)))
    R2_folds[test_fold] <- calc_R2(y_test, y_pred)
    # Get p-value for correlation test between predicted y and actual y.
    # If there was no model, y_pred will have var=0, so cor.test will yield NA.
    # In that case, give a random number from uniform distribution, which is what would
    # usually happen under the null.
    corr_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor(y_pred, y_test), 0)
    zscore_folds[test_fold] <- atanh(corr_folds[test_fold])*sqrt(length(y_test) - 3) # Fisher transformation
    pval_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor.test(y_pred, y_test)$p.value, runif(1))
  }
  R2_avg <- mean(R2_folds)
  R2_sd <- sd(R2_folds)
  rho_avg <- mean(corr_folds)
  rho_se <- sd(corr_folds)
  rho_avg_squared <- rho_avg**2
  # Stouffer's method for combining z scores.
  zscore_est <- sum(zscore_folds) / sqrt(n_train_test_folds)
  zscore_pval <- 2*pnorm(abs(zscore_est), lower.tail = FALSE)
  # Fisher's method for combining p-values: https://en.wikipedia.org/wiki/Fisher%27s_method
  pval_est <- pchisq(-2 * sum(log(pval_folds)), 2*n_train_test_folds, lower.tail = F)
  list(R2_avg=R2_avg, R2_sd=R2_sd, pval_est=pval_est, rho_avg=rho_avg, rho_se=rho_se, rho_zscore=zscore_est, rho_avg_squared=rho_avg_squared, zscore_pval=zscore_pval)
}

do_covariance <- function(gene_id, cis_gt, rsids, varIDs) {
  model_gt <- cis_gt[,varIDs, drop=FALSE]
  colnames(model_gt) <- rsids
  geno_cov <- cov(model_gt)
  geno_cov[lower.tri(geno_cov)] <- NA
  cov_df <- melt(geno_cov, varnames = c("rsid1", "rsid2"), na.rm = TRUE) %>%
              mutate(gene=gene_id) %>%
              select(GENE=gene, RSID1=rsid1, RSID2=rsid2, VALUE=value) %>%
              arrange(GENE, RSID1, RSID2)
  cov_df
}

main <- function(snp_annot_file, gene_annot, genotype_file, expr_df,
                 covariates_file, chrom, prefix, tissue, type, i, maf=0.01, n_folds=10, n_train_test_folds=5,
                 seed=20211114, cis_window=1e6, alpha=0.5, null_testing=FALSE) {
  #gene_annot <- get_gene_annotation(gene_annot_file, chrom)
  #expr_df <- get_gene_expression(expression_file, gene_annot)
  samples <- rownames(expr_df)
  n_samples <- length(samples)
  genes <- colnames(expr_df)
  n_genes <- length(expr_df)
  snp_annot <- get_filtered_snp_annot(snp_annot_file)
  gt_df <- get_maf_filtered_genotype(genotype_file, maf, samples)
  covariates_df <- get_covariates(covariates_file, samples)
  
  # Set seed----
  seed <- ifelse(is.na(seed), sample(1:1000000, 1), seed)
  set.seed(seed)
  
  # Prepare output data----
  model_summary_file <- '../../type/' %&% type %&% '/summary/' %&% tissue %&% "/" %&% chrom %&% "/" %&% type %&% '.' %&% tissue %&% '_' %&% prefix %&% '_chr' %&% chrom %&% '.' %&% genes[i] %&% '.model_summaries.txt'
  #model_summary_cols <- c('gene_id', 'gene_name', 'gene_type', 'alpha', 'n_snps_in_window', 'n_snps_in_model', 'lambda_min_mse',
                          #'test_R2_avg', 'test_R2_sd', 'cv_R2_avg', 'cv_R2_sd', 'in_sample_R2',
                          #'nested_cv_fisher_pval', 'rho_avg', 'rho_se', 'rho_zscore', 'rho_avg_squared', 'zscore_pval',
                          #'cv_rho_avg', 'cv_rho_se', 'cv_rho_avg_squared', 'cv_zscore_est', 'cv_zscore_pval', 'cv_pval_est')
  #write(model_summary_cols, file = model_summary_file, ncol = 24, sep = '\t')
  
  weights_file <- '../../type/' %&% type %&% '/weights/' %&% tissue %&% "/" %&% chrom %&% "/" %&% type %&% '.' %&% tissue %&% '_' %&% prefix %&% '_chr' %&% chrom %&% '.' %&% genes[i] %&% '.weights.txt'
  #weights_col <- c('gene_id', 'rsid', 'varID', 'ref', 'alt', 'beta')
  #write(weights_col, file = weights_file, ncol = 6, sep = '\t')
  
  tiss_chr_summ_f <- '../../type/' %&% type %&% '/summary/' %&% tissue %&% "/" %&% chrom %&% "/" %&% type %&% '.' %&% tissue %&% '_' %&% prefix %&% '_chr' %&% chrom %&% '_summary.txt'
  tiss_chr_summ_col <- c('n_samples', 'chrom', 'cv_seed', 'n_genes')
  tiss_chr_summ <- data.frame(n_samples, chrom, seed, n_genes)
  colnames(tiss_chr_summ) <- tiss_chr_summ_col
  write.table(tiss_chr_summ, file = tiss_chr_summ_f, quote = FALSE, row.names = FALSE, sep = '\t')
  
  covariance_file <- '../../type/' %&% type %&% '/covariances/' %&% tissue %&% "/" %&% chrom %&% "/" %&% type %&% '.' %&% tissue %&% '_' %&% prefix %&% '_chr' %&% chrom %&% '.' %&% genes[i] %&% '.covariances.txt'
  #covariance_col <- c('GENE', 'RSID1', 'RSID2', 'VALUE')
  #write(covariance_col, file = covariance_file, ncol = 4, sep = ' ')
  
  # Attempt to build model for each gene----
 
    cat(i, "/", n_genes, "\n")
    gene <- genes[i]
    gene_name <- gene_annot$gene_name[gene_annot$gene_id == gene]
    gene_type <- get_gene_type(gene_annot, gene)
    coords <- get_gene_coords(gene_annot, gene)
    cis_gt <- get_cis_genotype(gt_df, snp_annot, coords, cis_window)
    if (all(is.na(cis_gt))) {
      # No snps within window for gene.
      model_summary <- c(gene, gene_name, gene_type, alpha, 0, 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      write(model_summary, file = model_summary_file, ncol = 24, sep = '\t')
      #next
    } else {
      model_summary <- c(gene, gene_name, gene_type, alpha, ncol(cis_gt), 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      if (ncol(cis_gt) >= 2) {
        expression_vec <- expr_df[,i]
        adj_expression <- adjust_for_covariates(expression_vec, covariates_df)
        adj_expression <- as.matrix(adj_expression[(rownames(adj_expression) %in% rownames(cis_gt)),])
        if (null_testing)
          adj_expression <- sample(adj_expression)
        perf_measures <- nested_cv_elastic_net_perf(cis_gt, adj_expression, n_samples, n_train_test_folds, n_folds, alpha, samples)
        R2_avg <- perf_measures$R2_avg
        R2_sd <- perf_measures$R2_sd
        pval_est <- perf_measures$pval_est
        rho_avg <- perf_measures$rho_avg
        rho_se <- perf_measures$rho_se
        rho_zscore <- perf_measures$rho_zscore
        rho_avg_squared <- perf_measures$rho_avg_squared
        zscore_pval <- perf_measures$zscore_pval
        # Fit on all data
        cv_fold_ids <- generate_fold_ids(length(adj_expression), n_folds)
        fit <- tryCatch(cv.glmnet(cis_gt, adj_expression, nfolds = n_folds, alpha = 0.5, type.measure='mse', foldid = cv_fold_ids, keep = TRUE),
                        error = function(cond) {message('Error'); message(geterrmessage()); list()})
        if (length(fit) > 0) {
          cv_R2_folds <- rep(0, n_folds)
          cv_corr_folds <- rep(0, n_folds)
          cv_zscore_folds <- rep(0, n_folds)
          cv_pval_folds <- rep(0, n_folds)
          best_lam_ind <- which.min(fit$cvm)
          for (j in 1:n_folds) {
            fold_idxs <- which(cv_fold_ids == j)
            adj_expr_fold_pred <- fit$fit.preval[fold_idxs, best_lam_ind]
            cv_R2_folds[j] <- calc_R2(adj_expression[fold_idxs], adj_expr_fold_pred)
            cv_corr_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor(adj_expr_fold_pred, adj_expression[fold_idxs]), 0)
            cv_zscore_folds[j] <- atanh(cv_corr_folds[j])*sqrt(length(adj_expression[fold_idxs]) - 3) # Fisher transformation
            cv_pval_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor.test(adj_expr_fold_pred, adj_expression[fold_idxs])$p.value, runif(1))
          }
          cv_R2_avg <- mean(cv_R2_folds)
          cv_R2_sd <- sd(cv_R2_folds)
          adj_expr_pred <- predict(fit, as.matrix(cis_gt), s = 'lambda.min')
          training_R2 <- calc_R2(adj_expression, adj_expr_pred)
          
          cv_rho_avg <- mean(cv_corr_folds)
          cv_rho_se <- sd(cv_corr_folds)
          cv_rho_avg_squared <- cv_rho_avg**2
          # Stouffer's method for combining z scores.
          cv_zscore_est <- sum(cv_zscore_folds) / sqrt(n_folds)
          cv_zscore_pval <- 2*pnorm(abs(cv_zscore_est), lower.tail = FALSE)
          cv_pval_est <- pchisq(-2 * sum(log(cv_pval_folds)), 2*n_folds, lower.tail = F)
          
          if (fit$nzero[best_lam_ind] > 0) {
            
            weights <- fit$glmnet.fit$beta[which(fit$glmnet.fit$beta[,best_lam_ind] != 0), best_lam_ind]
            weighted_snps <- names(fit$glmnet.fit$beta[,best_lam_ind])[which(fit$glmnet.fit$beta[,best_lam_ind] != 0)]
            weighted_snps_info <- snp_annot %>% filter(varID %in% weighted_snps) %>% select(rsid, varID, ref_vcf, alt_vcf)
            weighted_snps_info$gene <- gene
            weighted_snps_info <- weighted_snps_info %>%
              merge(data.frame(weights = weights, varID=weighted_snps), by = 'varID') %>%
              select(gene, rsid, varID, ref_vcf, alt_vcf, weights)
            write.table(weighted_snps_info, file = weights_file, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')
            covariance_df <- do_covariance(gene, cis_gt, weighted_snps_info$rsid, weighted_snps_info$varID)
            write.table(covariance_df, file = covariance_file, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")
            model_summary <- c(gene, gene_name, gene_type, alpha, ncol(cis_gt), fit$nzero[best_lam_ind], fit$lambda[best_lam_ind], R2_avg, R2_sd, cv_R2_avg, cv_R2_sd, training_R2, pval_est,
                               rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval, cv_rho_avg, cv_rho_se, cv_rho_avg_squared, cv_zscore_est, cv_zscore_pval, cv_pval_est)
          } else {
            model_summary <- c(gene, gene_name, gene_type, alpha, ncol(cis_gt), 0, fit$lambda[best_lam_ind], R2_avg, R2_sd,
                               cv_R2_avg, cv_R2_sd, training_R2, pval_est, rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval,
                               cv_rho_avg, cv_rho_se, cv_rho_avg_squared, cv_zscore_est, cv_zscore_pval, cv_pval_est)
          }
        } else {
          model_summary <- c(gene, gene_name, gene_type, alpha, ncol(cis_gt), 0, NA, R2_avg, R2_sd, NA, NA, NA, pval_est, rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval,
                             NA, NA, NA, NA, NA, NA)
        }
      }
      write(model_summary, file = model_summary_file, ncol = 24, sep = '\t')
    }
	gc()
}

