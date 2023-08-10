
# fit lasso to each gene list
# simulate from lasso fit

library(dplyr)
library(glmnet)

fit_lasso <- function(X, gene_list){
  y <- as.integer(rownames(X) %in% gene_list)
  lasso_fit <- glmnet::cv.glmnet(X, y, family = "binomial")
  return(lasso_fit)
}

fit_lasso_all_gene_lists <- function(X, gene_lists){
  lasso_fits <- purrr::map(gene_lists, ~fit_lasso(X, .x))
  names(lasso_fits) <- names(gene_lists)
  return(lasso_fits)
}

simulate_from_lasso <- function(X, lasso_fit, reps){
  coef_1se <- coef(lasso_fit, s = 'lambda.1se')[,1]
  logit_1se <- (X %*% tail(coef_1se, -1))[,1] + coef_1se[1]
  ysim_1se <- purrr::map(1:reps, ~rbinom(length(logit_1se), 1, 1/(1 + exp(-logit_1se))))
  res <- list(y = ysim_1se, logit = logit_1se, intercept = coef_1se[1], beta = tail(coef_1se, -1))
  return(res)
}

simulate_from_lasso_map <- function(X, lasso_fits, reps){
  simulations <- purrr::map(lasso_fits, ~simulate_from_lasso(X, .x, reps))
  names(simulations) <- names(lasso_fits)
  return(simulations)
}

pdac_msigdb_X = readRDS('results/pdac/pdac_msigdb_X.rds')
pdac_data = readRDS('results/pdac/pdac_data.rds')

pdac_lasso_fits = fit_lasso_all_gene_lists(pdac_msigdb_X, pdac_data$gene_lists)
pdac_lasso_sims = simulate_from_lasso_map(pdac_msigdb_X, pdac_lasso_fits, 20)

saveRDS(pdac_lasso_fits, 'results/pdac/pdac_lasso_fits.rds')
saveRDS(pdac_lasso_sims, 'results/pdac/pdac_lasso_sims.rds')
