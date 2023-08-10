# fit lasso to lasso simso
library(tictoc)
library(glmnet)

run <- function(factor, output){
    data <- readRDS('results/pdac/pdac_data.rds')
    sims <- readRDS('results/pdac/pdac_lasso_sims.rds')
    X <- readRDS('results/pdac/pdac_msigdb_X.rds')

    beta <- sims[[factor]]$beta
    idx <- names(beta[beta!=0])

    # fit lasso to each simulation
    tic()
    message('fitting lasso')
    lasso_fits <- purrr::map(sims[[factor]]$y, ~cv.glmnet(X, .x, family='binomial'))
    toc()

    # count number of gene sets recovered
    refit_beta <- purrr::map_dfr(lasso_fits, ~tail(coef(.x, s='lambda.1se')[,1], -1))
    selected = colSums(refit_beta != 0)  # tally which variables are selected by lasso across simulations
    recovered = selected[idx]  # which among these are TP?

    res = list(
        factor = factor,
        beta = beta,
        idx = idx,
        lasso_fits = lasso_fits,
        refit_beta = refit_beta,
        selected = selected,
        recovered = recovered
    )
    saveRDS(res, output)
}


run(snakemake@wildcards[['factor']], snakemake@output[[1]])