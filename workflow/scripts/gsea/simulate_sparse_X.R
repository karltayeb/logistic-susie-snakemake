sim_X_sparse <- function(seed=0, ...){
    set.seed(seed)
    as.matrix(logisticsusie:::sim_X_sparse(...))
}

params <- snakemake@params$params
# params = list(n=500, p=100, pi1=0.1, p_flip=0.02, seed=1)
X <- rlang::exec(sim_X_sparse, !!!params)
saveRDS(X, snakemake@output[[1]])
saveRDS(cor(X), snakemake@output[[2]])