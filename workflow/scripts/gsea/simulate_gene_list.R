library(glue)
library(dplyr)

Xid = snakemake@wildcards[['x_id']]
SIM_ID = snakemake@wildcards[['sim_id']]
print(glue('X = {Xid}, sim_id = {SIM_ID}'))

manifest <- readr::read_delim('config/gsea_sim/gsea_sim_manifest.tsv')
X <- readRDS(glue('results/gsea/{Xid}/X.rds')) # called XX to avoid conflict with column name `X`

sigmoid <- function(x){1 / (1 + exp(-x))}

simulate <- function(beta0, beta, q, seed){
    set.seed(seed)
    n <- nrow(X)
    p <- ncol(X)

    idx <- sort(sample(p, q, replace=F))
    logit <- drop(X[, idx, drop = F] %*% rep(beta, q)) + beta0
    y <- rbinom(length(logit), 1, sigmoid(logit))
    sim <- list(idx = idx, logit = logit, y = y, n=n, p=p, X = Xid, beta0=beta, beta=beta, q=q)
    return(sim)
}

sim <- manifest %>% 
    dplyr::filter(X == Xid, sim_id == SIM_ID) %>%
    rowwise() %>%
    mutate(sim = list(simulate(beta0, beta, q, seed)))

saveRDS(sim$sim[[1]], snakemake@output$rds)
