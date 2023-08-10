library(susieR)
library(reticulate)
library(Matrix)
library(glue)
library(tictoc)

region <- snakemake@wildcards[['region']]
sim_id <- snakemake@wildcards[['region']]

# 1. load ser fit, it already has the effect sizes and standard errors
print('loading SER fit')
ser_path <- snakemake@input[['ser_path']]

np <- import("numpy")
ser <- py_load_object(ser_path)
betahat <- np$array(ser$betahat)
shat2 <- np$array(ser$shat2)
zhat <- as.vector(betahat / sqrt(shat2))

# 2. load ld
print('loading LD matrix')
ld_path <- snakemake@input[['ld_path']]
R <- readRDS(ld_path)


# 3. fit susie_rss
n = snakemake@params[['n']]
estimate_prior_variance = (snakemake@params[['estimate_prior_variance']] == 'true')
print(glue('fitting SuSiE RSS: n = {n}, estimate_prior_variance = {estimate_prior_variance}'))

tic('fitting RSS')
rss_fit <- susie_rss(zhat, R, n=n, estimate_prior_variance=estimate_prior_variance)
toc()

# add cs with consistent interfact to other methods
rss_fit$cs <- logisticsusie::compute_cs(rss_fit$alpha)

saveRDS(rss_fit, snakemake@output[[1]])
