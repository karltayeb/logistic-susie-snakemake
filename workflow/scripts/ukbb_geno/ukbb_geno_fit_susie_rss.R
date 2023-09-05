library(susieR)
library(reticulate)
library(Matrix)
library(glue)
library(tictoc)

region <- snakemake@wildcards[['region']]
sim_id <- snakemake@wildcards[['region']]
fit_params <- snakemake@params$fit_params

# 0. load simulation, compute hyperparameters
sim_path <- snakemake@input[['sim_path']]
sim <- readRDS(sim_path)[[1]]
phi <- mean(sim$y)
var_y <- 1 / (phi * (1-phi))
n <- length(sim$y)

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
print('fitting SuSiE RSS')


if(is.null(fit_params$var_y)){
    fit_params$var_y <- var_y
}
print(fit_params)


tic('fitting RSS')
data <- list(z=zhat, R=R, n=n)
args <- c(data, fit_params)
rss_fit <- rlang::exec(susie_rss, !!!args) #(zhat, R, n=n, estimate_prior_variance=estimate_prior_variance, L=L)
toc()


# add cs with consistent interfact to other methods
rss_fit$cs <- logisticsusie::compute_cs(rss_fit$alpha)
saveRDS(rss_fit, snakemake@output[[1]])
