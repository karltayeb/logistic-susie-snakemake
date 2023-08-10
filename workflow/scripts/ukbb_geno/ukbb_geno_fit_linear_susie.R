library(susieR)
library(reticulate)
library(Matrix)
library(glue)
library(tictoc)

region <- snakemake@wildcards[['region']]
sim_id <- snakemake@wildcards[['region']]

# 1. load data

print('loading data')
sim_path <- snakemake@input[['sim_path']]
sims <- readRDS(sim_path)

genotype_path <- snakemake@input[['genotype_path']]
X <- readRDS(genotype_path)

estimate_prior_variance = (snakemake@params[['estimate_prior_variance']] == 'true')

for(sim_id in names(sims)){
    print(sim_id)

    y <- as.vector(sims[[sim_id]]$y)

    # 3. fit susie_rss
    print(glue('\tfitting linear SuSiE: estimate_prior_variance = {estimate_prior_variance}'))

    tic()
    susie_fit <- susie(X, y, L=5)
    toc()

    # add cs with consistent interfact to other methods
    susie_fit$cs <- logisticsusie::compute_cs(susie_fit$alpha)

    save_path <- glue('results/ukbb_geno/{region}/{sim_id}/linear_susie.rds')
    print(glue('\t saving to: {save_path}'))
    saveRDS(susie_fit, save_path)
}

