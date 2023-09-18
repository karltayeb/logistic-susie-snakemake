library(susieR)
library(reticulate)
library(Matrix)
library(glue)
library(tictoc)

tic('loading data')
X <- readRDS(snakemake@input$X)
sim <- readRDS(snakemake@input$sim)
params <- snakemake@params$fit_params
toc()

tic('fitting susie')
args <- c(list(X=X, y=sim$y), snakemake@params$fit_params)
susie_fit <- rlang::exec(susie, !!!args)
toc()

save_path <- snakemake@output$fit
tic(glue('saving to: {save_path}'))
saveRDS(susie_fit, save_path)
toc()