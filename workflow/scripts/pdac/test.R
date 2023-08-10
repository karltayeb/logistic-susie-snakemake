a <- rnorm(1000)
saveRDS(a, file = 'results/rnorm3.rds')
saveRDS(a, file = snakemake@output[[1]])
# saveRDS(a, file = 'results/rnorm.rds') 