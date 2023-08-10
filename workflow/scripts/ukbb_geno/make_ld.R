library(glue)
library(Matrix)

# 2. load genotype matrix, make LD matrix
region <- snakemake@wildcards[['region']]
genotype_path <- glue('results/ukbb_geno/{region}/genotype.rds')
genotype <- readRDS(genotype_path)

sparse.cor <- function(x){
    n <- nrow(x)

    cMeans <- colMeans(x)
    cSums <- colSums(x)

    # Calculate the population covariance matrix.
    # There's no need to divide by (n-1) as the std. dev is also calculated the same way.
    # The code is optimized to minize use of memory and expensive operations
    covmat <- tcrossprod(cMeans, (-2*cSums+n*cMeans))
    crossp <- as.matrix(crossprod(x))
    covmat <- covmat+crossp

    sdvec <- sqrt(diag(covmat)) # standard deviations of columns
    covmat/crossprod(t(sdvec)) # correlation matrix
}

library(tictoc)

tic()
print(dim(genotype))
R <- cor(as.matrix(genotype))
toc()

saveRDS(R, snakemake@output[[1]])