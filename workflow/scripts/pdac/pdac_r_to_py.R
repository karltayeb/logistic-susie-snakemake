# convert pdac data and simulations to python using reticulate
library(reticulate)
library(glmnet)

pdac_data = readRDS('results/pdac/pdac_data.rds')
pdac_msigdb_X = readRDS('results/pdac/pdac_msigdb_X.rds')
pdac_lasso_fits = readRDS('results/pdac/pdac_lasso_fits.rds')
pdac_lasso_sims = readRDS('results/pdac/pdac_lasso_sims.rds')

# extract coefficients for each factor
np <- import('numpy')
coefs <- purrr::map(pdac_lasso_fits, ~np$array(coef(.x, s='lambda.1se')[,1]))
names(coefs) <- names(pdac_lasso_fits)

# send data to py
X_dense <- as.matrix(pdac_msigdb_X)
pdac_data_py <- r_to_py(pdac_data)

genes = rownames(pdac_msigdb_X)
genesets = colnames(pdac_msigdb_X)

factors = purrr::map(names(coefs), ~list(
    list = pdac_data$gene_lists[[.x]], 
    y = np$array(as.integer(genes %in% pdac_data$gene_lists[[.x]])),
    Ysim = do.call(rbind, pdac_lasso_sims[[.x]]$y) + 0.0,
    coef = coefs[[.x]]))
names(factors) <- names(coefs)

to_py <- list(
    genes = genes,
    genesets = genesets,
    X = X_dense,
    factors = factors
)
py <- r_to_py(to_py)
py_save_object(py, 'results/pdac/pdac_data_py.pkl')
