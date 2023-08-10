
# for static branching over factors/gene lists
make_pdac_tbl <- function(){
  # get file names of gene lists
  files <- list.files('data/yusha_sc_tumor/pdac/')
  keep <- purrr::map_lgl(files, ~stringr::str_detect(.x, 'factor'))
  files <- files[keep]
  names <- purrr::map_chr(files, ~strsplit(.x, '\\.')[[1]][1])

  pdac <- dplyr::as_tibble(data.frame(study='pdac', factor = names))
  return(pdac)
}

geneSet2X <- function(geneSet){
  # construct sparse matrix
  genes <- unique(geneSet$gene)
  sets <- unique(geneSet$geneSet)
  
  geneSet$row <- match(geneSet$gene, genes)
  geneSet$col <- match(geneSet$geneSet, sets)
  
  X <- Matrix::sparseMatrix(
    i = geneSet$row, 
    j = geneSet$col,
    x = rep(1, nrow(geneSet)),
    dimnames = list(genes, sets)
  )
}

load_msigdb_X <- function(background){
  msigdb <- gseasusie::load_gene_sets(c('c2', 'h')) # load mSigDB "C2" and "H"
  geneSet <- 
    purrr::map_dfr(msigdb, ~purrr::pluck(.x, "geneSet")) %>% # concat "C2" and "H"
    dplyr::select(geneSet, gene) %>%
    dplyr::distinct() %>% 
    dplyr::filter(gene %in% background) %>% # subset to genes in study
    dplyr::group_by(geneSet) %>%
    dplyr::mutate(n = length(gene)) %>%
    dplyr::filter(n >= 10) %>% # filter to gene sets of a minimum size
    dplyr::select(-n)
  X <- geneSet2X(geneSet) # make gene set matrix
  return(X)
}

fit_lasso <- function(X, gene_list){
  y <- as.integer(rownames(X) %in% gene_list)
  lasso_fit <- glmnet::cv.glmnet(X, y, family = "binomial")
  return(lasso_fit)
}

fit_lasso_all_gene_lists <- function(X, gene_lists){
  lasso_fits <- purrr::map(gene_lists, ~fit_lasso(X, .x))
  names(lasso_fits) <- names(gene_lists)
  return(lasso_fits)
}

simulate_from_lasso <- function(X, lasso_fit, reps){
  coef_1se <- coef(lasso_fit, s = 'lambda.1se')[,1]
  logit_1se <- (X %*% tail(coef_1se, -1))[,1] + coef_1se[1]
  ysim_1se <- purrr::map(1:reps, ~rbinom(length(logit_1se), 1, 1/(1 + exp(-logit_1se))))
  res <- list(y = ysim_1se, logit = logit_1se, intercept = coef_1se[1], beta = tail(coef_1se, -1))
  return(res)
}

simulate_from_lasso_map <- function(X, lasso_fits, reps){
  simulations <- purrr::map(lasso_fits, ~simulate_from_lasso(X, .x, reps))
  names(simulations) <- names(lasso_fits)
  return(simulations)
}

pdac_msigdb_X = load_msigdb_X(pdac_data$background)

pdac_tbl <- make_pdac_tbl()
pdac_targets <- list(
    tar_target(pdac_data, load_pdac_example()),
    tar_target(),
    tar_map(pdac_tbl, 
      tar_target(ibss_fit_marginal_only, fit_ibss_marginal_only(pdac_data$gene_lists[[factor]], pdac_msigdb_X)),
      tar_target(ibss_fit, fit_ibss_select_L(pdac_data$gene_lists[[factor]], pdac_msigdb_X)),
      tar_target(ibss_fit_lasso_init, fit_ibss_lasso_init(pdac_data$gene_lists[[factor]], pdac_msigdb_X))
    ),
    tar_target(pdac_lasso_fits, fit_lasso_all_gene_lists(pdac_msigdb_X, pdac_data$gene_lists)),
    tar_target(pdac_lasso_sims, simulate_from_lasso_map(pdac_msigdb_X, pdac_lasso_fits, 20))
)