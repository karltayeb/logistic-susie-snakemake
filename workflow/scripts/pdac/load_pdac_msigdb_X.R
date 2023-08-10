library(dplyr)

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

# hard coded
pdac_data = readRDS('results/pdac/pdac_data.rds')
pdac_msigdb_X = load_msigdb_X(pdac_data$background)
saveRDS(pdac_msigdb_X, file = 'results/pdac/pdac_msigdb_X.rds')

# pdac_data = readRDS(snakemake@input[[1]])
# pdac_msigdb_X = load_msigdb_X(pdac_data$background)
# saveRDS(pdac_msigdb_X, file = snakemake@output[[1]])