# Yusha example -------
# this is an example where binsusie essentially fits the null model.
# in the below analysis we fit various SERs and SuSiEs
# to see if our new strategies improve the situation
library(dplyr)


#' extract file name from path-- without extention
path2file <- function(path, keep_extention=FALSE){
  file <- stringr::str_split_1(path, '/') %>%
    tail(1)
    
  if(!keep_extention){
    file <- file %>%
      stringr::str_split_1('\\.') %>%
      head(1)
  }
  return(file)
}


#' Load background genes for pdac example
load_pdac_example <- function(background_path, list_paths){
  # load background genes
  genes_in_study <- read.csv2(background_path, header = F)$V1
  hs <- org.Hs.eg.db::org.Hs.eg.db


  # map to entrez IDs
  idmap <- AnnotationDbi::select(hs, keys = genes_in_study,
                                columns = c('SYMBOL', 'ENTREZID'),
                                keytype = 'SYMBOL')
  genes_in_study_entrez <- idmap$ENTREZID[!is.na(idmap$ENTREZID)]


  # load gene lists, convert to ENTREZ
  gene_lists <- list()
  gene_lists_symbol <- list()
  for(path in list_paths){
    name = path2file(path)
    gene_list <- read.csv2(path, header=F)$V1
    gene_list_entrez <- idmap %>%
      dplyr::filter(!is.na(ENTREZID)) %>%
      dplyr::filter(SYMBOL %in% gene_list) %>% 
      {.$ENTREZID}
    gene_lists[[name]] <- gene_list_entrez
    gene_lists_symbol[[name]] <- gene_list
  }

  data <- list(background = genes_in_study_entrez, gene_lists = gene_lists)
  return(data)
}

# hard coded data paths
background = 'resources/yusha_sc_tumor/pdac/gene_list.txt'
files = list.files('resources/yusha_sc_tumor/pdac/', full.names=T)
lists = files[grepl('factor', files)]

# load data, convert to ensemblIDs
data <- load_pdac_example(background, lists)
saveRDS(data, file = 'results/pdac/pdac_data.rds')


