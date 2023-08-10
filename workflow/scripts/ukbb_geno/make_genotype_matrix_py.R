library(data.table)
library(Matrix)

read.geno.raw <- function (geno.file) {
  geno <- fread(geno.file,sep = "\t",header = TRUE,stringsAsFactors = FALSE,
		showProgress = TRUE)
  class(geno)    <- "data.frame"
  ids            <- with(geno,paste(FID,IID,sep = "_"))
  geno           <- geno[-(1:6)]
  rownames(geno) <- ids
  geno           <- as.matrix(geno)
  storage.mode(geno) <- "double"
  return(geno)
}

# region <- 'chr1_100794065_101983457'
region <- snakemake@wildcards[['region']]

print(region)
region2 <- stringr::str_replace_all(region, '_', '.')
file <- paste0('bloodcells_', region2, '.raw.gz')
path <- '/project2/mstephens/yuxin/ukb-bloodcells/genotypes'
geno <- read.geno.raw(paste(path, file, sep='/'))

# select 50k random samples
set.seed(1)
ind_idx <- sort(sample(nrow(geno), 50000))
X <- geno[ind_idx,]

# filter to SNPs with MAF > 0.01
allele_frequencies <- Matrix::colMeans(X) / 2
maf <- pmin(allele_frequencies, 1 - allele_frequencies)
maf_filter <- (maf > 0.01)
X <- X[, maf_filter]

Xsp <- as(X, "sparseMatrix")
saveRDS(Xsp, snakemake@output[[2]])

# organize data
ind_ids <- rownames(X)
snp_ids <- colnames(X)
X <- as(X, "sparseMatrix")
res <- list(ind_ids = ind_ids, snp_ids = snp_ids, X = X)
res_py <- reticulate::r_to_py(res)
reticulate::py_save_object(res_py, snakemake@output[[1]])

