source(".Rprofile")
library(dplyr)
library(glue)
library(Matrix)
library(tictoc)

#' Simulate case control data
#'
#' Normalized effects sizes are i.i.d normal
#'    that is, variance inversely proportional to MAF.
#' We model an unliked genetic component with a normal random effect--
#'    consistent with a highly polygenic background.
#' Genetic liability (on log odds scale) is scaled and shifted
#'    to achieve desired "h2 of liability" and disease prevalance.
#' Simulation is meant to reflect unbiased sampling from a population
#'    i.e. unbalanced case/control ratio
#'
#' @param X genotype matrix of causal SNPs
#' @param q number of casual variants
#' @param k disease prevelance in cohort
#' @param h2 heritability of liability on the log-odds scale
#' @param rho fraction of heritabaility of liability in locus, split equally among causal variants
simulate_case_control <- function(X, k = 0.05, h2 = 0.2, rho = 0.1) {
  n <- nrow(X)
  q <- ncol(X)

  # compute prior variance
  nu <- pi^2 / 3 * h2 / (1 - h2) # total variance of the genetic component
  var_res <- nu * (1 - rho) # total variance of unlinked genetic component
  sigma2 <- rho * nu / q # prior variance of effect

  Xscaled <- scale(X)
  beta <- rnorm(q) * sqrt(sigma2) #
  genetic_liability <- (Xscaled %*% beta)[, 1] + rnorm(n) * sqrt(var_res)

  # set intercept to achieved desired prevalence
  liability <- genetic_liability + rlogis(n, location = 0, scale = 1)
  beta0 <- -quantile(liability, 1 - k) |> unname() # get intercept

  # compute normalized effects (being lazy, just run the regression)
  logit <- beta0 + genetic_liability

  sim <- list(
    k = k,
    q = q,
    h2 = h2,
    rho = rho,
    prior_variance = sigma2,
    intercept = beta0,
    beta = beta,
    logit = logit,
    y = rbinom(length(logit), 1, 1 / (1 + exp(-logit)))
  )
  return(sim)
}


get_region <- function(pvar) {
  chrom <- pgenlibr::GetVariantChrom(pvar, 1)
  lower <- pgenlibr::GetVariantPos(pvar, 1)
  upper <- pgenlibr::GetVariantPos(pvar, pgenlibr::GetVariantCt(pvar))
  region <- glue::glue("chr{chrom}:{lower}-chr{chrom}:{upper}")
  return(region)
}

params <- snakemake@params$sim_params
print(params)

# SIM_ID <- snakemake@params$sim_id
set.seed(as.integer(params$seed))

# load PGEN
pgen_path <- glue::glue("{params$pfile}.pgen")
pvar_path <- glue::glue("{params$pfile}.pvar")
psam_path <- glue::glue("{params$pfile}.psam")
pvar <- pgenlibr::NewPvar(pvar_path)
pgen <- pgenlibr::NewPgen(pgen_path, pvar = pvar)

# q <- snakemake@params$num_causal_variables
num_causal_variables <- as.integer(params$num_causal_variables)
num_snps <- pgenlibr::GetVariantCt(pvar)
variant_idx <- sample(1:num_snps, num_causal_variables)
causal_rsids <- purrr::map_chr(variant_idx, ~ pgenlibr::GetVariantId(pvar, .x))
X <- pgenlibr::ReadList(pgen, variant_idx)

# params <- snakemake@params
# params <- list(k = 0.05, h2 = 0.2, rho = 0.1)
sim <- with(params, simulate_case_control(X, k, h2, rho))

# add details to sim
sim <- c(sim, list(
  num_causal_variables = num_causal_variables,
  region = get_region(pvar),
  pvar_path = pvar_path,
  pgen_path = pgen_path,
  causal_idx = variant_idx,
  causal_rsids = causal_rsids
))

saveRDS(sim, snakemake@output$rds)

library(data.table)
psam <- fread(psam_path)
psam[, Binary := (sim$y + 1)]
psam[, Continuous := as.numeric(sim$y - mean(sim$y))]
psam[, SEX := NULL]
fwrite(psam, file = snakemake@output$pheno, sep = "\t", quote = F)

sam <- read.table(psam_path, header = T)
# to python
# tic("\t Python")
# sim_py <- reticulate::r_to_py(sim)
# reticulate::py_save_object(sim_py, snakemake@output$pkl)
# toc()
