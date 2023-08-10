library(dplyr)
library(glue)
library(Matrix)

#' Simulate case control data
#' 
#' Normalized effects sizes are i.i.d normal 
#'    that is, variance inversely proportional to MAF.
#' We model an unliked genetic component with a normal random effect-- 
#'    consistent with a highly polygenic background.
#' Genetic liability (on log odds scale) is scaled and shifted 
#'    to achieve desired "h2 of liability" and disease prevalance.
#' Simulation is meant to reflect unbiased sampling from a population
#'    which is reasonable for biobank data?
#'    
#' @param X genotype matrix
#' @param q number of casual variants
#' @param k disease prevelance in cohort
#' @param h2 heritability of liability on the log-odds scale
#' @param rho fraction of heritabaility of liability for each causal variant 
#'    NOTE: rho * q < 1
simulate_case_control <- function(X, q=3, k=0.05, h2=0.2, rho=0.1, sim_id){
  set.seed(strtoi(sim_id, 16))

  n <- nrow(X) 
  p <- ncol(X)
  idx <- sample(p, q, replace=F) # sample causal variants
  beta <- abs(rnorm(q)) * sqrt(rho / (1 - q * rho))
  Xsub <- as.matrix(X[, idx])
  Xscaled <- scale(X[, idx])
  l1 <- (Xscaled %*% beta)[, 1] # risk of causal variants
  l2 <-  rnorm(n) # risk due to other variants
  genetic_liability <- scale(l1 + l2)[,1] *  sqrt((h2 * pi^2/3) / (1 - h2))
  
  # set intercept empirically to achieved desired prevalence
  liability <- genetic_liability + rlogis(n, location = 0, scale=1)
  beta0 <-  - quantile(liability, 1 - k) # get intercept
  
  # compute normalized effects (being lazy, just run the regression)
  logit <- beta0 + genetic_liability
  beta <- unname(lm(logit ~ Xscaled)$coef)
  
  sim <- list(
    k = k,
    q = q,
    h2 = h2,
    rho = rho,
    idx = idx,
    intercept = beta[1],
    beta = tail(beta, -1),
    logit = logit,
    y = rbinom(length(logit), 1, 1 / (1 + exp(-logit))),
    sim_id = sim_id
  )
  return(sim)
}

region = snakemake@wildcards[['region']]
print(region)

manifest <- readr::read_delim('config/ukbb_sim/ukbb_sim_manifest.tsv')
manifest <- manifest %>% filter(region == region)

X <- readRDS(glue('results/ukbb_geno/{region}/genotype.rds'))
sim <- manifest %>%
  rowwise() %>%
  mutate(sim = list(simulate_case_control(X=X, q=q, k=k, h2=h2, rho=rho, sim_id=sim_id)))

sims <- sim$sim
names(sims) <- sim$sim_id
saveRDS(sims, glue('results/ukbb_geno/{region}/sim.rds'))

sim_py <- reticulate::r_to_py(sims)
reticulate::py_save_object(sim_py, glue('results/ukbb_geno/{region}/sim.pkl'))

