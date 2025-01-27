library(dplyr)

#' to run the RSS comparison pipeline
parameter_table <- tibble::tribble(
  ~sim_id, ~pfile, ~n_causal_variants, ~disease_prevalence, ~h2, ~causal_snp_h2, ~seed,
  "sim1", "resources/ukb-geno/bloodcells_chr14.100926769.101431881", 1, 0.2, 0.4, 0.005, 1L,
  "sim2", "resources/ukb-geno/bloodcells_chr14.100926769.101431881", 2, 0.2, 0.4, 0.01, 2L,
  "sim3", "resources/ukb-geno/bloodcells_chr14.100926769.101431881", 3, 0.2, 0.4, 0.02, 3L,
)
