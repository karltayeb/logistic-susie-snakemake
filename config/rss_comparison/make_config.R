library(dplyr)

#' to run the RSS comparison pipeline
parameter_table <- tibble::tribble(
  ~sim_id, ~pfile, ~num_causal_variables, ~k, ~h2, ~rho, ~seed,
  "sim1", "resources/ukb-geno/bloodcells_chr14.100926769.101431881", 1, 0.2, 0.4, 0.05, 1L,
  "sim2", "resources/ukb-geno/bloodcells_chr14.100926769.101431881", 2, 0.2, 0.4, 0.05, 1L,
  "sim3", "resources/ukb-geno/bloodcells_chr14.100926769.101431881", 3, 0.2, 0.4, 0.05, 1L,
)
