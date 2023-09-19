library(dplyr)

manifest <- tidyr::crossing(
    X = c('X_bin_strong', 'X_bin_med', 'X_bin_weak'),
    rep = 1:10,
    beta0 = c(-5, -4, -3, -2),
    beta = c(0.5, 1., 2., 4.),
    q = c(1, 3, 5)) %>%
  rowwise() %>%
  mutate(
    sim_hash = rlang::hash(c(X, rep, beta0, beta, q)),
    sim_id = stringr::str_extract(sim_hash, "^.{5}"),
    seed = strtoi(sim_id, 16)
  )

manifest %>% readr::write_delim('config/gsea_sim/gsea_constant_manifest.tsv', delim='\t')

manifest %>% readr::write_delim('config/gsea_sim/gsea_sim_manifest.tsv', delim='\t')
