library(dplyr)
make_region <- function(path){
  path[1] %>%
    stringr::str_split_1('bloodcells_') %>% tail(1) %>%
    stringr::str_split_1('\\.') %>% head(3) %>%
    paste0(collapse='_')
}
regions <- list.files('/project2/mstephens/yuxin/ukb-bloodcells/genotypes') %>%
  purrr::map_chr(make_region) %>%
  unique()
chr1_regions <- regions[stringr::str_detect(regions, 'chr1_')] %>% head(1)

manifest <- tidyr::crossing(
    region = chr1_regions,
    rep = 1:5,
    q = 1:3, 
    k = c(0.01, 0.05, 0.1, 0.2),
    h2 = c(0.5), 
    rho = c(0.005, 0.01, 0.02, 0.04, 0.1),
    name = 'main'
  ) %>%
  mutate(gamma = rho*h2) %>%
  filter(rho * q < 1.) %>%
  rowwise() %>%
  mutate(
    sim_hash = rlang::hash(c(region, rep, q, k, h2, rho)),
    sim_id = stringr::str_extract(sim_hash, "^.{5}"),
    seed = strtoi(sim_id, 16)
  )
manifest %>% readr::write_delim('config/ukbb_sim/ukbb_main_manifest.tsv', delim='\t')

ser_manifest <- tidyr::crossing(
    region = 'chr1_100794065_101983457',
    rep = 1:5,
    q = 1, 
    k = c(0.005, 0.01, 0.05, 0.1, 0.2),
    gamma = cumprod(c(0.001, rep(2, 9))),
    h2 = cumprod(c(0.01, rep(2, 5))),
    name = 'ser'
  ) %>%
  mutate(rho = gamma/h2) %>%
  filter(rho < 1) %>%
  rowwise() %>%
  mutate(
    sim_hash = rlang::hash(c(region, rep, q, k, h2, rho)),
    sim_id = stringr::str_extract(sim_hash, "^.{5}"),
    seed = strtoi(sim_id, 16)
  )
ser_manifest %>% 
  readr::write_delim('config/ukbb_sim/ukbb_ser_manifest.tsv', delim='\t')

bind_rows(manifest, ser_manifest) %>%
  readr::write_delim('config/ukbb_sim/ukbb_sim_manifest.tsv', delim='\t')
