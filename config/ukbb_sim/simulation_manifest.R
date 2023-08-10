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
    rep = 1:10,
    q = 1:3, 
    k = c(0.2),
    h2 = c(0.5), 
    rho = c(0.02, 0.04, 0.1, 0.2)
  ) %>%
  filter(rho * q < 1.) %>%
  rowwise() %>%
  mutate(
    sim_hash = rlang::hash(c(region, rep, q, k, h2, rho)),
    sim_id = stringr::str_extract(sim_hash, "^.{5}"),
    seed = strtoi(sim_id, 16)
  )

manifest %>% readr::write_delim('config/ukbb_sim/ukbb_sim_manifest.tsv', delim='\t')
