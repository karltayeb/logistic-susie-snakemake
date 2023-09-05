library(dplyr)
library(reticulate)
library(tidyr)
library(tictoc)
library(ggplot2)

manifest <- readr::read_delim('config/ukbb_sim/ukbb_sim_manifest.tsv', delim='\t')


np_cast <- function(obj){
    np <- import('numpy')
    items <- names(obj)
    for(item in items){
        if(typeof(obj[[item]]) == 'environment'){
            obj[[item]] <- np$array(obj[[item]])
        }
    }
    return(obj)
}

np_cast_ser <- function(ser){
    ser <- np_cast(ser)
    ser$cs <- np_cast(ser$cs)
    return(ser)
}

extract_sim_id <- function(path){
    path %>%
        stringr::str_split_1('/') %>%
        tail(2) %>%
        head(1)
}


tic('load simulation')
sim_paths <- Sys.glob('results/ukbb_geno/chr1*/*/sim.rds')
sim <- purrr::map(sim_paths, readRDS) %>% purrr::flatten()
sim_tbl <- tibble(sim_id = names(sim), sim = sim)
toc()

tic('load LD')
R <- readRDS(snakemake@input$ld)
toc()

# make SER ABF-- we've already stored alpha and lbf approximated via ABF
ser_use_abf <- function(ser){
    ser$lbf <- ser$log_abf
    ser$alpha <- ser$alpha_abf
    ser$cs <- ser$cs_abf
    return(ser)
}

#' wrap CS in a list
encapsulate_cs <- function(ser){
    ser$cs <- list(L0 = ser$cs)
    return(ser)
}

#' Make class of CS a list
fix_cs <- function(cs){
    cs$cs_size <- cs$size
    cs$target_coverage <- cs$requested_coverage
    cs$size <- NULL
    cs$requested_coverage <- NULL
    class(cs) <- 'list'
    return(cs)
}

#' fix susie so that it plays well with tidy operations
list_susie <- function(fit){
    new_cs <- purrr::map(fit$cs, fix_cs)
    names(new_cs) <- names(fit$cs)
    fit$cs <- new_cs
    fit$lbf_ser <- fit$lbf  # rename lbf of each SER
    class(fit) <- 'list'
    return(fit)
}

#' Go from one indexing to zero indexing
one_index_cs <- function(cs){
    cs$cs <- cs$cs + 1
    return(cs)
}

#' One index a list of CSs
one_index_css <- function(fit){
    new_cs <- purrr::map(fit$cs, one_index_cs)
    names(new_cs) <- names(fit$cs)
    fit$cs <- new_cs
    return(fit)
}

load_sers <- function(paths){
    # load SERs
    tic('load sers')
    sers <- purrr::map(paths, ~np_cast_ser(py_load_object(.x)))
    sim_ids <- purrr::map_chr(paths, extract_sim_id)
    ser_tbl <- tibble(sim_id = sim_ids, ser = sers)

    ser_tbl <- ser_tbl %>%
        rowwise() %>%
        mutate(ser_abf = list(ser_use_abf(ser))) %>%
        mutate(
            ser = list(one_index_css(encapsulate_cs(ser))), 
            ser_abf = list(one_index_css(encapsulate_cs(ser_abf)))) %>%
        ungroup()
    toc()
    return(ser_tbl)
}

load_gibss <- function(paths){
    # load GIBSS
    tic('loading GIBSS')
    gibss <- purrr::map(paths, ~py_load_object(.x))
    sim_ids <- purrr::map_chr(paths, extract_sim_id)
    gibss_tbl <- tibble(sim_id = sim_ids, gibss = gibss) %>%
        rowwise() %>%
        mutate(gibss = list(one_index_css(gibss))) %>%
        ungroup()
    toc()
    return(gibss_tbl)
}

load_rss <- function(paths){
    # load RSS
    tic('load RSS')
    rss <- purrr::map(paths, readRDS)
    sim_ids <- purrr::map_chr(paths, extract_sim_id)
    rss_tbl <- tibble(sim_id = sim_ids, rss = rss) %>%
        rowwise() %>%
        mutate(rss = list(list_susie(rss))) %>%
        ungroup()
    toc()
    return(rss_tbl)
}

load_linear_susie <- function(paths){
    # load linear susie
    linear_susie <- purrr::map(paths, readRDS)
    sim_ids <- purrr::map_chr(paths, extract_sim_id)
    linear_susie_tbl <- tibble(sim_id = sim_ids, linear_susie = linear_susie) %>%
        rowwise() %>%
        mutate(linear_susie = list(list_susie(linear_susie))) %>%
        ungroup()
    return(linear_susie_tbl)
}

ser_paths <- Sys.glob('results/ukbb_geno/chr1*/*/ser_default.pkl')
labf_paths <- Sys.glob('results/ukbb_geno/chr1*/*/gibss_labf_default.pkl')
abf_paths <- Sys.glob('results/ukbb_geno/chr1*/*/gibss_abf_default.pkl')
rss_l1_eb_paths <- Sys.glob('results/ukbb_geno/chr1*/*/rss_L1-eb.rds')
rss_eb_paths <- Sys.glob('results/ukbb_geno/chr1*/*/rss_L5-eb.rds')
susie_paths <- Sys.glob('results/ukbb_geno/chr1*/*/linear_susie_default.rds')
susie_l1_paths <- Sys.glob('results/ukbb_geno/chr1*/*/linear_susie_L1.rds')

ser_tbl <- load_sers(ser_paths)
gibbs_abf_tbl <- load_gibss(abf_paths)
gibbs_labf_tbl <- load_gibss(labf_paths)
rss_eb_tbl <- load_rss(rss_eb_paths)
rss_l1_eb_tbl <- load_rss(rss_l1_eb_paths)
linear_susie_tbl <- load_linear_susie(susie_paths)
linear_susie_l1_tbl <- load_linear_susie(susie_l1_paths)

# load
res <- manifest %>%
    left_join(sim_tbl) %>%
    left_join(ser_tbl) %>%
    left_join(rename(gibbs_labf_tbl, gibss_labf = gibss)) %>%
    left_join(rename(gibbs_abf_tbl, gibss_abf = gibss)) %>%
    left_join(rss_eb_tbl) %>%
    left_join(rename(rss_l1_eb_tbl, rss_l1 = rss)) %>%
    left_join(linear_susie_tbl) %>%
    left_join(rename(linear_susie_l1_tbl, linear_susie_l1 = linear_susie))

compute_purity <- function(cs, R){
    min(abs(R[cs, cs]))
}

not_methods <- c(intersect(colnames(res), colnames(manifest)), 'sim')
res_cs <- res %>%
    pivot_longer(-any_of(not_methods), names_to = 'method', values_to = 'fit') %>%
    hoist(sim, 'idx') %>%
    hoist(fit, 'cs') %>%
    hoist(fit, 'lbf_ser') %>%
    select(-c(sim, fit)) %>%
    unnest_longer(c(cs, lbf_ser)) %>%
    unnest_wider(cs) %>%
    rowwise() %>%
    mutate(found = list(idx %in% cs), covered = any(found), n_found = sum(found), purity = compute_purity(cs, R)) %>%
    ungroup()

# Power/FDR
compute_pip2 <- function(alpha, keep){
    if(length(dim(alpha)) < 2){
        alpha <- matrix(alpha, nrow=1)
        keep <- TRUE
    }
    pips <- logisticsusie::compute_pip(alpha[keep,])
    return(pips)
}

# filter out CSs with large size and then compute PIPs
res_pip <- res %>%
    pivot_longer(-any_of(not_methods), names_to = 'method', values_to = 'fit') %>%
    hoist(sim, 'idx') %>%
    hoist(fit, alpha = 'alpha', cs = 'cs') %>%
    mutate(is_null = purrr::map_lgl(alpha, ~any(is.null(.x)))) %>%
    dplyr::filter(!is_null) %>%
    rowwise() %>%
    mutate(
        cs_size = list(purrr::map_dbl(cs, ~.x$cs_size)), 
        keep = list(cs_size < 200),
        pip = list(compute_pip2(alpha, keep)),
        causal = list(as.integer(1:length(pip) %in% idx))
    ) %>%
    select(-c(sim, alpha, cs, fit, is_null, cs_size, keep, idx))

saveRDS(res_cs, snakemake@output$cs_summary)
saveRDS(res_pip, snakemake@output$pip_summary)