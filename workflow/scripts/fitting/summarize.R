library(dplyr)
library(stringr)
library(reticulate)
library(tidyr)
library(tictoc)
library(ggplot2)



extract_sim_id <- function(path){
    path %>%
        stringr::str_split_1('/') %>%
        tail(2) %>%
        head(1)
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
rename_lbf <- function(fit){
    fit$lbf_ser <- fit$lbf  # rename lbf of each SER
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

load_gibss <- function(paths){
    # load GIBSS
    tic('loading GIBSS')
    gibss <- purrr::map(paths, readRDS)
    sim_ids <- purrr::map_chr(paths, extract_sim_id)
    gibss_tbl <- tibble(sim_id = sim_ids, fit = gibss)
    toc()
    return(gibss_tbl)
}

load_rss <- function(paths){
    # load RSS
    tic('load RSS')
    rss <- purrr::map(paths, readRDS)
    sim_ids <- purrr::map_chr(paths, extract_sim_id)
    rss_tbl <- tibble(sim_id = sim_ids, fit = fit) %>%
        rowwise() %>%
        mutate(rss = list(list_susie(rss))) %>%
        ungroup()
    toc()
    return(rss_tbl)
}

load_linear_susie <- function(paths){
    # load linear susie
    linear_susie <- purrr::map(paths, readRDS)
    linear_susie <- purrr::map(linear_susie, unclass)
    sim_ids <- purrr::map_chr(paths, extract_sim_id)
    linear_susie_tbl <- tibble(sim_id = sim_ids, fit = linear_susie) %>%
        rowwise() %>%
        mutate(fit = list(rename_lbf(fit))) %>%
        ungroup()
    return(linear_susie_tbl)
}

load_fits <- function(paths){
    path <- paths[1]
    method <- path %>%
        str_split_1('/') %>%
        tail(1) %>%
        str_split_1('[.]') %>%
        head(1)

    if(str_detect(method, 'gibss')){
        loader <- load_gibss
    } else if(str_detect(method, 'linear_susie')){
        loader <- load_linear_susie
    } else if(str_detect(method, 'rss')){
        loader <- load_rss
    } else{
        warning("couldn't match method name to loader")
    }

    fit_tbl <- loader(paths) %>% mutate(method = method)
    return(fit_tbl)  
}

manifest <- readr::read_delim(snakemake@input$manifest, delim='\t')

tic('load simulation')
sim_paths <- snakemake@input$sims
print(sim_paths)
sim <- purrr::map(sim_paths, readRDS) #%>% purrr::flatten()
sim_id <- purrr::map_chr(sim_paths, extract_sim_id)
sim_tbl <- tibble(sim_id = sim_id, sim = sim)
toc()

tic('load LD')
R <- readRDS(snakemake@input$ld)
toc()

tic('loading models')
fit_paths <- snakemake@input$fits
print(fit_paths)
fit_tbl <- load_fits(fit_paths)
toc()


tic('make full table')
res <- manifest %>%
    right_join(fit_tbl) %>%
    left_join(sim_tbl)
toc()

unclass <- function(obj){
    class(obj) <- 'list'
    return(obj)
}

compute_cs <- function(alpha){
    cs <- logisticsusie::compute_cs(alpha)
    cs2 <- purrr::map(cs, unclass)
    return(cs2)
}

compute_purity <- function(cs, R){
    min(abs(R[cs, cs]))
}

tic('compute cs summary')
res_cs <- res %>%
    hoist(sim, 'idx') %>%
    hoist(fit, 'lbf_ser') %>%
    rowwise() %>% mutate(cs = list(compute_cs(fit$alpha))) %>% ungroup() %>%
    select(-c(sim, fit)) %>%
    unnest_longer(c(cs, lbf_ser)) %>%
    unnest_wider(cs) %>%
    rowwise() %>%
    mutate(found = list(idx %in% cs), covered = any(found), n_found = sum(found), purity = compute_purity(cs, R)) %>%
    ungroup()
toc()

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
tic('compute pip summary')
res_pip <- res %>%
    hoist(sim, 'idx') %>%
    rowwise() %>% mutate(alpha = list(fit$alpha), cs = list(compute_cs(fit$alpha))) %>% ungroup() %>%
    mutate(is_null = purrr::map_lgl(alpha, ~any(is.null(.x)))) %>%
    dplyr::filter(!is_null) %>%
    rowwise() %>%
    mutate(
        cs_size = list(purrr::map_dbl(cs, ~.x$size)), 
        keep = list(cs_size < 200),
        pip = list(compute_pip2(alpha, keep)),
        causal = list(as.integer(1:length(pip) %in% idx))
    ) %>%
    select(-c(sim, alpha, cs, fit, is_null, cs_size, keep, idx))
toc()

saveRDS(res_cs, snakemake@output$cs_summary)
saveRDS(res_pip, snakemake@output$pip_summary)