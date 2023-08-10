library(dplyr)
library(reticulate)
library(tidyr)

manifest <- readr::read_delim('config/ukbb_sim/ukbb_sim_manifest.tsv', delim='\t')

sim_path <- Sys.glob('results/ukbb_geno/chr1*/sim.rds')
ser_paths <- Sys.glob('results/ukbb_geno/chr1*/*/ser.pkl')
labf_paths <- Sys.glob('results/ukbb_geno/chr1*/*/gibss.pkl')
abf_paths <- Sys.glob('results/ukbb_geno/chr1*/*/gibss_abf.pkl')
rss_paths <- Sys.glob('results/ukbb_geno/chr1*/*/rss.rds')
susie_paths <- Sys.glob('results/ukbb_geno/chr1*/*/linear_susie.rds')

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


sim <- readRDS(sim_path)
sim_tbl <- tibble(sim_id = names(sim), sim = sim)

# load SERs
sers <- purrr::map(ser_paths, ~np_cast_ser(py_load_object(.x)))
sim_ids <- purrr::map_chr(ser_paths, extract_sim_id)
ser_tbl <- tibble(sim_id = sim_ids, ser = sers)

# make SER ABF-- we've already stored alpha and lbf approximated via ABF

ser_use_abf <- function(ser){
    ser$lbf <- ser$log_abf
    ser$alpha <- ser$alpha_abf
    ser$cs <- ser$cs_abf
    return(ser)
}

ser_tbl <- ser_tbl %>%
    rowwise() %>%
    mutate(ser_abf = list(ser_use_abf(ser))) %>%
    ungroup()

# load GIBSS
gibss <- purrr::map(labf_paths, ~py_load_object(.x))
sim_ids <- purrr::map_chr(labf_paths, extract_sim_id)
gibss_tbl <- tibble(sim_id = sim_ids, gibss = gibss)

# load GIBSS - abf
gibss_abf <- purrr::map(abf_paths, ~py_load_object(.x))
sim_ids <- purrr::map_chr(abf_paths, extract_sim_id)
gibss_abf_tbl <- tibble(sim_id = sim_ids, gibss_abf = gibss_abf)


# load RSS
rss <- purrr::map(rss_paths, readRDS)
sim_ids <- purrr::map_chr(rss_paths, extract_sim_id)
rss_tbl <- tibble(sim_id = sim_ids, rss = rss)

# load linear susie
linear_susie <- purrr::map(susie_paths, readRDS)
sim_ids <- purrr::map_chr(susie_paths, extract_sim_id)
linear_susie_tbl <- tibble(sim_id = sim_ids, linear_susie = linear_susie)

# load
res <- manifest %>%
    left_join(sim_tbl) %>%
    left_join(ser_tbl) %>%
    left_join(gibss_tbl) %>%
    left_join(gibss_abf_tbl) %>%
    left_join(rss_tbl) %>%
    left_join(linear_susie_tbl)

# CS -------

# fix ser cs so it has the same structure as SuSiE CSs
encapsulate_cs <- function(ser){
    ser$cs <- list(L0 = ser$cs)
    return(ser)
}

one_index_cs <- function(cs){
    cs$cs <- cs$cs + 1
    return(cs)
}

one_index_css <- function(fit){
    new_cs <- purrr::map(fit$cs, one_index_cs)
    names(new_cs) <- names(fit$cs)
    fit$cs <- new_cs
    return(fit)
}

list_cs <- function(cs){
    cs$cs_size <- cs$size
    cs$target_coverage <- cs$requested_coverage
    cs$size <- NULL
    cs$requested_coverage <- NULL
    class(cs) <- 'list'
    return(cs)
}

list_susie <- function(fit){
    new_cs <- purrr::map(fit$cs, fix_cs)
    names(new_cs) <- names(fit$cs)
    fit$cs <- new_cs
    class(fit) <- 'list'
    return(fit)
}

res <- res %>%
    rowwise() %>%
    mutate(
        gibss = list(one_index_css(gibss)),
        gibss_abf = list(one_index_css(gibss_abf)),
        ser = list(one_index_css(encapsulate_cs(ser))),
        ser_abf = list(one_index_css(encapsulate_cs(ser_abf))),
        linear_susie = list(list_susie(linear_susie)),
        rss = list(list_susie(rss))
    ) %>%
    ungroup()

res_cs <- res %>%
    # select(-c(linear_susie, rss)) %>%
    pivot_longer(c(ser, ser_abf, gibss, gibss_abf, linear_susie, rss), names_to = 'method', values_to = 'fit') %>%
    hoist(sim, 'idx') %>%
    hoist(fit, 'cs') %>%
    select(-c(sim, fit)) %>%
    unnest_longer(cs) %>%
    unnest_wider(cs) %>%
    rowwise() %>%
    mutate(found = list(idx %in% cs), covered = any(found), n_found = sum(found)) %>%
    ungroup()

res_cs %>%
    filter(cs_size < 200) %>%
    {hist(.$cs_size)}

res_cs %>%
    filter(cs_size < 200) %>%
    group_by(method) %>%
    summarize(coverage = mean(covered), n_found = list(table(n_found)), mean_cs_size = mean(cs_size)) %>%
    unnest_wider(n_found)

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
    pivot_longer(c(ser, ser_abf, gibss, gibss_abf, linear_susie, rss), names_to = 'method', values_to = 'fit') %>%
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


# source('workflow/scripts/plotting/color_key.R')
source('workflow/scripts/plotting/plot_functions.R')
# colors <- list(ser = 'red', ses_abf = 'pink', gibss = 'blue', gibss_abf = 'green')
methods <- unique(res_pip$method)

# fdp vs power
fdp_plot <- res_pip %>%
    #dplyr::filter(method %in% methods) %>%
    unnest_longer(c(pip, causal)) %>%
    ggplot(aes(pip, causal, color=method)) + 
    stat_fdp_power(geom='path', max_fdp=1.0) +
    labs(x = 'FDP', y = 'Power')
fdp_plot


# PIP Calibration
pip_calibration_plot <- res_pip %>%
    unnest_longer(c(pip, causal)) %>%
    ggplot(aes(pip, causal)) +
    stat_pip_calibration() +
    geom_abline(intercept = 0, slope = 1, color='red') +
    theme_bw() +
    labs(x='PIP', y='Frequency') + 
    facet_wrap(vars(method))
pip_calibration_plot
