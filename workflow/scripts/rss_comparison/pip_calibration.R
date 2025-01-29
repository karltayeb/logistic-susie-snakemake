library(dplyr)
library(ggplot2)
library(tidyr)
source("workflow/scripts/plotting/plot_functions.R")

sims <- purrr::map(1:3, ~ readRDS(glue::glue("results/rss_comparison/simulation/sim{.x}/sim{.x}.sim.rds")))
fits <- purrr::map(1:3, ~ readRDS(glue::glue("results/rss_comparison/simulation/sim{.x}/sim{.x}.logistic_rss_fit.rds")))

compute_pip <- function(alpha) {
  apply(alpha, 2, function(x) {
    1 - prod(1 - x)
  })
}

make_causal_vector <- function(sim) {
  # causal <- rep(0, sim$n_snps)
  causal <- rep(0, 2582)
  causal[sim$causal_idx] <- 1
  return(causal)
}

sims <- tibble::tibble(fit = fits, sim = sims)
pips <- sims %>%
  rowwise() %>%
  mutate(
    pip = list(compute_pip(fit$alpha)),
    causal = list(make_causal_vector(sim))
  )
plot <- pips %>% make_pip_calibration_plot()
ggsave("results/rss_comparison/simulation/pip_calibration.pdf", plot)

cs_coverage_plot <- function(cs){

}

cs_coverage_plot <- function(cs, max_cs_size=1e10){
  plot <- cs %>% 
    flatten_cs() %>%
    filter(cs_size <= max_cs_size) %>%
    hoist(sim, 'L') %>%
    ggplot(aes(x=as.factor(L), y=as.numeric(covered), color=fit_method)) + 
    stat_bs_mean(position=position_dodge(0.1)) +
    scale_color_manual(values=colors) + 
    geom_hline(yintercept = 0.95, color='red') +
    theme_bw()
  
  return(plot)
}
sims2 <- sims %>%
  rowwise() %>%
  mutate(cs = list(fit$sets$cs), causal_idx = list(sim$causal_idx)) %>%
  select(c(sim, cs, causal_idx)) %>%
  unnest_longer(cs) %>%
  rowwise() %>%
  mutate(covered = any(causal_idx %in% cs), size=length(cs), n_causal_variants = sim$n_causal_variants) 
  
plot <- sims2 %>%
  ggplot(aes(x=n_causal_variants, y=as.numeric(covered))) +
  stat_bs_mean() + theme_bw()
ggsave('results/rss_comparison/simulation/cs_coverage.pdf', plot)
