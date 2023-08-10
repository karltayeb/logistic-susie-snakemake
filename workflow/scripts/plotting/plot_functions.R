
aes = ggplot2::aes

# Bootstrap mean + CI of y ------
bootstrap_mean_ci <- function(y, reps=50){
  res <- quantile(purrr::map_dbl(1:reps, ~mean(sample(y, replace = T))), c(0.05, 0.5, 0.95))
  names(res) <- c('ymin', 'y', 'ymax')
  return(as.list(res))
}

StatBootstrapMean <- ggplot2::ggproto(
  "StatBootstrapMean", ggplot2::Stat,
  compute_group = function(data, scales){
    # x = pip, y = causal
    res <- data %>%
      group_by(x) %>%
      summarise(
        size = n(),
        y_range = list(bootstrap_mean_ci(y))
      ) %>% 
      unnest_wider(y_range) %>%
      ungroup()
    with(res, data.frame(x = x, y=y, ymin=ymin, ymax=ymax))
  },
  required_aes = c("x", "y"),
  default_aes = aes(ymin =stat(ymin), ymax=stat(ymax))
)

stat_bs_mean <- function(mapping = NULL, data = NULL, geom = "pointrange",
                         position = "identity", na.rm = FALSE, show.legend = NA,
                         inherit.aes = TRUE, ...){
  ggplot2::layer(
    stat = StatBootstrapMean, data=data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

# PIP Calibration ---------
StatPIPCal <- ggplot2::ggproto(
  "StatPIPCal", ggplot2::Stat,
  compute_group = function(data, scales){
    pip_bins <- seq(0, 1.0, by=0.05)

    # x = pip, y = causal
    res <- data %>%
      mutate(pip_bin = cut(x, pip_bins, labels=F)) %>%
      mutate(pip_bin = seq(0.05, 1.0, by=0.05)[pip_bin]) %>%
      group_by(pip_bin) %>%
      summarise(
        size = n(),
        y_range = list(bootstrap_mean_ci(y))
      ) %>% 
      unnest_wider(y_range) %>%
      ungroup()
    with(res, data.frame(x = pip_bin, y=y, ymin=ymin, ymax=ymax))
  },
  required_aes = c("x", "y"),
  default_aes = aes(ymin =stat(ymin), ymax=stat(ymax))
)

stat_pip_calibration <- function(mapping = NULL, data = NULL, geom = "pointrange",
                     position = "identity", na.rm = FALSE, show.legend = NA,
                     inherit.aes = TRUE, ...){
  ggplot2::layer(
    stat = StatPIPCal, data=data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

make_pip_calibration_plot <- function(pips){
  pips %>%
  unnest_longer(c(pip, causal)) %>%
  ggplot(aes(pip, causal)) + 
  stat_pip_calibration() +
  geom_abline(intercept = 0, slope = 1, color='red') +
  theme_bw() +
  labs(x='PIP', y='Frequency')
}

# PIPs Power vs FDR --------

# take a list of pips and their causal status
# return power, FDP, and threshold
compute_power_fdr_thresh <- function(pip, causal){
  n_causal <- sum(causal)

  sorter <- order(pip, decreasing = T)
  sorted_pips <- pip[sorter]
  sorted_causal <- causal[sorter]
  n <- 1:length(sorted_causal)

  power <- cumsum(sorted_causal) / n_causal
  fdr <- (n - cumsum(sorted_causal)) / n
  thresh <- sorted_pips

  # prune
  include <- logical(length(thresh))
  include[1] <- T
  last <- 1
  for(i in 2:length(thresh)){
    if(((fdr[i] - fdr[last]) > 0.005) | ((power[i] - power[last]) > 0.005)){
      include[i] <- T
      last <- i
    }
  }

  power <- power[include]
  fdr <- fdr[include]
  thresh <- thresh[include]
  return(list(power=power, fdr=fdr, thresh=thresh))
}

StatFDPvP <- ggplot2::ggproto(
  "StatFDP", ggplot2::Stat,
  compute_group = function(data, scales, max_fdp=1, max_power=1){
    # x = pip, y = causal
    res <- compute_power_fdr_thresh(data$x, data$y)
    res <- data.frame(x = res$fdr, y = res$power) %>%
      filter(x <= max_fdp)
    return(res)
  },
  required_aes = c('x', 'y')
)

stat_fdp_power <- function(mapping = NULL, data = NULL, geom = "path",
                     position = "identity", na.rm = FALSE, show.legend = NA,
                     inherit.aes = TRUE, ...){
  ggplot2::layer(
    stat = StatFDPvP, data=data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

# takes an unnested tbl of pips (1 row, 1 pip) and produces the fdr plot
make_fdp_plot <- function(pips, methods, colors, ...){
  colors <- colors[methods]
  pips %>%
    unnest_longer(c(pip, causal)) %>%
    ggplot(aes(pip, causal, color=fit_method)) + 
    stat_fdp_power(geom='path', ...) +
    labs(x='False Discovery Proportion', y='Power') +
    scale_color_manual(values=colors)
}

# CS Coverage ---------
flatten_cs <- function(cs){
  # 1 row, 1 credible set
  cs <- cs %>%
    unnest_longer(c(cs)) %>%
    rowwise() %>%
    mutate(covered = any(cs$covered), cs_size=cs$size, which_covered = list(cs$which_covered))
  return(cs)
}

make_cs_coverage_plot <- function(cs, max_cs_size=1e10){
  methods <- unique(cs$fit_method)
  colors <- plot_colors()[methods]
  
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

# CS Size ---------
#' plot cs coverage by cs size
make_coverage_by_cs_size <- function(cs){
  methods <- unique(cs$fit_method)
  colors <- plot_colors()[methods]

  cs %>%
    flatten_cs() %>%
    mutate(
      covered = as.numeric(covered),
      size_bin = cut(cs_size, c(0, 5, 10, 25, 50, 100))
    ) %>%
    ggplot(aes(x=size_bin, y=covered, color=fit_method)) +
    stat_bs_mean(position=position_dodge(0.3)) +
    geom_hline(yintercept=0.95) +
    scale_color_manual(values=colors) + 
    theme_bw()
}

#' plot cs size on histogram
make_cs_size_histogram <- function(cs){
  methods <- unique(cs$fit_method)
  colors <- plot_colors()[methods]

  cs %>%
    flatten_cs() %>%
    ggplot(aes(x=cs_size, fill=fit_method)) +
    geom_histogram() + 
    facet_wrap(vars(fit_method)) +
    scale_fill_manual(values=colors) + 
    theme_bw()
}

scatter2 <- function(data, mapping, ...){
  ggplot(data=data, mapping = mapping) +
  geom_point() +
  geom_abline(yintercept=0, slope=1, color='red') +
  theme_bw()
}

bin2 <- function(data, mapping, ...){
  ggplot(data=data, mapping = mapping) +
  geom_bin2d() +
  geom_abline(yintercept=0, slope=1, color='red') +
  theme_bw()
}

#' match CSs across simulations and plot between pairs of methods
make_pairwise_cs_size_plot <- function(cs){
  # directly plot CSs against eachother
  pairs <- cs %>%
    flatten_cs() %>%
    unnest_longer(which_covered) %>%
    group_by(fit_method, hash, which_covered) %>%
    summarize(cs_size = min(cs_size)) %>%
    ungroup() %>%
    pivot_wider(names_from=fit_method, values_from=cs_size) %>%
    select(-c(hash, which_covered))

  pairs %>%
    GGally::ggpairs(lower = list(continuous = bin2))
}

# Helpers --------
# helper functions to make filtring down to different methods less annoying
filter_susie <- function(tbl){
  tbl %>%
    filter(stringr::str_detect(fit_method, 'susie|ibss'))
}

filter_ser <- function(tbl){
  tbl %>%
    filter(stringr::str_detect(fit_method, 'ser'))
}

