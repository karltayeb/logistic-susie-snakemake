source(".Rprofile")

# load data
sumstats <- read.csv(snakemake@input$sumstat_path, sep = "\t")
ld <- read.table(snakemake@input$ld_path, sep = "\t")

# logistic has column Z_STAT, linear has column T_STAT. Either way, treat as Z-score.
if ("Z_STAT" %in% colnames(sumstats)) {
  z <- sumstats$Z_STAT
} else {
  z <- sumstats$T_STAT
}
n <- sumstats$OBS_CT[1]
rssfit <- susieR::susie_rss(z, as.matrix(ld), n)
saveRDS(rssfit, snakemake@output$rss)
