## R package notes

I found that the midway installation of R was having trouble finding important headers when in the conda environment.
So we are forced to use an install of R within conda. 
In this case, it is best practice to install as many packages as we can via conda. 
I try to only install from github repos within R (although, this also triggers install of any dependencies that are not already installed via conda).

```
mamba install -c conda-forge r-base r-essentials r-devtools r-reticulate
mamba install -c r  r-tidyverse
mamba install -c bioconda r-bioconductor bioconductor-annotationsdbi bioconductor-org.hs.eg.db 
mamba install -c bioconda bioconductor-org.hs.eg.db
mamba install -c conda-forge r-susier r-tictoc
```

```{r}
devtools::install_github("kartayeb/gseasusie@cleanup")
```


### RSS Comparison

It's common practice to apply summary statistics based fine-mapping methods to summary statistics from GWAS using logistic regression, GLMs, generalized linear mixed models, etc.
Most summary statistic based methods exploit the fact that in a Gaussian linear model the effect estimates of the joint model can be recovered from the marginal summary statistics and the in-sample LD matrix.
In other regression models, it is implicity assumed that these relationship at least hold approximately.
Here we explore the behavior of SuSiE-RSS on summary statistics obtained from logistic regression.

