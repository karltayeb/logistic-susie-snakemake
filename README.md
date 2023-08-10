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
```

```{r}
devtools::install_github("kartayeb/gseasusie@cleanup")
```
