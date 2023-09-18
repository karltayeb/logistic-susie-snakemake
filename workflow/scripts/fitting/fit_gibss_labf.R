# Fit Logistic GIBSS-LABF
library(reticulate)
np <- import("numpy")
logistic_susie <- import("susiepy.logistic_susie")

augment <- function(X, y){
    n <- nrow(X)
    p <- ncol(X)
    yaug <- c(y, c(1, 1, 0, 0))
    tmp <- do.call(cbind, purrr::map(1:p, ~c(1, 0, 1, 0)))
    Xaug <- rbind(X, tmp)
    return(list(X=Xaug, y=yaug))
}

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

fit_ser <- function(X, y, ...){
    fit <- logistic_susie$logistic_ser(X, np$array(y), ...)
    fit_np <- np_cast(fit)
    fit_np$cs <- logisticsusie::compute_cs(as.vector(fit_np$alpha))
    return(fit_np)
}


fit_gibss <- function(X, y, L, ...){
    fit <- logistic_susie$logistic_gibss(X, np$array(y), as.integer(L), ...)
    fit$ser_fits <- NULL
    fit_np <- np_cast(fit)
    #fit_np$cs <- logisticsusie::compute_cs(fit_np$alpha)
    return(fit_np)
}

X <- readRDS(snakemake@input$X)
sim <- readRDS(snakemake@input$sim)
params <- snakemake@params$fit_params

if('maxit' %in% names(params)){
    params$maxit <- as.integer(params$maxit)
}

data <- augment(X, sim$y)
#params = list(L = 5)
print(params)
fit_args <- c(data, params)
fit <- rlang::exec(fit_gibss, !!!fit_args)
saveRDS(fit, snakemake@output$fit)