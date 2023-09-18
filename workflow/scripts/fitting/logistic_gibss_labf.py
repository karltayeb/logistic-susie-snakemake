import numpy as np
from scipy.stats import norm
from scipy.special import logsumexp
from susiepy.logistic_susie import logistic_ser as fit_ser
from susiepy.generalized_ibss import gibss_generator

def compute_alpha_abf(ser_fit):
    """
    Modify SER so it uses ABF instead of Laplace ABF
    """
    betahat = ser_fit['betahat']
    shat2 = ser_fit['shat2']
    prior_variance = ser_fit['prior_variance']
    log_abf = norm.logpdf(betahat, 0, scale = np.sqrt(shat2 + prior_variance)) \
        - norm.logpdf(betahat, 0, scale = np.sqrt(shat2))
    alpha = np.exp(log_abf - logsumexp(log_abf))
    
    ser_fit['lbf'] = log_abf
    ser_fit['alpha'] = alpha
    return ser_fit

def fit_ser_abf(X, y, *args, **kwargs):
    """
    Fit SER using ABF rather than Laplace ABF
    """
    ser_fit = fit_ser(X, y, *args, **kwargs)
    ser_fit = compute_alpha_abf(ser_fit)
    ser_fit['psi'] = X @ (ser_fit['post_mean'] * ser_fit['alpha']) + np.inner(ser_fit['intercept'], ser_fit['alpha'])
    return ser_fit

logistic_gibss_abf = gibss_generator(fit_ser_abf)
