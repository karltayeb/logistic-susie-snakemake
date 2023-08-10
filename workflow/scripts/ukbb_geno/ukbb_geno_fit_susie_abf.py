import pickle
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.special import logsumexp
import susiepy
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

def _compute_cs(alpha, target_coverage = 0.95):
    sorter = np.argsort(alpha)[::-1]
    cumalpha = np.cumsum(alpha[sorter])
    cs_size = np.argmax(cumalpha >= 0.95) + 1

    cs = dict(cs_size = cs_size, 
            cs = sorter[:cs_size], 
            target_coverage = target_coverage, 
            coverage = cumalpha[cs_size - 1],
            alpha = alpha[sorter[:cs_size]])
    return cs

def compute_cs(alpha, target_coverage = 0.95):
    """
    compute cs
    """
    if len(alpha.shape) == 1:
        cs = _compute_cs(alpha, target_coverage)
    else:
        cs = {f'L{l}': _compute_cs(a) for l, a in enumerate(alpha)}
    return cs

def fit_sim(sim, L, **kwargs):
    """
    fit simulation, return dict with simulation and fit
    """
    print(f'fitting SuSiE (L = {L}) using ABF...')
    n_chunks = int(X.shape[1] / 10)
    logistic_gibss_abf_fit = logistic_gibss_abf(X.toarray(), np.array(sim['y']), L = L, n_chunks = n_chunks, **kwargs)
    logistic_gibss_abf_fit['cs'] = compute_cs(logistic_gibss_abf_fit['alpha'])
    print("elapsed time = {elapsed_time}".format(**logistic_gibss_abf_fit))
    return(logistic_gibss_abf_fit)

# unpack wildcards
region = snakemake.wildcards.region
sim_id = snakemake.wildcards.sim_id


# load X
data = pickle.load(open(f'results/ukbb_geno/{region}/genotype.pkl', 'rb'))
n, p = data['X'].shape
X = data['X']


# load sims
manifest = pd.read_csv('config/ukbb_sim/ukbb_sim_manifest.tsv', sep='\t')
sim_path = f'results/ukbb_geno/{region}/sim.pkl'
sims = pickle.load(open(sim_path, 'rb'))

fit = fit_sim(sims[sim_id], L=5, maxit=50, tol=1e-10)
pickle.dump(fit, open(f'results/ukbb_geno/{region}/{sim_id}/gibss_abf.pkl', 'wb'))

# fit model, save each loop so we can use partially completed output
# ser_fits = {}
# for i in sims.keys():
#     print(i)
#     ser_fits[i] = fit_sim(sims[i], L=5, maxit=20) 
#     pickle.dump(ser_fits, open(f'results/ukbb_geno/{region}/ser_sims_gibss_abf_fit.pkl', 'wb'))