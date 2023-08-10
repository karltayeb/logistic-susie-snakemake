import pickle
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.special import logsumexp
import susiepy
from susiepy.logistic_susie import logistic_ser as fit_ser

def compute_abf_alpha(ser_fit):
    betahat = ser_fit['betahat']
    shat2 = ser_fit['shat2']
    prior_variance = ser_fit['prior_variance']
    log_abf = norm.logpdf(betahat, 0, scale = np.sqrt(shat2 + prior_variance)) \
        - norm.logpdf(betahat, 0, scale = np.sqrt(shat2))
    alpha = np.exp(log_abf - logsumexp(log_abf))
    
    ser_fit['log_abf'] = log_abf
    ser_fit['alpha_abf'] = alpha
    return(ser_fit)


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


def fit_ser2(X, y, **kwargs):
    """
    fit ser, also record log_abf and approximation of alpha from log_abf
    """
    # fit model
    n_chunks = int(X.shape[1] / 10) # about right chunk size for 50k individuals
    ser_fit = fit_ser(X, y, n_chunks=n_chunks, **kwargs)
    
    # add ABF computations
    ser_fit = compute_abf_alpha(ser_fit)
    
    # compute cs
    ser_fit['cs'] = compute_cs(ser_fit['alpha'])
    ser_fit['cs_abf'] = compute_cs(ser_fit['alpha_abf'])
    return(ser_fit)


def fit_sim(sim):
    """
    fit simulation, return dict with simulation and fit
    """
    print('fitting...')
    ser_fit = fit_ser2(X.toarray(), np.array(sim['y']))
    return(ser_fit)


# region = 'chr1_100794065_101983457'
# sim_id = '59f82'

region = snakemake.wildcards['region']
sim_id = snakemake.wildcards['sim_id']

# load X
data = pickle.load(open(f'results/ukbb_geno/{region}/genotype.pkl', 'rb'))
n, p = data['X'].shape
X = data['X']


# load sims
manifest = pd.read_csv('config/ukbb_sim/ukbb_sim_manifest.tsv', sep='\t')
sim_path = f'results/ukbb_geno/{region}/sim.pkl'
sims = pickle.load(open(sim_path, 'rb'))

# fit model, save each loop so we can use partially completed output
print(sim_id)
ser_fit = fit_sim(sims[sim_id])
pickle.dump(ser_fit, open(f'results/ukbb_geno/{region}/{sim_id}/ser.pkl', 'wb'))

# ser_fits = {}
# print(sim_id)
# ser_fits[sim_id] = fit_sim(sims[sim_id])
# pickle.dump(ser_fits, open(f'results/ukbb_geno/{region}/{sim_id}/ser.pkl', 'wb'))

# ser_fits = {i: fit_sim(sims[i]) for i in sims.keys()}
# pickle.dump(ser_fits, open('results/ukbb_geno/chr1_100794065_101983457/ser_sims_ser_fit.pkl', 'wb'))