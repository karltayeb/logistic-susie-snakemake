import pickle
import numpy as np
import pandas as pd
from susiepy.logistic_susie import logistic_gibss

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
    print(f'fitting SuSiE (L = {L}) using Laplace ABF...')
    n_chunks = int(X.shape[1] / 10)
    logistic_gibss_fit = logistic_gibss(X.toarray(), np.array(sim['y']), L = L, n_chunks = n_chunks, **kwargs)
    logistic_gibss_fit['cs'] = compute_cs(logistic_gibss_fit['alpha'])
    print("elapsed time = {elapsed_time}".format(**logistic_gibss_fit))
    return logistic_gibss_fit

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
pickle.dump(fit, open(f'results/ukbb_geno/{region}/{sim_id}/gibss.pkl', 'wb'))