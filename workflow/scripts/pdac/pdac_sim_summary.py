import glob

import pickle
import numpy as np
import pandas as pd
import time

def load_pickle(path):
    return pickle.load(open(path, 'rb'))

def save_pickle(obj, path):
    return pickle.dump(obj, open(path, 'wb'))

def extract_filename(path):
    return path.split('/')[-1].split('.')[0]

def pluck(d, query):
    return {k: d[k][query] for k in d.keys()}

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
    if len(alpha.shape) == 1:
        cs = _compute_cs(alpha, target_coverage)
    else:
        cs = {f'L{l}': _compute_cs(a) for l, a in enumerate(alpha)}
    return cs


def _annotate_cs(cs, idx=None, X=None):
    if idx is not None:
        cs['causal_in_cs'] = np.intersect1d(cs['cs'], idx)
        cs['n_causal_in_cs'] = cs['causal_in_cs'].size
    if X is not None:
        vars = cs['cs']
        if(vars.size > 200):
            # save some time if cs is too big
            vars = np.random.choice(vars, 200, replace=False)
        cs['purity'] = np.abs(np.min(np.corrcoef(X[:, vars].T)))
    return cs

def annotate_cs(cs, **kwargs):
    if 'cs' in cs.keys():
        cs =  _annotate_cs(cs, **kwargs)
    else:
        cs = {k: _annotate_cs(v, **kwargs) for k, v in cs.items()}
    return cs

def compute_pips(alpha):
    if alpha.shape[0] == 0:
        return np.zeros(alpha.shape[1])
    else:
        log_alpha = np.log(1 - alpha + 1e-12)
        return 1 - np.exp(log_alpha.sum(0))


data = load_pickle('results/pdac/pdac_data_py.pkl')
X = data['X']

fit_paths = glob.glob('results/pdac/fits/*.pkl')
fits = {extract_filename(p)[:-6]: load_pickle(p) for p in fit_paths}

tic = time.perf_counter()
sim_paths = glob.glob('results/pdac/sims/*.pkl')
sims = {extract_filename(p): load_pickle(p) for p in sim_paths}
toc = time.perf_counter()
elapsed = toc - tic
print(f'Time to load simulations = {elapsed}')

def summarize_lasso_sim(factor):
    print(factor)
    
    beta0 = data['factors'][factor]['coef'][0]
    beta = data['factors'][factor]['coef'][1:]
    gs_size = X.sum(0)
    idx = np.where(beta != 0 )[0]

    discovered = 0
    rows = []
    for i in range(20):
        fit = sims.get(f'{factor}_{i}_gibss')
        if fit is not None:
            Alpha = fit['alpha']
            cs = annotate_cs(compute_cs(Alpha), idx=idx, X=X)
            for l in cs.keys():
                if cs[l]['cs'].size < 1000:
                    discovered += np.isin(idx, cs[l]['causal_in_cs']).astype(int)
                    rows.append(dict(
                        factor = factor,
                        sim_id = i, 
                        ser = l,
                        cs = cs[l]['cs'],
                        cs_alpha = cs[l]['alpha'],
                        cs_size = cs[l]['cs_size'],
                        cs_causal = cs[l]['causal_in_cs'],
                        cs_n_causal = cs[l]['n_causal_in_cs'],
                        cs_purity = cs[l]['purity'],
                        ser_lbf = fit['lbf_ser'][int(l[1:])],
                        ser_prior_variance = fit['prior_variance'][int(l[1:])] 
                    ))
    
    print(discovered)
    ser_summary = pd.DataFrame(rows)
    return ser_summary


def make_effect_summary(factor):
    beta0 = data['factors'][factor]['coef'][0]
    beta = data['factors'][factor]['coef'][1:]
    gs_size = X.sum(0)
    idx = np.where(beta != 0 )[0]
    
    effect_summary = pd.DataFrame(dict(
        factor = factor,
        idx = idx, 
        beta = beta[idx], 
        gs_size = gs_size[idx]
    ))
    return effect_summary

def make_pip_tbl(factor):
    causal = np.array(data['factors'][factor]['coef'][1:] != 0).astype(int)
    sim_ids = [s for s in sims.keys() if factor in s]
    res = []
    for sim in sim_ids:
        cs = compute_cs(sims[sim]['alpha'])
        cs_size = pluck(cs, 'cs_size')
        include = np.array([cs_size[f'L{l}'] < 100 for l in range(len(cs_size))])
        pip = compute_pips(sims[sim]['alpha'][include])
        res.append(dict(factor = factor, sim_id = int(sim.split('_')[1]), pip = pip, causal = causal))
    res = pd.DataFrame(res)  
    return res

factors = np.unique([x.split('_')[0] for x in list(sims.keys())])
ser_summary = pd.concat([summarize_lasso_sim(factor) for factor in factors])
pip_tbl = pd.concat([make_pip_tbl(factor) for factor in factors])

effect_summary = pd.concat([make_effect_summary(factor) for factor in factors])
discovered = ser_summary[['factor', 'sim_id', 'cs_causal']]\
    .explode('cs_causal')\
    .drop_duplicates()\
    .groupby(['factor', 'cs_causal'])\
    .size()\
    .reset_index()\
    .rename(columns={'cs_causal': 'idx', 0: 'n_discovered'})\
    .set_index(['factor', 'idx'])
effect_summary = effect_summary\
    .set_index(['factor','idx'])\
    .join(discovered, how='left')\
    .reset_index()\
    .fillna(0)
    
save_pickle(
    dict(sim_summary_ser = ser_summary,
         sim_summary_effect = effect_summary,
         pips = pip_tbl),
    'results/pdac/sim_summary.pkl'
)