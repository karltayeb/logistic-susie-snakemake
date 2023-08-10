import pickle
import pandas as pd

fits = pickle.load(open('results/ukbb_geno/chr1_100794065_101983457/ser_sims_ser_fit.pkl', 'rb'))

dfs = []

for i in range(len(fits)):
    fit = fits[i]
    pip1 = fit['ser_fit']['alpha'] 
    pip2 = fit['ser_fit']['alpha_abf']
    causal = np.zeros(pip2.size)
    causal[fit['idx']] = 1.

    # add other simulation details here if necessary
    dfs.append(pd.DataFrame(dict(sim = i, method = 'laplace_abf', pip = pip1, causal = causal)))
    dfs.append(pd.DataFrame(dict(sim = i, method = 'abf', pip = pip2, causal = causal)))
    

res = pd.concat(dfs)
pickle.dump(res, 'results/ukbb_geno/chr1_100794065_101983457/ser_sim_ser_fit_pip_summary.pkl')