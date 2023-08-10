import pickle
import numpy as np
from susiepy.generalized_ibss import generalized_ibss

factor = snakemake.wildcards.factor
data = pickle.load(open('results/pdac/pdac_data_py.pkl', 'rb'))

genes = data['genes']
genesets = data['genesets']
X = data['X']
dat = data['factors'][factor]

# fit original data
fit = generalized_ibss(X, dat['y'], L=20, estimate_prior_variance = True, maxit=50)
pickle.dump(fit, open(snakemake.output.fit, 'wb'))

# fit sim data
# M = data['Ysim'].shape[0]
# sim_fits = {m : generalized_ibss(X, dat['Ysim'][m], L=20) for m in range(M)}
# pickle.dump(fit, open(snakemake.output.sim_fits, 'wb'))