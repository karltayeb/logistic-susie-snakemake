import pickle
import numpy as np
from susiepy.generalized_ibss import generalized_ibss, fit_ser

factor = snakemake.wildcards.factor
sim_id = int(snakemake.wildcards.sim_id)
data = pickle.load(open('results/pdac/pdac_data_py.pkl', 'rb'))

genes = data['genes']
genesets = data['genesets']
X = data['X']
dat = data['factors'][factor]
y = dat['Ysim'][sim_id]

# fit sim data
sim_fit = generalized_ibss(X, y, L=20, maxit=50)
pickle.dump(sim_fit, open(snakemake.output.simfit, 'wb'))