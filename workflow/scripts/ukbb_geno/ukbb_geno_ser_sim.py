import pickle
import numpy as np

# load data
# data = pickle.load(open(snakemake.input[0], 'rb'))
data = pickle.load(open('results/ukbb_geno/chr1_100794065_101983457/genotype.pkl', 'rb'))
n, p = data['X'].shape
X = data['X']

allele_freq = np.array(X.mean(0)).flatten() / 2
maf = np.minimum(allele_freq, 1 - allele_freq)

def make_sim(lower=0.2, upper=0.25, b0=-4, n_causal=1):
    # pick causal variant withing range
    idxs = np.intersect1d(np.where(maf > lower)[0], np.where(maf < upper)[0])
    idx = np.random.choice(idxs, n_causal, replace=False).flatten() # sample causal variants

    # simulate y
    x = (X[:, idx]).toarray()
    beta =  np.abs(np.random.normal(loc=0, scale=1, size=idx.size))

    # simulate fixed intercept
    beta0 = b0
    logit = beta0 + x @ beta

    # simulate sample specific intercept
    # beta0 = np.random.normal(loc = b0, scale=0, size=n)
    # logit = beta0 + x @ beta

    y = np.random.binomial(1, 1 / (1 + np.exp(-logit)), logit.size)
    sim = dict(
        idx = idx,
        beta = beta,
        beta0 = beta0,
        y = y,
        ybar = y.mean()
    )
    return sim

np.random.seed(42)
sims = {i: make_sim() for i in range(50)}
sim_path = 'results/ukbb_geno/chr1_100794065_101983457/ser_sims.pkl'
pickle.dump(sims, open(sim_path, 'wb'))