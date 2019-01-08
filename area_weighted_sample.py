import numpy as np
from scipy.stats import binom

avg_maldi = [0.420247211,
	0.046956096,
	0.036893222,
	0.031744461,
	0.044888963,
	0.419270046]

stdev_maldi = [0.036706515,
	0.027800178,
	0.020484287,
	0.025120514,
	0.030945762,
	0.048171729]

def get_random_maldi(avg_maldi,stdev_maldi):
	random_maldi = [np.random.normal(avg,stdev) for avg,stdev in zip(avg_maldi,stdev_maldi)]
	return random_maldi

def get_ssr(maldi,ddt_conc):
	binomial = [binom.pmf(i,5,ddt_conc) for i in range(6)]
	ssr = np.sum(np.square(np.array(maldi)-np.array(binomial)))
	return ssr

num_samples = 2000
ssrs = np.empty(num_samples)

for i in range(num_samples):
	maldi = get_random_maldi(avg_maldi,stdev_maldi)
	ssrs[i] = get_ssr(maldi,0.5)

print(np.mean(ssrs))
print(np.std(ssrs))
