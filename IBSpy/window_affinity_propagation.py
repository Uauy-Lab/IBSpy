from audioop import bias
from cmath import nan
import math
import statistics
import sys
from numpy import NaN, ndarray
import pandas as pd
import numpy as np
import warnings
from sklearn.cluster import AffinityPropagation
from sklearn.exceptions import ConvergenceWarning
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
import scipy.stats as stats
from scipy.stats import kurtosis
from scipy.stats import skew
import random


class AffinityRunResults:
	def __init__(self, predicted, score, damping, random_state, varieties) -> None:
		self.predicted = predicted
		self.score = score
		self.damping = damping
		self.random_state = random_state
		self.varieties  = varieties 
		self.mutual_info_score = 0
		self.number_of_runs = 0
		self.stdev = 0
		self.chromosome = None
		self.start = None
		self.end = None
		self.descriptive_stats = {}
		self.norm_descriptive_stats = {}
		self.samples = 0
		self.windows = 0
		

	@property
	def failed(self):
		return math.isnan(self.score)

	@property
	def n(self):
		return self.samples * self.windows

	@property
	def n_haplotypes(self):
		return len(set(self.predicted))

	def as_df(self):

		df = pd.DataFrame({
			'Chromosome': self.chromosome,
			'Start': self.start, 
			'End': self.end, 
			'genotype': self.varieties, #TODO Maybe change this to genotype. 
			'group':self.predicted,
			'n_haplotypes': self.n_haplotypes,
			'mutual_info_score':self.mutual_info_score,
			'sc_score':self.score,
			'number_of_runs':self.number_of_runs,
			'damping':self.damping,
			'n_windows': self.windows,
			'n_samples': self.samples,
			'n': self.n, 
			'median': self.descriptive_stats["median"],
			'mean': self.descriptive_stats["mean"],
			'stdev': self.descriptive_stats["stdev"],
			'skew': self.descriptive_stats["skew"],
			'kurt': self.descriptive_stats["kurt"],
			'median_scaled': self.norm_descriptive_stats["median"],
			'mean_scaled': self.norm_descriptive_stats["mean"],
			'stdev_scaled': self.norm_descriptive_stats["stdev"],
			'skew_scaled': self.norm_descriptive_stats["skew"],
			'kurt_scaled': self.norm_descriptive_stats["kurt"]
			}
			)
		# print(df)
		return df

def select_best_cluster(affinity_runs ):
	best_score = 0
	best_cluster = None
	for i, run_1 in enumerate(affinity_runs):
		p = run_1.predicted
		scores =  list()
		for j, run_2, in enumerate(affinity_runs):
			if j == i:
				next
			info_score = metrics.adjusted_mutual_info_score(p, run_2.predicted)
			scores.append(info_score)
		mean_score = statistics.mean(scores)
		run_1.mutual_info_score = mean_score
		run_1.number_of_runs    = len(affinity_runs)
		if(len(scores) > 2):
			run_1.stdev = statistics.stdev(scores)
		else:
			run_1.stdev = NaN
		if mean_score > best_score:
			best_score   = mean_score
			best_cluster = run_1
	return best_cluster

def single_affinity_run(x, random_state=42, damping=0.5):
	af = AffinityPropagation(random_state=random_state, damping=damping).fit(x)
	predicted = af.predict(x)
	centers = len(af.cluster_centers_indices_) 
	if centers > 2 and centers < len(x):
		sc_score = metrics.silhouette_score(x, af.labels_, metric='sqeuclidean')
	else:
		sc_score = math.nan
	return (predicted, sc_score)

def descriptive_stats(gr: pd.DataFrame):
	ret = {}
	melted  = gr.melt(id_vars=["Chromosome", "Start", "End"])
	ret['stdev'] = statistics.stdev(melted['value'])
	ret['median'] = statistics.median(melted['value'])
	ret['mean'] = statistics.mean(melted['value'])
	ret['skew'] = skew(melted['value'])
	ret['kurt'] = kurtosis(melted['value'], bias=True)
	#print(ret)
	return ret 


def descriptive_stats_normalized(x: ndarray):
	ret = {}
	flat = x.flatten()
	ret['stdev'] = statistics.stdev(flat)
	ret['median'] = statistics.median(flat)
	ret['mean'] = statistics.mean(flat)
	ret['skew'] = skew(flat)
	ret['kurt'] = kurtosis(flat, bias = True)
	return ret 


def cluster_by_haplotype (gr, dampings=[0.5,0.6,0.7,0.8,0.9], iterations=1000, seed=42, max_missing=5, min_iterations=10):
	with warnings.catch_warnings():
		warnings.filterwarnings("ignore", category=ConvergenceWarning)

		if not isinstance(gr, pd.DataFrame):
			gr = gr.as_df()
		n_windows = gr.shape[0]
		n_samples = gr.shape[1] - 3 #Minus three as the columns Chromosomo, Start and End must not be called. 
		t_df = gr.set_index(['Chromosome', 'Start', 'End']).T
		varieties = t_df.index
		x:ndarray = StandardScaler(with_mean=True).fit_transform(t_df)

		runs = list()
		best_dmp_score = 0
		best_dmp = dampings[0]
		random.seed(seed)
		for dmp in dampings:			
			tmp_score = 0
			tmp_runs = multiple_affinity_run(min_iterations, max_missing, min_iterations, varieties, x, dmp )
			if len(tmp_runs) > 0:
				tmp_score = statistics.mean(map(lambda r:r.score, tmp_runs))
			if tmp_score > best_dmp_score:
				best_dmp = dmp
				runs = tmp_runs
				best_dmp_score = tmp_score
		best_cluster = select_best_cluster(runs)

		if best_cluster is not None and  best_cluster.stdev > 0.001:
			tmp_runs = multiple_affinity_run(iterations - len(runs), max_missing, min_iterations, varieties, x, best_dmp )
			runs.extend(tmp_runs)
		desc_stats = descriptive_stats(gr)
		desc_stats_norm = descriptive_stats_normalized(x)
		for run in runs:
			run.descriptive_stats = desc_stats
			run.norm_descriptive_stats = desc_stats_norm
			run.windows = n_windows
			run.samples = n_samples
			#TODO: Add number of haplotypes in window. 
	return runs

def multiple_affinity_run(iterations, max_missing, min_iterations, varieties, x, dmp):
	tmp_runs = list()
	missing = 0
	max_seed = 2**32 - 1
	for i in range(iterations):
		random_state = random.randint(1, max_seed)
		x_predicted, sc_score = single_affinity_run(x, random_state=random_state, damping=dmp)
		tmp = AffinityRunResults(x_predicted, sc_score, dmp, random_state, varieties)
		if tmp.failed:
			missing += 1
		else:
			tmp_runs.append(tmp)
		if missing > max_missing and i < min_iterations:
			tmp_runs.clear()
			break
	return tmp_runs