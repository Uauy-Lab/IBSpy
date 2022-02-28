from cmath import nan
import math
import statistics
import sys
from numpy import NaN
import pandas as pd
import warnings
from sklearn.cluster import AffinityPropagation
from sklearn.exceptions import ConvergenceWarning
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
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

	@property
	def failed(self):
		return math.isnan(self.score)

	def as_df(self):

		df = pd.DataFrame({
			'Chromosome': self.chromosome,
			'Start': self.start, 
			'End': self.end, 
			'variety': self.varieties, 
			'group':self.predicted,
			'mutual_info_score':self.mutual_info_score,
			'sc_score':self.score,
			'number_of_runs':self.number_of_runs,
			'damping':self.damping
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

def cluster_by_haplotype (gr, dampings=[0.5,0.6,0.7,0.8,0.9], iterations=100, seed=42, max_missing=5, min_iterations=5):
	with warnings.catch_warnings():
		warnings.filterwarnings("ignore", category=ConvergenceWarning)
		t_df = gr.as_df().set_index(['Chromosome', 'Start', 'End']).T
		# print(t_df.columns)
		varieties = t_df.index
		x = StandardScaler(with_mean=True).fit_transform(t_df)
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