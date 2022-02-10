from cmath import nan
import math
import sys
import warnings
from sklearn.cluster import AffinityPropagation
from sklearn.exceptions import ConvergenceWarning
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
import random

def single_affinity_run(x, random_state=42, damping=0.5):
	af = AffinityPropagation(random_state=random_state, damping=damping).fit(x)
	predicted = af.predict(x)
	centers = len(af.cluster_centers_indices_) 
	if centers > 2 and centers < len(x):
		sc_score = metrics.silhouette_score(x, af.labels_, metric='sqeuclidean')
	else:
		sc_score = math.nan
	return (predicted, sc_score)

def cluster_by_haplotype (gr, dampings=[0.5,0.6,0.7,0.8,0.9], iterations=100, seed=42):
	with warnings.catch_warnings():
		warnings.filterwarnings("ignore", category=ConvergenceWarning)
		t_df = gr.as_df().set_index(['Chromosome', 'Start', 'End']).T
		# print(t_df.columns)
		varieties = t_df.index
		x = StandardScaler(with_mean=True).fit_transform(t_df)
		predicted = list()
		scores = list()
		used_dampings = list()
		random_states = list()
		random.seed(seed)
		max_seed = 2**32 - 1
		for dmp in dampings:
			for i in range(iterations):
				random_state = random.randint(1, max_seed)
				X_predicted, sc_score = single_affinity_run(x, random_state=random_state, damping=dmp)
				predicted.append(X_predicted)
				scores.append(sc_score)
				used_dampings.append(dmp)
				random_states.append(random_state)		
	return predicted, scores, used_dampings, random_states, varieties