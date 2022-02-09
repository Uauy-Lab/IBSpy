import sys
from sklearn.cluster import AffinityPropagation
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
import random


def single_affinity_run(x, random_state=42, damping=0.5):
	af = AffinityPropagation(random_state=random_state, damping=damping).fit(x)
	predicted = af.predict(x)
	if len(af.cluster_centers_indices_) > 2:
		sc_score = '%0.2f'%metrics.silhouette_score(x, af.labels_, metric='sqeuclidean')
	else:
		sc_score = 0
	return (predicted, sc_score)



def cluster_by_haplotype (gr, dampings=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9], iterations=100, seed=42):
	t_df = gr.as_df.set_index(['Chromosome', 'Start', 'End']).T
	x = StandardScaler(with_mean=True).fit_transform(t_df)
	predicted = list()
	scores = list()
	used_dampings = list()
	random_states = list()
	random.seed(seed)
	for i in range(iterations):
		for dmp in dampings:
			random_state = random.randint(1, sys.maxsize)
			X_predicted, sc_score = single_affinity_run(x, random_state=random_state, damping=dmp)
			predicted.append(X_predicted)
			scores.append(sc_score)
			used_dampings.append(dmp)
			random_states.append(random_state)
	return predicted, scores, used_dampings, random_states