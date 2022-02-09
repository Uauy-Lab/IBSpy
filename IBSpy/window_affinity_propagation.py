from sklearn.cluster import AffinityPropagation
from sklearn.preprocessing import StandardScaler
from sklearn import metrics


def single_affinity_run(x, random_state=42, damping=0.5):
	af = AffinityPropagation(random_state=random_state, damping=damping).fit(x)
	predicted = af.predict(x)
	if len(af.cluster_centers_indices_) > 2:
		sc_score = '%0.2f'%metrics.silhouette_score(x, af.labels_, metric='sqeuclidean')
	else:
		sc_score = 0
	return (predicted, sc_score)



def cluster_by_haplotype (df, dampings=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9], iterations=100):
	t_df = df.set_index(['Chromosome', 'start', 'end']).T
	X = StandardScaler(with_mean=True).fit_transform(t_df)
	for dmp in dampings:
		af = AffinityPropagation(random_state=0, damping=dmp).fit(X)
		X_predicted = af.fit_predict(X)
		if len(af.cluster_centers_indices_) > 2:
			sc_score = '%0.2f'%metrics.silhouette_score(X, af.labels_, metric='sqeuclidean')
		else:
			sc_score = 0
		t_df.insert(loc=0,  column=f'sc_{dmp}', value=sc_score, allow_duplicates=True)
		t_df.insert(loc=0,  column=f'dmp_{dmp}', value=X_predicted, allow_duplicates=True)
	return t_df