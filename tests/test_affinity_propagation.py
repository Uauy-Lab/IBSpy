import math
import os
import unittest
import pandas as pd
from pyranges import PyRanges
from sklearn.preprocessing import StandardScaler
from IBSpy.window_affinity_propagation import single_affinity_run, cluster_by_haplotype
import numpy as np
import pandas.testing as pdt

def clustering_to_df(predicted, varieties):
	df = pd.DataFrame(predicted).T
	df.insert(loc = 0,
          column = 'variety',
          value = varieties)
	return df

class TestAffinityPropagation(unittest.TestCase):
	def setUp(self) -> None:
		self.test_prefix="./tests/data/affinity/windows"
		self.good_region="./tests/data/affinity/good_region_test.tsv"
		self.bad_region="./tests/data/affinity/bad_region_test.tsv" 
		self.good_predicted_expected="./tests/data/affinity/predicted_good.tsv"
		self.bad_predicted_expected="./tests/data/affinity/predicted_bad.tsv" 
		self.out_folder="./tests/data/affinity/out/"
		try:
			os.mkdir(self.out_folder)
		except OSError as error:
			pass

	def test_single_affinity_run(self):
		df = pd.read_csv(self.test_prefix + "/0.tsv", delimiter="\t")
		gr = PyRanges(df)
		t_df = gr.as_df().set_index(['Chromosome', 'Start', 'End']).T
		x = StandardScaler(with_mean=True).fit_transform(t_df)
		predicted, score = single_affinity_run(x)
		expected_predicted=[0,0,1,1]
		# expected_score = math.nan
		np.testing.assert_array_almost_equal(expected_predicted, predicted)
		self.assertTrue(math.isnan(score))
	
	def test_cluster_by_haplotype_short(self):
		df = pd.read_csv(self.test_prefix + "/0.tsv", delimiter="\t")
		gr = PyRanges(df)
		predicted, scores, used_dampings, random_states, varieties = cluster_by_haplotype(gr, seed=42, iterations=10, dampings=[0.6, 0.7])
		expected_predicted=[0,0,1,1]
		expected_predicted2=[0,0,0,0]
		for i, p in enumerate(predicted):
			ep = expected_predicted
			if i == 3:
				ep = expected_predicted2
			np.testing.assert_array_almost_equal(ep, p)
		for s in scores:
			self.assertTrue(math.isnan(s))
		# df2 = pd.DataFrame(predicted)
		# print(df2) 	
	
	def test_cluster_by_haplotype_long(self):
		
		df = pd.read_csv(self.good_region, delimiter="\t")
		gr = PyRanges(df)
		predicted, scores, used_dampings, random_states, varieties = cluster_by_haplotype(gr, seed=42, iterations=3, dampings=[0.5, 0.6, 0.7, 0.8, 0.9])
		df = clustering_to_df(predicted, varieties)
		df.to_csv("./tests/data/affinity/out/predicted_good.tsv", sep="\t", index=False)
		expected=pd.read_csv(self.good_predicted_expected, delimiter="\t")
		df = pd.read_csv("./tests/data/affinity/out/predicted_good.tsv",delimiter="\t")
		pdt.assert_frame_equal(expected, df)
		
		df = pd.read_csv(self.bad_region, delimiter="\t")
		gr = PyRanges(df)
		predicted, scores, used_dampings, random_states, varieties = cluster_by_haplotype(gr, seed=42, iterations=3, dampings=[0.5, 0.6, 0.7, 0.8, 0.9])
		df = clustering_to_df(predicted, varieties)
		df.to_csv("./tests/data/affinity/out/predicted_bad.tsv", sep="\t", index=False)
		expected=pd.read_csv(self.bad_predicted_expected, delimiter="\t")
		df = pd.read_csv("./tests/data/affinity/out/predicted_bad.tsv",delimiter="\t")
		pdt.assert_frame_equal(expected, df)

if __name__ == '__main__':
    unittest.main()