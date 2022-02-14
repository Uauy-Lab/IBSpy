import math
import os
import unittest
import pandas as pd
from pyranges import PyRanges
from sklearn.preprocessing import StandardScaler
from IBSpy.window_affinity_propagation import select_best_cluster, single_affinity_run, cluster_by_haplotype
import numpy as np
import pandas.testing as pdt

def clustering_to_df(runs):
	predicted = map(lambda r: r.predicted, runs )
	varieties = runs[0].varieties
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
		self.unconverged_region="./tests/data/affinity/unconverged_region_test.tsv" 
		self.good_predicted_expected="./tests/data/affinity/predicted_good.tsv"
		self.bad_predicted_expected="./tests/data/affinity/predicted_bad.tsv" 
		self.unconverged_predicted_expected="./tests/data/affinity/unconverged_bad.tsv" 
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
		runs = cluster_by_haplotype(gr, seed=42, iterations=10, dampings=[0.6, 0.7], max_missing=10)
		expected_predicted=[0,0,1,1]
		expected_predicted2=[0,0,0,0]
		for i, p in enumerate(runs):
			ep = expected_predicted
			if i == 3:
				ep = expected_predicted2
			np.testing.assert_array_almost_equal(ep, p.predicted )
			self.assertTrue(math.isnan(p.score))
		
		# df2 = pd.DataFrame(predicted)
		# print(df2) 	

	def run_single_hap_run(self, input, predicted_expected, tmp_out):
		df = pd.read_csv(input, delimiter="\t")
		# print(df)
		gr = PyRanges(df)
		runs = cluster_by_haplotype(gr, seed=42, iterations=50, dampings=[0.5, 0.6, 0.7, 0.8, 0.9], max_missing=5)
		df = clustering_to_df(runs)
		# df.to_csv(predicted_expected, sep="\t", index=False)
		df.to_csv(tmp_out, sep="\t", index=False)
		expected=pd.read_csv(predicted_expected, delimiter="\t")
		df = pd.read_csv(tmp_out,delimiter="\t")
		pdt.assert_frame_equal(expected, df)
		return select_best_cluster(runs)
		
	
	def test_cluster_by_haplotype_long(self):
		best =  self.run_single_hap_run(self.unconverged_region, self.unconverged_predicted_expected,"./tests/data/affinity/out/predicted_unconverged.tsv" )
		print("....")
		print(best.score)
		print(best.damping)
		print(best.random_state)
		print(best.mutual_info_score)
		print(best.number_of_runs)
		print(best.stdev)
		self.assertAlmostEqual(best.score, 0.9891135390808048)
		self.assertAlmostEqual(best.damping, 0.8)
		self.assertAlmostEqual(best.random_state, 2536146026)
		self.assertAlmostEqual(best.mutual_info_score, 1.0)
		self.assertAlmostEqual(best.number_of_runs, 5)
		self.assertAlmostEqual(best.stdev, 0.0)
		
		best =  self.run_single_hap_run(self.good_region, self.good_predicted_expected,"./tests/data/affinity/out/predicted_good.tsv" )
		self.assertAlmostEqual(best.score, 0.9061312647192928)
		self.assertAlmostEqual(best.damping, 0.7)
		self.assertAlmostEqual(best.random_state, 2906402158)
		self.assertAlmostEqual(best.mutual_info_score, 1.0)
		self.assertAlmostEqual(best.number_of_runs, 5)
		self.assertAlmostEqual(best.stdev, 0.0)
		best =  self.run_single_hap_run(self.bad_region, self.bad_predicted_expected,"./tests/data/affinity/out/predicted_bad.tsv" )
		self.assertAlmostEqual(best.score, 0.2153727591312509)
		self.assertAlmostEqual(best.damping, 0.7)
		self.assertAlmostEqual(best.random_state, 2906402158)
		self.assertAlmostEqual(best.mutual_info_score, 1.0)
		self.assertAlmostEqual(best.number_of_runs, 5)
		self.assertAlmostEqual(best.stdev, 1.2161883888976234e-16)
		

if __name__ == '__main__':
    unittest.main()