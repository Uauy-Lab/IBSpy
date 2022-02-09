import unittest
import numpy

import pandas as pd
from pyranges import PyRanges
from sklearn.preprocessing import StandardScaler
from IBSpy.window_affinity_propagation import single_affinity_run


class TestAffinityPropagation(unittest.TestCase):
	def setUp(self) -> None:
		self.test_prefix="./tests/data/affinity/windows"
		pass
	
	def test_single_affinity_run(self):
		df = pd.read_csv(self.test_prefix + "/0.tsv", delimiter="\t")
		gr = PyRanges(df)
		# df = gr.as_df()
		t_df = gr.as_df().set_index(['Chromosome', 'Start', 'End']).T
		x = StandardScaler(with_mean=True).fit_transform(t_df)
		predicted, score = single_affinity_run(x)
		expected_predicted=[0,0,1,1]
		expected_score = 0
		numpy.testing.assert_array_almost_equal(expected_predicted, predicted)
		self.assertEqual(expected_score, score)
		pass

if __name__ == '__main__':
    unittest.main()