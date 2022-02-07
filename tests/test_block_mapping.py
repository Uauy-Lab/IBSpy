from operator import delitem
import unittest

import pandas as pd
import pandas.testing as pdt
import IBSpy


class TestBlockMapping(unittest.TestCase):
	def setUp(self) -> None:
		self.mapping_file = "./tests/data/affinity/jag_chi_test_windows.tsv"
		self.mapping_file_expected_1 = "./tests/data/affinity/jag_chi_test_windows_expected1.tsv"
		self.mapping_file_expected_2 = "./tests/data/affinity/jag_chi_test_windows_expected2.tsv"
		self.mapped_result_1 = "./tests/data/affinity/out/jag_chi_test_windows_1.tsv"
		self.mapped_result_2 = "./tests/data/affinity/out/jag_chi_test_windows_2.tsv"
		self.merged_file_expected = "./tests/data/affinity/jag_chi_test_merge.tsv"
		self.merged_reult = "./tests/data/affinity/out/jag_chi_test_merge.tsv"
		self.block_mapping = IBSpy.BlockMapping(self.mapping_file)

	def compare_dfs(self, path_1, path_2, extras=True):
		mapped =   pd.read_csv(path_1, delimiter="\t")
		expected = pd.read_csv(path_2, delimiter="\t")
		pdt.assert_series_equal(mapped["Chromosome"], expected["Chromosome"], check_names=False)
		pdt.assert_series_equal(mapped["Start"], expected["Start"], check_names=False)
		pdt.assert_series_equal(mapped["End"], expected["End"], check_names=False)
		if extras:
			pdt.assert_series_equal(mapped["reference"], expected["reference"], check_names=False)
			pdt.assert_series_equal(mapped["assembly"], expected["assembly"], check_names=False)
			pdt.assert_series_equal(mapped["block_no"], expected["block_no"], check_names=False)

	def test_mapping_for_range(self):
		mapped = self.block_mapping.mapping_for_range("chr1A__chi", 100001, 300000)
		mapped.to_csv(path=self.mapped_result_1, sep="\t")
		self.compare_dfs(self.mapped_result_1,self.mapping_file_expected_1  )
		
		mapped = self.block_mapping.mapping_for_range("chr1A__chi", 100001, 300000, assembly="chinese")
		mapped.to_csv(path=self.mapped_result_2, sep="\t") 
		self.compare_dfs(self.mapped_result_2,self.mapping_file_expected_2  )

	def test_all_regions_for(self):
		merged = self.block_mapping.all_regions_for("chr1A__chi", 100001, 300000, assembly="chinese")
		merged.to_csv(path=self.merged_reult, sep="\t")
		self.compare_dfs(self.merged_reult,self.merged_file_expected, extras=False  )

if __name__ == '__main__':
    unittest.main()