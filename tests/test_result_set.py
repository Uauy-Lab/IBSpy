from asyncio.log import logger
from fileinput import filename
import logging
import os
import unittest
import pandas as pd
import pandas.testing as pdt
import shutil
import IBSpy



class TestResultSet(unittest.TestCase):
    logger = logging.getLogger("test")
    logger.setLevel(logging.DEBUG)

    def setUp(self):
        self.out_folder="./tests/data/affinity/out/"
        try:
            os.mkdir(self.out_folder)
        except OSError as error:
            pass
        self.samples_metadata="./tests/data/affinity/samples_metadata.tsv"
        self.combine_resut="./tests/data/affinity/chr1A_variations_5000w.tsv"
        self.combine_resut_renamed="./tests/data/affinity/chr1A_variations_5000w_renamed.tsv"
        self.combine_resut_out="./tests/data/affinity/out/chr1A_variations_5000w.tsv"
        self.combine_resut_renamed_out="./tests/data/affinity/out/chr1A_variations_5000w_renamed.tsv"
        self.mapping_path = "./tests/data/affinity/chi_jag_mapping.tsv"
        self.options = IBSpy.IBSpyOptions()
        self.options.metadata = self.samples_metadata
        self.block_mapping_file = "./tests/data/affinity/jag_chi_test_windows.tsv"
        self.mapped_window_result1 = "./tests/data/affinity/mapped_window_1.tsv"
        self.mapped_window_result2 = "./tests/data/affinity/mapped_window_2.tsv"
        self.mapped_window_result1_out = "./tests/data/affinity/out/mapped_window_1.tsv"
        self.mapped_window_result2_out = "./tests/data/affinity/out/mapped_window_2.tsv"

        try:
            shutil.rmtree("./tests/data/affinity/out")
        except OSError as error:
            pass
        os.mkdir("./tests/data/affinity/out")

    def compare_dfs(self, path_1, path_2, extras=[]):
        mapped =   pd.read_csv(path_1, delimiter="\t")
        expected = pd.read_csv(path_2, delimiter="\t")
        pdt.assert_series_equal(mapped["Chromosome"], expected["Chromosome"], check_names=False)
        pdt.assert_series_equal(mapped["Start"], expected["Start"], check_names=False)
        pdt.assert_series_equal(mapped["End"], expected["End"], check_names=False)
        for e in extras:
            pdt.assert_series_equal(mapped[e], expected[e], check_names=False)


    def test_load_data(self):
        
        ibspy_results = IBSpy.IBSpyResultsSet(options= self.options)
        self.assertEqual(len(ibspy_results.samples_df.index), 8)
        self.options.references = ["chinese"]
        self.options.samples=["Agron", "ability"]
        ibspy_results = IBSpy.IBSpyResultsSet(options=self.options )
        self.assertEqual(len(ibspy_results.samples_df.index), 2)
    
    def test_values_matrix(self):
        expected = pd.read_csv(self.combine_resut, delimiter='\t')
        ibspy_results = IBSpy.IBSpyResultsSet(options=self.options)
        matrix = ibspy_results.values_matrix
        matrix.to_csv(self.combine_resut_out, sep='\t')
        matrix = pd.read_csv(self.combine_resut_out, delimiter='\t')
        pdt.assert_series_equal(matrix["Chromosome"], expected["Chromosome"], check_names=False)
    
    def test_values_matrix_rename(self):
        expected = pd.read_csv(self.combine_resut_renamed, delimiter='\t')
        self.options.chromosome_mapping = self.mapping_path
        ibspy_results = IBSpy.IBSpyResultsSet(options=self.options)
        matrix = ibspy_results.values_matrix
        matrix.to_csv(self.combine_resut_renamed_out, sep='\t')
        matrix = pd.read_csv(self.combine_resut_renamed_out, delimiter='\t')
        pdt.assert_series_equal(matrix["Chromosome"], expected["Chromosome"], check_names=False)

    def test_values_matrix_chromosome_iterator(self):
        ibspy_results = IBSpy.IBSpyResultsSet(options=self.options)
        expected_seqname = ['chr1A', 'chr1A_WhJag']
        expected_len = [109,109]
        for i, (chr_map, chr )  in enumerate(ibspy_results.values_matrix_seqname_iterator()):
            self.assertEqual(chr, expected_seqname[i])
            self.assertAlmostEqual(len(chr_map), expected_len[i])

    def test_values_matrix_iterator(self):
        ibspy_results = IBSpy.IBSpyResultsSet(options=self.options)
        expected_len = [20, 20, 20, 20, 20, 9, 20, 20, 20, 20, 20, 9]
        for i, block in enumerate(ibspy_results.values_matrix_iterator()):
            self.assertAlmostEqual(len(block), expected_len[i])

    def test_mapped_window(self):
        self.options.chromosome_mapping = self.mapping_path
        ibspy_results = IBSpy.IBSpyResultsSet(options=self.options)
        mapped_window = ibspy_results.mapped_window("chr1A__chi",  100001, 300000, assembly="chinese")
        mapped_window.to_csv(self.mapped_window_result1_out, sep="\t")
        self.compare_dfs(self.mapped_window_result1, self.mapped_window_result1_out, extras=["Agron", "ability", "WATDE0010", "WATDE0020"])

        self.options.block_mapping = self.block_mapping_file
        ibspy_results = IBSpy.IBSpyResultsSet(options=self.options)
        mapped_window = ibspy_results.mapped_window("chr1A__chi",  100001, 300000, assembly="chinese")
        mapped_window.to_csv(self.mapped_window_result2_out, sep="\t")
        self.compare_dfs(self.mapped_window_result2, self.mapped_window_result2_out, extras=["Agron", "ability", "WATDE0010", "WATDE0020"])

    def test_window_iterator(self):
        self.options.chromosome_mapping = self.mapping_path
        self.options.block_mapping = self.block_mapping_file
        ibspy_results = IBSpy.IBSpyResultsSet(options=self.options)

        


if __name__ == '__main__':
    unittest.main()