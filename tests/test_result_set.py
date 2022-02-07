from asyncio.log import logger
from fileinput import filename
import logging
import os
import unittest
import pandas
import pandas.testing as pdt
import shutil
import IBSpy



class TestResultSet(unittest.TestCase):
    logger = logging.getLogger("test")
    logger.setLevel(logging.DEBUG)

    def setUp(self):
        self.samples_metadata="./tests/data/affinity/samples_metadata.tsv"
        self.combine_resut="./tests/data/affinity/chr1A_variations_5000w.tsv"
        self.combine_resut_renamed="./tests/data/affinity/chr1A_variations_5000w_renamed.tsv"
        self.combine_resut_out="./tests/data/affinity/out/chr1A_variations_5000w.tsv"
        self.combine_resut_renamed_out="./tests/data/affinity/out/chr1A_variations_5000w_renamed.tsv"
        self.mapping_path = "./tests/data/affinity/chi_jag_mapping.tsv"
        self.options = IBSpy.IBSpyOptions()
        self.options.metadata = self.samples_metadata

        try:
            shutil.rmtree("./tests/data/affinity/out")
        except OSError as error:
            pass
        os.mkdir("./tests/data/affinity/out")

    def test_load_data(self):
        
        ibspy_results = IBSpy.IBSpyResultsSet(options= self.options)
        self.assertEqual(len(ibspy_results.samples_df.index), 8)
        self.options.references = ["chinese"]
        self.options.samples=["Agron", "ability"]
        ibspy_results = IBSpy.IBSpyResultsSet(options=self.options )
        self.assertEqual(len(ibspy_results.samples_df.index), 2)
    
    def test_values_matrix(self):
        expected = pandas.read_csv(self.combine_resut, delimiter='\t')
        ibspy_results = IBSpy.IBSpyResultsSet(options=self.options)
        matrix = ibspy_results.values_matrix
        matrix.to_csv(self.combine_resut_out, sep='\t', index=False )
        matrix = pandas.read_csv(self.combine_resut_out, delimiter='\t')
        pdt.assert_series_equal(matrix["Chromosome"], expected["Chromosome"], check_names=False)
    
    def test_values_matrix_rename(self):
        expected = pandas.read_csv(self.combine_resut_renamed, delimiter='\t')
        self.options.chromosome_mapping = self.mapping_path
        ibspy_results = IBSpy.IBSpyResultsSet(options=self.options)
        matrix = ibspy_results.values_matrix
        matrix.to_csv(self.combine_resut_renamed_out, sep='\t', index=False )
        matrix = pandas.read_csv(self.combine_resut_renamed_out, delimiter='\t')
        pdt.assert_series_equal(matrix["Chromosome"], expected["Chromosome"], check_names=False)

    def test_values_matrix_chromosome_iterator(self):
        ibspy_results = IBSpy.IBSpyResultsSet(options=self.options)
        expected_seqname = ['chr1A', 'chr1A_WhJag']
        expected_len = [109,109]
        for i, chr_map in enumerate(ibspy_results.values_matrix_seqname_iterator()):
            chromosomes = chr_map['Chromosome'].unique()
            self.assertEqual(chromosomes[0], expected_seqname[i])
            self.assertAlmostEqual(len(chr_map), expected_len[i])

    def test_values_matrix_iterator(self):
        ibspy_results = IBSpy.IBSpyResultsSet(options=self.options)
        expected_len = [20, 20, 20, 20, 20, 9, 20, 20, 20, 20, 20, 9]
        for i, block in enumerate(ibspy_results.values_matrix_iterator()):
            self.assertAlmostEqual(len(block), expected_len[i])

    

if __name__ == '__main__':
    unittest.main()