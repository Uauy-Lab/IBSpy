from fileinput import filename
import logging
import os
import unittest
import pandas
import pandas.testing as pdt
import shutil
from IBSpy import IBSpyResultsSet


class TestResultSet(unittest.TestCase):
    logger = logging.getLogger("test")
    logger.setLevel(logging.DEBUG)

    def setUp(self):
        self.samples_metadata="./tests/data/affinity/samples_metadata.tsv"
        self.combine_resut="./tests/data/affinity/chr1A_variations_5000w.tsv"
        self.combine_resut_out="./tests/data/affinity/out/chr1A_variations_5000w.tsv"
        shutil.rmtree("./tests/data/affinity/out")
        os.mkdir("./tests/data/affinity/out")

    def test_load_data(self):
        ibspy_results = IBSpyResultsSet(filename=self.samples_metadata)
        self.assertEqual(len(ibspy_results.samples_df.index), 8)

        ibspy_results = IBSpyResultsSet(filename=self.samples_metadata, references=["chinese"], samples=["Agron", "ability"])
        self.assertEqual(len(ibspy_results.samples_df.index), 2)
    
    def test_values_matrix(self):
        expected = pandas.read_csv(self.combine_resut, delimiter='\t')
        ibspy_results = IBSpyResultsSet(filename=self.samples_metadata)
        matrix = ibspy_results.values_matrix
        matrix.to_csv(self.combine_resut_out, sep='\t', index=False )
        matrix = pandas.read_csv(self.combine_resut_out, delimiter='\t')
        pdt.assert_series_equal(matrix["seqname"], expected["seqname"], check_names=False)

if __name__ == '__main__':
    unittest.main()