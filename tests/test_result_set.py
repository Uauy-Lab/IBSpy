from fileinput import filename
import logging
import unittest

from IBSpy import IBSpyResultsSet


class TestResultSet(unittest.TestCase):
    logger = logging.getLogger("test")
    logger.setLevel(logging.DEBUG)

    def setUp(self):
        self.samples_metadata="./tests/data/affinity/samples_metadata.tsv"

    def test_load_data(self):
        ibspy_results = IBSpyResultsSet(filename=self.samples_metadata)
        self.assertEqual(len(ibspy_results.samples_df.index), 8)

        ibspy_results = IBSpyResultsSet(filename=self.samples_metadata, references=["chinese"], samples=["Agron", "ability"])
        self.assertEqual(len(ibspy_results.samples_df.index), 2)
    
    def test_values_matrix(self):
        ibspy_results = IBSpyResultsSet(filename=self.samples_metadata)
        ibspy_results.values_matrix
        matrix = ibspy_results.values_matrix
        print(matrix)

if __name__ == '__main__':
    unittest.main()