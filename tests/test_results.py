import os
import sys
import numpy as np


import logging
import operator
import unittest
from IBSpy import IBSpyResults 

class TestResults(unittest.TestCase):
    logger = logging.getLogger("test")
    logger.setLevel(logging.DEBUG)

    def setUp(self):

        db = "./tests/plots_data/test_kmeribs-Wheat_Jagger-Flame.tsv" 
        windows = 200000
        filter_counts = 500
        self.stitch_number = 3
        self.n_components = 3
        self.covariance_type = 'full'
        self.results = IBSpyResults(db, windows, filter_counts)
        # print(self.results.count_by_windows())
    def test_windows(self):
        pd = self.results.count_by_windows()

        self.assertEqual(pd.iloc[0]['variations'], 243)
        self.assertEqual(pd.iloc[1]['variations'], 180)
        self.assertEqual(pd.iloc[2]['variations'],  41)
        self.assertEqual(pd.iloc[3]['variations'], 64)
        self.assertEqual(pd.iloc[4]['variations'],  48)
        self.assertEqual(pd.iloc[5]['variations'],  103)
        self.assertEqual(pd.iloc[6]['variations'],  66)

        self.assertEqual(pd.iloc[0]['seqname'], "chr1A_WhJag")
        self.assertEqual(pd.iloc[1]['seqname'], "chr1A_WhJag")
        self.assertEqual(pd.iloc[2]['seqname'], "chr1D_WhJag")
        self.assertEqual(pd.iloc[3]['seqname'], "chr1D_WhJag")
        self.assertEqual(pd.iloc[4]['seqname'], "chr2A_WhJag")
        self.assertEqual(pd.iloc[5]['seqname'], "chr2A_WhJag")
        self.assertEqual(pd.iloc[6]['seqname'], "chr2A_WhJag")

        self.assertEqual(pd.iloc[0]['window'], 200000)
        self.assertEqual(pd.iloc[1]['window'], 400000)
        self.assertEqual(pd.iloc[2]['window'], 800000)
        self.assertEqual(pd.iloc[3]['window'], 1000000)
        self.assertEqual(pd.iloc[4]['window'], 400000)
        self.assertEqual(pd.iloc[5]['window'], 600000)
        self.assertEqual(pd.iloc[6]['window'], 800000)


    def test_transform_counts_to_log(self):
        counts = self.results.count_by_windows()
        log_test, pd = self.results.transform_counts_to_log(counts)
#         print(log_test.tolist())
        self.assertEqual(log_test[0], 5.493061443340548)
        self.assertEqual(log_test[6], 4.189654742026425)

    def test_build_model(self):
        counts = self.results.count_by_windows()
        log_test, pd = self.results.transform_counts_to_log(counts)
        model = self.results.build_gmm_model(log_test, pd,self.n_components, self.covariance_type) 
        print("....\n")
        print(model)
        # self.assertEqual(model.iloc[0]['v_gmm'],  0)
        # self.assertEqual(model.iloc[1]['v_gmm'],  0)
        # self.assertEqual(model.iloc[2]['v_gmm'],  1)
        # self.assertEqual(model.iloc[3]['v_gmm'],  1)
        # self.assertEqual(model.iloc[4]['v_gmm'],  1)
        # self.assertEqual(model.iloc[5]['v_gmm'],  0)
        # self.assertEqual(model.iloc[6]['v_gmm'],  0)
        #TODO: Finish ths test

    def test_stitch_gmm_haplotypes(self):
        counts = self.results.count_by_windows()
        log_test, pd = self.results.transform_counts_to_log(counts)
        model = self.results.build_gmm_model(log_test, pd,self.n_components, self.covariance_type) 
        hap_pd = self.results.stitch_gmm_haplotypes(model, self.stitch_number)
        # print(hap_pd)
        
        self.assertEqual(hap_pd.iloc[0]['mean'], 48.60)
        self.assertEqual(hap_pd.iloc[1]['mean'], 60.00)
        self.assertEqual(hap_pd.iloc[2]['mean'], 20.50)
        self.assertEqual(hap_pd.iloc[3]['mean'], 32.00)
        self.assertEqual(hap_pd.iloc[4]['mean'], 48.00)
        self.assertEqual(hap_pd.iloc[5]['mean'], 25.75)
        self.assertEqual(hap_pd.iloc[6]['mean'], 33.00)


if __name__ == '__main__':
    unittest.main()