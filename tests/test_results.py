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
        db = "./tests/data/test_kmeribs-Wheat_Jagger-Flame.tsv" 
        #lens= "./tests/data/test_chr_sizes_jagger.tsv"
        windows = 200000
        v_filter = 500
        self.stitch_number = 3
        self.n_components = 3
        self.covariance_type = 'full'
        self.results = IBSpyResults(db, windows, v_filter)

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

    def test_normalize_data(self):
        log_test, pd = self.results.normalize_data()
#         print(log_test.tolist())
        self.assertEqual(log_test[0], np.array([5.493061443340548]))
        self.assertEqual(log_test[6], np.array([4.189654742026425]))

    def test_stitch_haploBlocks(self):
        hap_pd = self.results.stitch_haploBlocks(self.n_components,self.covariance_type, self.stitch_number)
        print(hap_pd)
        
        self.assertEqual(hap_pd.iloc[0]['v_mean'], 48)
        self.assertEqual(hap_pd.iloc[1]['v_mean'], 60)
        self.assertEqual(hap_pd.iloc[2]['v_mean'], 20)
        self.assertEqual(hap_pd.iloc[3]['v_mean'], 32)
        self.assertEqual(hap_pd.iloc[4]['v_mean'], 48)
        self.assertEqual(hap_pd.iloc[5]['v_mean'], 25)
        self.assertEqual(hap_pd.iloc[6]['v_mean'], 33)


		#print(pd.iloc[0]['variations'])

#  seqname    start  variations
# 0  chr1A_WhJag   200000         192
# 1  chr1A_WhJag   400000         180
# 0  chr2A_WhJag   400000          48
# 1  chr2A_WhJag   600000         103
# 2  chr2A_WhJag   800000          66
# 0  chr1D_WhJag   800000          41
# 1  chr1D_WhJag  1000000          64

