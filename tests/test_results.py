from locale import normalize
import os
import sys
import numpy as np


import logging
import operator
import unittest

from pandas import options
from IBSpy import IBSpyResults
from IBSpy.IBSpy_options import IBSpyOptions 

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
        opts1 = IBSpyOptions()
        opts1.normalize = False
        opts1.score = 'variations'
        opts1.window_size = windows
        self.results = IBSpyResults(db, opts1)
        opts = IBSpyOptions()
        opts.score = "observed_kmers"
        opts.normalize = True
        opts.filter_counts = filter_counts
        opts.window_size = windows
        self.results_norm = IBSpyResults(db, opts)
        # print(self.results.count_by_windows())

    def test_windows(self):
        pd = self.results.count_by_windows()
        # print(pd)
        self.assertEqual(pd.iloc[0]['variations'], 188)
        self.assertEqual(pd.iloc[1]['variations'], 235)
        self.assertEqual(pd.iloc[2]['variations'],  22)
        self.assertEqual(pd.iloc[3]['variations'],  83)
        self.assertEqual(pd.iloc[4]['variations'],  130)
        self.assertEqual(pd.iloc[5]['variations'],  87)


        self.assertEqual(pd.iloc[0]['seqname'], "chr1A_WhJag")
        self.assertEqual(pd.iloc[1]['seqname'], "chr1A_WhJag")
        self.assertEqual(pd.iloc[2]['seqname'], "chr1D_WhJag")
        self.assertEqual(pd.iloc[3]['seqname'], "chr1D_WhJag")
        self.assertEqual(pd.iloc[4]['seqname'], "chr2A_WhJag")
        self.assertEqual(pd.iloc[5]['seqname'], "chr2A_WhJag")

        self.assertEqual(pd.iloc[0]['end'], 200000)
        self.assertEqual(pd.iloc[1]['end'], 400000)
        self.assertEqual(pd.iloc[2]['end'], 800000)
        self.assertEqual(pd.iloc[3]['end'], 1000000)
        self.assertEqual(pd.iloc[4]['end'], 600000)
        self.assertEqual(pd.iloc[5]['end'], 800000)
    

    def test_windows_normalised(self):
        # print(self.results_norm.db)
        # print(self.results_norm.count_by_windows())
        pd = self.results_norm.count_by_windows() 
        pd['mean'] = pd['mean'].round(decimals = 6) 
        self.assertEqual(pd.iloc[0]['mean'], 0.982182)
        self.assertEqual(pd.iloc[1]['mean'], 0.976105)
        self.assertEqual(pd.iloc[2]['mean'], 0.993488)
        self.assertEqual(pd.iloc[3]['mean'], 0.992799)
        self.assertEqual(pd.iloc[4]['mean'], 0.989477)
        self.assertEqual(pd.iloc[5]['mean'], 0.991207)

        self.assertEqual(pd.iloc[0]['seqname'], "chr1A_WhJag")
        self.assertEqual(pd.iloc[1]['seqname'], "chr1A_WhJag")
        self.assertEqual(pd.iloc[2]['seqname'], "chr1D_WhJag")
        self.assertEqual(pd.iloc[3]['seqname'], "chr1D_WhJag")
        self.assertEqual(pd.iloc[4]['seqname'], "chr2A_WhJag")
        self.assertEqual(pd.iloc[5]['seqname'], "chr2A_WhJag")

        self.assertEqual(pd.iloc[0]['end'], 200000)
        self.assertEqual(pd.iloc[1]['end'], 400000)
        self.assertEqual(pd.iloc[2]['end'], 800000)
        self.assertEqual(pd.iloc[3]['end'], 1000000)
        self.assertEqual(pd.iloc[4]['end'], 600000)
        self.assertEqual(pd.iloc[5]['end'], 800000)
    


    def test_transform_counts_to_log(self):
        counts = self.results.count_by_windows()
        log_test, pd = self.results.transform_counts_to_log(counts)
        # print(log_test.tolist())
        self.assertEqual( round( log_test[0][0], 2),  round(5.236441962829949, 2 ) )
        self.assertEqual( round( log_test[1][0], 2),  round(5.459585514144159, 2))
        self.assertEqual( round( log_test[2][0], 2),  round(3.091042453358316, 2))
        self.assertEqual( round( log_test[3][0], 2),  round(4.418840607796598, 2))
        self.assertEqual( round( log_test[4][0], 2),  round(4.867534450455582, 2))
        self.assertEqual( round( log_test[5][0], 2),  round(4.465908118654584, 2 ) )

    # def test_build_model(self):
    #     counts = self.results.count_by_windows()
    #     log_test, pd = self.results.transform_counts_to_log(counts)
    #     model = self.results.build_gmm_model(log_test, pd,self.n_components, self.covariance_type) 
    #     # print("....\n")
    #     # print(model)
    #     # self.assertEqual(model.iloc[0]['v_gmm'],  0)
    #     # self.assertEqual(model.iloc[1]['v_gmm'],  0)
    #     # self.assertEqual(model.iloc[2]['v_gmm'],  1)
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
        hap_pd['mean'] = hap_pd['mean'].round(decimals = 6) 
        # print("\n")
        # print(hap_pd)
        
        self.assertEqual(hap_pd.iloc[0]['mean'], 47.000000)
        self.assertEqual(hap_pd.iloc[1]['mean'], 58.750000)
        self.assertEqual(hap_pd.iloc[2]['mean'], 22.000000)
        self.assertEqual(hap_pd.iloc[3]['mean'], 27.666667)
        self.assertEqual(hap_pd.iloc[4]['mean'], 32.500000)
        self.assertEqual(hap_pd.iloc[5]['mean'], 29.000000)
      


if __name__ == '__main__':
    unittest.main()