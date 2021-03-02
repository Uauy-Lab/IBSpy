import os
import sys

import logging
import operator
import unittest
from IBSpy import IBSpyResults 

class TestResults(unittest.TestCase):
	logger = logging.getLogger("test")
	logger.setLevel(logging.DEBUG)

	def setUp(self):
		db = "./tests/data/test_kmeribs-Wheat_Jagger-Flame.tsv" 
		lens= "./tests/data/test_chr_sizes_jagger.tsv"
		windows = 200000
		self.results = IBSpyResults(db, lens, windows)

	def test_windows(self):
		pd = self.results.count_by_windows()

		self.assertEqual(pd.iloc[0]['variations'], 243)
		self.assertEqual(pd.iloc[1]['variations'], 180)
		self.assertEqual(pd.iloc[2]['variations'],  48)
		self.assertEqual(pd.iloc[3]['variations'], 103)
		self.assertEqual(pd.iloc[4]['variations'],  66)
		self.assertEqual(pd.iloc[5]['variations'],  41)
		self.assertEqual(pd.iloc[6]['variations'],  64)

		self.assertEqual(pd.iloc[0]['seqname'], "chr1A_WhJag")
		self.assertEqual(pd.iloc[1]['seqname'], "chr1A_WhJag")
		self.assertEqual(pd.iloc[2]['seqname'], "chr2A_WhJag")
		self.assertEqual(pd.iloc[3]['seqname'], "chr2A_WhJag")
		self.assertEqual(pd.iloc[4]['seqname'], "chr2A_WhJag")
		self.assertEqual(pd.iloc[5]['seqname'], "chr1D_WhJag")
		self.assertEqual(pd.iloc[6]['seqname'], "chr1D_WhJag")

		self.assertEqual(pd.iloc[0]['start'], 200000)
		self.assertEqual(pd.iloc[1]['start'], 400000)
		self.assertEqual(pd.iloc[2]['start'], 400000)
		self.assertEqual(pd.iloc[3]['start'], 600000)
		self.assertEqual(pd.iloc[4]['start'], 800000)
		self.assertEqual(pd.iloc[5]['start'], 800000)å
		self.assertEqual(pd.iloc[6]['start'], 1000000)

		#print(pd.iloc[0]['variations'])

#  seqname    start  variations
# 0  chr1A_WhJag   200000         192
# 1  chr1A_WhJag   400000         180
# 0  chr2A_WhJag   400000          48
# 1  chr2A_WhJag   600000         103
# 2  chr2A_WhJag   800000          66
# 0  chr1D_WhJag   800000          41
# 1  chr1D_WhJag  1000000          64

