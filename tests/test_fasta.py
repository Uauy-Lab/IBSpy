import os
import sys


import logging
import operator
import unittest
from IBSpy import FastaChunkReader


class TestFasta(unittest.TestCase):

	logger = logging.getLogger("test")
	logger.setLevel(logging.DEBUG)
	
	def setUp(self):
		self.data_path="./tests/data/"

	def test_fasta_iterator(self):
		fasta_iter = FastaChunkReader(self.data_path + "/short_test.fa", chunk_size = 50, kmer_size=31)
		for seq in fasta_iter:
			#print(seq)
			print(seq['seq'])
			print(len(seq['seq']))