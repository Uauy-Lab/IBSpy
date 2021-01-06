import os
import sys


import logging
import operator
import unittest
from IBSpy import FastaChunkReader, KmerGWASDBBuilder


class TestFasta(unittest.TestCase):

	logger = logging.getLogger("test")
	logger.setLevel(logging.DEBUG)
	
	def setUp(self):
		self.data_path="./tests/data/"

	def test_fasta_iterator(self):
		fasta_iter = FastaChunkReader(self.data_path + "/short_test.fa", chunk_size = 50, kmer_size=31)
		#sizes = []
		seqs=[
			'ATGGGCAAGGGCAGGAATAAGAAGAAGAAGGCGCCCGCCGCCGCCGCCCC',
			'AGAAGAAGAAGGCGCCCGCCGCCGCCGCCCCCTCCCATGTAGCGCTCAAC',
			'CGCCGCCGCCCCCTCCCATGTAGCGCTCAACGGCGACGAGCCGCAGCAGG',
			'GTAGCGCTCAACGGCGACGAGCCGCAGCAGGTGGATTCGGCAGCGGCCAG',
			'AGCCGCAGCAGGTGGATTCGGCAGCGGCCAGCGGAGACGCGATT',
			'ATGGCGCCGACTGCCCAGGTACCACCTCCCTCACCAACCCCCGCTGNNNN',
			'TACCACCTCCCTCACCAACCCCCGCTGNNNNNNNNNNNNNNCAGTTATTA',
			'CCCCGCTGNNNNNNNNNNNNNNCAGTTATTAATGTGCATGCCTCATTTCC',
			'NNNCAGTTATTAATGTGCATGCCTCATTTCCTGTTTCTTGGTTCTTCAGA',
			'TGCCTCATTTCCTGTTTCTTGGTTCTTCAGATTGCGGCGCATGA']

		starts=[0,19,38, 57, 76, 0, 19, 38, 57, 76]
		ends=[50,69,88,107,120,50,69,88,107,120]
		seqnames=[
			'chr4B__ari:150492206-150497405',
			'chr4B__ari:150492206-150497405',
			'chr4B__ari:150492206-150497405',
			'chr4B__ari:150492206-150497405',
			'chr4B__ari:150492206-150497405',
			'chr4B__jag:153229422-153234621',
			'chr4B__jag:153229422-153234621',
			'chr4B__jag:153229422-153234621',
			'chr4B__jag:153229422-153234621',
			'chr4B__jag:153229422-153234621']

		i = 0
		for seq in fasta_iter:
			self.assertEqual(seqs[i], seq['seq'])
			self.assertEqual(starts[i], seq['start'])
			self.assertEqual(ends[i], seq['end'])
			self.assertEqual(seqnames[i], seq['seqname'])
			i += 1
			#print("'" + str(seq['seqname']) + "',")

	def test_fasta_to_kmers(self):

		builder = KmerGWASDBBuilder(31)
		fasta_iter = FastaChunkReader(self.data_path + "/short_test.fa", chunk_size = 50, kmer_size=31)
		sizes = [20, 20, 20, 20, 14, 16, 0, 0, 17, 14]
		i = 0
		for seq in fasta_iter: 
			kmers = builder.sequence_to_kmers(seq['seq'])
			self.assertEqual(sizes[i], len(kmers))
			i += 1
			
			


