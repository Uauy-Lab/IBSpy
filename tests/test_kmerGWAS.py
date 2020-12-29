import os
import sys
#import pathlib
#add_path = os.path.dirname(os.path.abspath(__file__))+ "/.." 
#print(add_path)
#sys.path.insert(0, add_path)

import logging
import operator
import unittest
from IBSpy import KmerGWASDBBuilder 

class TestKmerGWAS(unittest.TestCase):

	logger = logging.getLogger("test")
	logger.setLevel(logging.DEBUG)

	def setUp(self):
		self.assertEqual(1, 1)
		

	def test_build_kmer(self):
		test = [
		'AAAAAAAAA',
		'AAAAAAAAT',
		'AAAAAAATA',
		'AAAAAATAA',
		'AAAAATAAA',
		'AAAATAAAA',
		'AAATAAAAA',
		'AATAAAAAA',
		'ATAAAAAAA',
		'TAAAAAAAA',
		'TTTTTTTTT',
		'GGGGGGGGG',
		'GGGGGGGGT',
		'GGGGGGGTG',
		'GGGGGGTGG',
		'GGGGGTGGG',
		'CCCCTCCCC',
		'AAATCCCCC',
		'AATCCCCCC',
		'ATCCCCCCC',
		'TCCCCACCC',
		#'TTTTTTTTT',
		]
		kmer_builder = KmerGWASDBBuilder(9)
		
		for x in test:
			tested = kmer_builder.string_to_kmer(x)
			back = kmer_builder.kmer_to_string(tested)
			self.assertEqual(x, back)
			#print (x +  "->"  + str(tested)  + "->"  + back) 
		test = [
		'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
		'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
		'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',
		'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
		'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
		'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
		'AAAAAAAAAGGGAAAGGGAAGAGGGAAAGGG',
		'AAAAAAAATGGAAAGGGAAGAGGGAAAGGGA',
		'AAAAAAATAGGAAAGGGAAGAGGGAAAGGGA',
		'AAAAAATAAGGAAAGGGAAGAGGGAAAGGGA',
		'AAAAATAAAGGAAAGGGAAGAGTTCTTGGGA',
		'AAAATAAAATTCTTGGGAAGAGTTCTTGGGA',
		'AAATAAAAAGGAAAGGGAAGAGGGAAAGGGA',
		'AATAAAAAAGGAAAGGGAAGAGGGAAAGGGA',
		'ATAAAAAAAGGAAAGGGAAGAGGGAAAGGGA',
		'TAAAAAAAAGGAAAGGGGGGCCTTAAAGGGA',
		'TTTTTTTTTGGAAAGGGGGGCCTTAAAGGGA',
		'GGGGGGGGGGGAAAGGGGGGCCTTAAAGGGA',
		'GGGGGGGGTGGAAAGGGAAGAGGGAAAGGGA',
		'GGGGGGGTGGGAAAGGGAAGAGGGAAAGGGA',
		'GGGGGGTGGGGAAAGGGAAGAGGGAAAGGGA',
		'GGGGGTGGGGGAAAGGGAAGAGGGAAAGGGA',
		'CCCCTCCCCGGAAAGGGAAGAGGGAAAGGGA',
		'AAATCCCCCGGAAAGGGAAGAGGGAAAGGGA',
		'AATCCCCCCGGAAAGGGAAGAGGGAAAGGGA',
		'ATCCCCCCCGGAAAGGGAAGAGGGAAAGGGA',
		'TCCCCACCCGGAAAGGGAAGAGGGAAAGGGA',
		'TTTTTTTTTGGAAAGGGAAGAGGGAAAGGGA'
		]

		kmer_builder2 = KmerGWASDBBuilder(31)
		
		for x in test:
			tested = kmer_builder2.string_to_kmer(x)
			back   = kmer_builder2.kmer_to_string(tested)
			#print(len(x))
			#print (x +  "->"  + str(tested)  + "->"  + back)
			self.assertEqual(x, back)
			 



if __name__ == '__main__':
	unittest.main()