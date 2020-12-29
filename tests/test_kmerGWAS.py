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
		self.assertEqual(1,1)
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
		'TTTTTTTTT'
		]

		kmer_builder = KmerGWASDBBuilder(9)
		for x in test:
			tested = kmer_builder.string_to_kmer(x)
			back = kmer_builder.kmer_to_string(tested)
			print (x +  "->"  + str(tested)  + "->"  + back) 



if __name__ == '__main__':
	unittest.main()