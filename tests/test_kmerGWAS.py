import os
import sys
#import pathlib


add_path = os.path.dirname(os.path.abspath(__file__))+ "/.." 

print(add_path)
sys.path.insert(0, add_path)

import logging
import operator
import unittest
import IBSpy 

class TestKmerGWAS(unittest.TestCase):

	logger = logging.getLogger("test")
	logger.setLevel(logging.DEBUG)

	def setUp(self):
		self.assertEqual(1, 1)

	def test_build_kmer(self):
		self.assertEqual(1,1)

if __name__ == '__main__':
	unittest.main()