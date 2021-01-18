import os
import sys
#import pathlib
#add_path = os.path.dirname(os.path.abspath(__file__))+ "/.." 
#print(add_path)
#sys.path.insert(0, add_path)

import logging
import operator
import unittest
from IBSpy import JellyfishSDB, FastaChunkReader 

loaded=False
try:
    import dna_jellyfish as jf
    loaded=True
except ImportError:
    pass

class TestJellyfish(unittest.TestCase):

	logger = logging.getLogger("test")
	logger.setLevel(logging.DEBUG)

	def setUp(self):
		self.data_path="./tests/data/"

			
	# def test_compare_kmers(self):
	# 	kmer_builder = KmerGWASDBBuilder(31)
	# 	tests = [
	# 	'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
	# 	'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
	# 	'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',
	# 	'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
	# 	'AAAAAAAAAGGGAAAGGGAAGAGGGAAAGGG']
	# 	self.assertEqual(kmer_builder.compare(tests[0], tests[1]),  0)
	# 	self.assertEqual(kmer_builder.compare(tests[2], tests[1]),  1)
	# 	self.assertEqual(kmer_builder.compare(tests[1], tests[2]), -1)

	def run_db_tests(self, kmerdb):
		references=['jagger', 'arinalrfor', 'julius', 'lancer', 'landmark', 'mace', 'norin61', 'spelta', 'stanley', 'sy_mattis']
		expected = [
		[2970, 2201], [2939, 2170], [2970, 2201], [2609],
		[2970, 2201], [2970, 2201], [2970, 2201], [2970, 2201],
		[2879, 2176], [2970, 2201]
		]
		i = 0
		for r in references:	
			path = self.data_path + "/test4B." + r + ".fa"
			windows = kmerdb.kmers_in_windows(path, window_size=3000)
			obs =  map(lambda w: w.observed_kmers, windows)
			j = 0
			for o in obs:
				self.assertEqual(expected[i][j], o)
				j += 1
			i += 1



	def test_jellyfish_db(self):
		if not loaded:
			print("JellyFish python module not installed", file=sys.stderr)
		else:
			kmerdb = JellyfishSDB(31)
			kmerdb.load(self.data_path + "/test4B.jagger.fa.k31.jf")
			self.run_db_tests(kmerdb)


if __name__ == '__main__':
	unittest.main()