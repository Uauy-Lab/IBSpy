import os
import sys
#import pathlib
#add_path = os.path.dirname(os.path.abspath(__file__))+ "/.." 
#print(add_path)
#sys.path.insert(0, add_path)

import logging
import operator
import unittest
from IBSpy import KmerGWASDBBuilder, KmerGWASDB, FastaChunkReader 

class TestKmerGWAS(unittest.TestCase):

	logger = logging.getLogger("test")
	logger.setLevel(logging.DEBUG)

	def setUp(self):
		self.data_path="./tests/data/"

	def test_build_kmer(self):
		test = [
		'AAAAAAAAA',
		'TAAAAAAAA',
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
		'TTTTTTTTT',
		]
		kmer_builder = KmerGWASDBBuilder(9)
		
		for x in test:
			tested = kmer_builder.string_to_kmer(x)
			back = kmer_builder.kmer_to_string(tested)
			self.assertEqual(x, back)
			tested = kmer_builder.string_to_kmer(x.lower())
			back = kmer_builder.kmer_to_string(tested)
			self.assertEqual(x, back)
			#print (x +  "->"  + str(tested)  + "->"  + back) 
		test = [
		'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
		'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
		'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',
		'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
		'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
		'TTTTTTTTTTTTCCCCCCCTTTTTTTTTTTT',
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

	def test_build_kmer_db(self):
		kmerdb = KmerGWASDB(31)
		kmer_builder = KmerGWASDBBuilder(31)
		kmerdb.load_from_fasta(self.data_path + "/short_test_duplicate.fa")
		#Appears twice: GGGCAGGAATAAGAAGAAGAAGGCGCCCGCC GTTTCTTGGTTCTTCAGATTGCGGCGCATGA
		#print("Size: " + str(len(kmerdb)))
		self.assertEqual(len(kmerdb), 120)

		#for i in range(0,len(kmerdb)):
			#print(kmerdb[i])
			#print(kmer_builder.kmer_to_string(kmerdb[i]))
		#	pass
		
		kmer_builder2 = KmerGWASDBBuilder(31)
		
		test = [
		"GGGCAGGAATAAGAAGAAGAAGGCGCCCGCC",
		"GTTTCTTGGTTCTTCAGATTGCGGCGCATGA",
		'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
		'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT']
		results = [True, True, False, False]
		i = 0
		for x in test:
			tested = kmer_builder2.string_to_kmer(x)
			#print(x)
			self.assertEqual (tested in kmerdb ,results[i])
			i+=1
			
	def test_compare_kmers(self):
		kmer_builder = KmerGWASDBBuilder(31)
		tests = [
		'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
		'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT',
		'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG',
		'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
		'AAAAAAAAAGGGAAAGGGAAGAGGGAAAGGG']
		self.assertEqual(kmer_builder.compare(tests[0], tests[1]),  0)
		self.assertEqual(kmer_builder.compare(tests[2], tests[1]),  1)
		self.assertEqual(kmer_builder.compare(tests[1], tests[2]), -1)

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
			obs =  map(lambda w: w['observed_kmers'], windows)
			j = 0
			for o in obs:
				self.assertEqual(expected[i][j], o)
				j += 1
			i += 1



	def test_kmer_db(self):
		kmerdb = KmerGWASDB(31)
		kmerdb.load_from_fasta(self.data_path + "/test4B.jagger.fa")
		self.run_db_tests(kmerdb)
		
		saved_file = self.data_path + "test4B.jagger.temp.kmerGWAS_k31"
		kmerdb.save(saved_file)
		kmerdb_read = KmerGWASDB(31)
		kmerdb_read.load(saved_file)
		self.run_db_tests(kmerdb_read)

		saved_file = self.data_path + "test4B.jagger.kmerGWAS_k31"
		kmerdb.save(saved_file)
		kmerdb_read = KmerGWASDB(31)
		kmerdb_read.load(saved_file)
		self.run_db_tests(kmerdb_read)

if __name__ == '__main__':
	unittest.main()