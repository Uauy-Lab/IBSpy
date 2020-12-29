from .kmer_database import KmerDB, KmerBuilder
#from ctypes import *
import os
import inspect
import kmerGWAS
import cython

class KmerGWASDB(KmerDB):
	pass

class KmerGWASDBBuilder(KmerBuilder):
	#cdef char* _builder
	def __init__(self, kmer_size):
		self.kmer_size = kmer_size
		self._builder = kmerGWAS.KmerGWAS_builder(kmer_size)
		#print(self._builder);

	def string_to_kmer(self, sequence):
		binary_kmer = self._builder.string_to_kmer(sequence)
		return binary_kmer


	def kmer_to_string(self, binary_kmer):
		return self._builder.kmer_to_string(binary_kmer).decode('UTF-8')
