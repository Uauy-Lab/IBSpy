from .kmer_database import KmerDB, KmerBuilder
from ctypes import *
#cdll.LoadLibrary("build/lib.macosx-11-x86_64-3.9/kmerGWAS.cpython-39-darwin.so") 
from .kmerGWAS import *

class KmerGWASDB(KmerDB):
	pass

class KmerGWASDBBuilder(KmerBuilder):
	def string_to_kmer(sequence):
		raise NotImplementedError

	def kmer_to_string(sequence):
		raise NotImplementedError
