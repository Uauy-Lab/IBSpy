from .kmer_database import KmerDB, KmerBuilder
from ctypes import *
import os

#libname = os.path.abspath(
#	os.path.join(os.path.dirname(__file__), "kmerGWAS.cpython-39-darwin.so"))
#print(libname)


#libk = cdll.LoadLibrary("build/lib.macosx-11-x86_64-3.9/kmerGWAS.cpython-39-darwin.so") 
#libk = cdll.LoadLibrary(libname) 
#from .kmerGWAS import kmer_string_alloc
import kmerGWAS

print(kmerGWAS)
print(kmerGWAS.__dict__)
#str_alloc = kmerGWAS.kmer_string_alloc

class KmerGWASDB(KmerDB):
	pass

class KmerGWASDBBuilder(KmerBuilder):
	def __init__(self, kmer_size):
		self.kmer_size = kmer_size;
		self.buffer = kmerGWAS.kmer_general.KmerGWAS_builder(kmer_size);
		print(self.buffer);

	def string_to_kmer(self, sequence):
		print(sequence)

	def kmer_to_string(self, sequence):
		raise NotImplementedError
