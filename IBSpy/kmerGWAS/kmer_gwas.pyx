cimport kmer_gwas

from libc.string cimport strcpy, strlen

cdef extern from "stdint.h":
	ctypedef unsigned long long uint64_t
	ctypedef unsigned short uint8_t

ctypedef uint64_t kmerGWAS_kmer

cdef extern from "kmer_db.h":
	struct KmerGwasTable:
		uint64_t number_of_kmers
		uint64_t capacity
		kmerGWAS_kmer * kmer
		uint8_t kmer_size
		uint8_t readonly


cdef class KmerGWAS_builder:
	cdef short _kmer_size
	cdef char * _buffer
	cdef kmerGWAS_kmer _kmer
	
	def __cinit__(self, kmer_size):
		self._kmer_size = kmer_size
		self._buffer = kmer_gwas.kmer_string_alloc(kmer_size)

	def string_to_kmer(self, sequence):
		py_byte_string = sequence.encode('UTF-8')
		strcpy(self._buffer, py_byte_string)
		self._kmer = kmer_gwas.kmer2bits(self._buffer, self._kmer_size)
		return self._kmer

	def kmer_to_string(self, kmer ):
		self._kmer = kmer
		kmer_gwas.bits2kmer31(self._kmer, self._kmer_size, self._buffer)
		return self._buffer

	def compare(self, a, b):
		return kmer_gwas.kmer_compare(a,b)

	def __dealloc__(self):
		if self._buffer is not NULL:
			kmer_gwas.kmer_string_free(self._buffer)

cdef class KmerGWAS_database: 
	cdef KmerGwasTable * _kgt

	def __cinit__(self, kmer_size): 
		self._kgt = kmer_gwas.kmer_gwas_table_new(kmer_size)

	def kmer_size(self):
		return self._kgt.kmer_size

	def add_kmers(self, sequence):
		py_byte_string = sequence.encode('UTF-8')
		kmer_gwas.kmer_gwas_table_add_kmers_from_string(py_byte_string, self._kgt)

	def read_mmap(self, filename):
		py_byte_string = filename.encode('UTF-8')
		kmer_gwas.kmer_gwas_table_mmap_read(py_byte_string, self._kgt)

	def save(self, filename):
		py_byte_string = filename.encode('UTF-8')
		kmer_gwas.kmer_gwas_table_save(py_byte_string, self._kgt)

	def sort_unique(self):
		return kmer_gwas.kmer_gwas_sort_and_filter_unique(self._kgt)

	def __contains__(self, kmer):
		ret = kmer_gwas.kmer_gwas_table_find(kmer, self._kgt) 
		return ret is not NULL

	def __len__(self):
		return self._kgt.number_of_kmers

	def  __getitem__(self, index):
		return kmer_gwas.kmer_gwas_table_get(index, self._kgt)

	def __dealloc__(self):
		kmer_gwas.kmer_gwas_table_free(&self._kgt)


