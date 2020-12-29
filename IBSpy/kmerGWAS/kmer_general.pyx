cimport kmer_general

from libc.string cimport strcpy, strlen

cdef extern from "stdint.h":
	ctypedef unsigned long long uint64_t
	ctypedef unsigned short uint8_t

ctypedef uint64_t kmerGWAS_kmer

cdef class KmerGWAS_builder:
	cdef short _kmer_size
	cdef char * _buffer
	cdef kmerGWAS_kmer _kmer
	
	def __cinit__(self, kmer_size):
		self._kmer_size = kmer_size
		self._buffer = kmer_general.kmer_string_alloc(kmer_size)


	def string_to_kmer(self, sequence):
		py_byte_string = sequence.encode('UTF-8')
		strcpy(self._buffer, py_byte_string)
		self._kmer = kmer_general.kmer2bits(self._buffer, self._kmer_size)
		return self._kmer

	def kmer_to_string(self, kmer ):
		self._kmer = kmer
		kmer_general.bits2kmer31(self._kmer, self._kmer_size, self._buffer)
		return self._buffer

	def __dealloc__(self):
		if self._buffer is not NULL:
			kmer_general.kmer_string_free(self._buffer)
