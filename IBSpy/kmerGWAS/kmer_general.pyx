# queue.pyx

cimport kmer_general

cdef class KmerGWAS_builder:
	cdef short _kmer_size
	cdef char * _buffer;
	def __cinit__(self, kmer_size):
		self._kmer_size = kmer_size
		self._buffer = kmer_general.kmer_string_alloc(kmer_size)


	def string_to_kmer(self, sequence):
		kmer_general.kmer_string_alloc(sequence)
