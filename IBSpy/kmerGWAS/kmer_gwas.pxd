cdef extern from "stdint.h":
	ctypedef unsigned long long uint64_t
	ctypedef unsigned short uint8_t

cdef extern from "kmer_general_c.h":
	ctypedef uint64_t kmerGWAS_kmer

	kmerGWAS_kmer kmer2bits(char * k, uint8_t kmer_size)
	char * kmer_string_alloc(const uint8_t kmer_size)
	void kmer_string_free(char * res)
	void bits2kmer31(kmerGWAS_kmer kmer, uint8_t kmer_size, char * res)
