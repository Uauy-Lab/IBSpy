//Minimal kmer functions from https://github.com/voichek/kmersGWAS
//Original code released under GPL-3.0

#ifndef KMER_GENERAL_H
#define KMER_GENERAL_H
#include <stdint.h>
#include <stdbool.h>

typedef uint64_t kmerGWAS_kmer;

// Func: kmer reverse complement
inline kmerGWAS_kmer kmer_reverse_complement(kmerGWAS_kmer x, const uint8_t k_len) {
	x = ((x & 0xFFFFFFFF00000000) >> 32) | ((x & 0x00000000FFFFFFFF) << 32);
	x = ((x & 0xFFFF0000FFFF0000) >> 16) | ((x & 0x0000FFFF0000FFFF) << 16);
	x = ((x & 0xFF00FF00FF00FF00) >> 8)  | ((x & 0x00FF00FF00FF00FF) << 8);
	x = ((x & 0xF0F0F0F0F0F0F0F0) >> 4)  | ((x & 0x0F0F0F0F0F0F0F0F) << 4);
	x = ((x & 0xCCCCCCCCCCCCCCCC) >> 2)  | ((x & 0x3333333333333333) << 2);
	return (~x) >> (64 - k_len - k_len);
}

kmerGWAS_kmer kmer2bits(char * k, uint8_t kmer_size);
char * kmer_string_alloc(const uint8_t kmer_size);
void kmer_string_free(char * res);
void bits2kmer31(kmerGWAS_kmer w, const uint8_t kmer_size, char * res);
bool kmer_is_canonical(kmerGWAS_kmer b, const uint8_t kmer_size);
#endif

