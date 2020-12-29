#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include "kmer_general_c.h"
#include "nucleotide.h"


/**
 * @brief  Given a k-mer representation in bits convert it to string
 * @param  64-bit k-mer representation 
 * @return  string with bp representation
 */
void bits2kmer31(kmerGWAS_kmer w, const uint8_t kmer_size, char * res) {
	const static char dict_bp[] = {'A','C','G','T'};
	const static kmerGWAS_kmer mask2bits = 0x0000000000000003;
	uint8_t i;
	if(!kmer_is_canonical(w, kmer_size))
		w = kmer_reverse_complement(w, kmer_size);
	for(i=0; i<kmer_size; i++) {
		res[kmer_size-i-1] = dict_bp[w & mask2bits];
		w = (w>>2);
	}
	res[kmer_size] = '\0';
	return;
}

bool kmer_is_canonical(kmerGWAS_kmer b, const uint8_t kmer_size){
	kmerGWAS_kmer flag1 = b & 0x4000000000000000;
	kmerGWAS_kmer flag2 = b & 0x8000000000000000;
	if(flag1)
		return true;
	if(flag2)
		return false;
	kmerGWAS_kmer bt = kmer_reverse_complement(b,  kmer_size);
	return b < bt;

}

char * kmer_string_alloc(const uint8_t kmer_size){
	char* res = malloc(sizeof(char) * (kmer_size+1));
	assert(res != NULL);
	return res;
}

void kmer_string_free(char * res){
	free(res);
}

kmerGWAS_kmer kmer2bits(char * kmer, uint8_t kmer_size) {
	kmerGWAS_kmer b  = 0;
	kmerGWAS_kmer bt = 0;
	
	for(uint64_t i=0; i < kmer_size; i++) {
		Nucleotide n = char_to_binary_nucleotide(kmer[kmer_size - i - 1]);
		b |= ((uint64_t) n << (i*2));
		// if(kmer_size > 10){
		// 	char * tmp = kmer_string_alloc(kmer_size);
		// 	bits2kmer31(b, kmer_size, tmp);
		// 	fprintf(stderr, "%i: %s %llu %llu\n", i, tmp, n, (n << (i*2)) );
		// 	kmer_string_free(tmp);
		// }
	}
	
	bt = kmer_reverse_complement(b,  kmer_size);
	if (bt < b)
		return bt | 0x8000000000000000;
	else
		return b  | 0x4000000000000000;
}