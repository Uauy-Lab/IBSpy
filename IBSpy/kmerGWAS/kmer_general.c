#include "kmer_general.h"
#include <assert.h>


/**
 * @brief  Given a k-mer representation in bits convert it to string
 * @param  64-bit k-mer representation 
 * @return  string with bp representation
 */
void bits2kmer31(kmerGWAS_kmer w, const uint8_t kmer_size, char * res) {
	const static char dict_bp[] = {'A','C','G','T'};
	const static kmerGWAS_kmer mask2bits = 0x0000000000000003;
 	uint8_t i;
	for(i=0; i<kmer_size; i++) {
		res[kmer_size-i-1] = dict_bp[w & mask2bits];
		w = (w>>2);
	}
	res[kmer_size] = NULL;
	return;
}

char * kmer_string_alloc(const uint8_t kmer_size){
	char* res = malloc(sizeof(char) * (k+1));
	assert(res != NULL);
	return res;
}

void kmer_string_free(char * res){
	free(res);
}

kmerGWAS_kmer kmer2bits(char * kmer, uint8_t kmer_size) {
	kmerGWAS_kmer b  = 0;
	kmerGWAS_kmer bt = 0;
	kmerGWAS_kmer cur_dubit;
	for(size_t i=0; i < kmer_size; i++) {
		Nucleotide n = char_to_binary_nucleotide(kmer[kmer_size - i - 1])
		
		b |= (n << (i*2));
	}
	//kmerGWAS_kmer bt = kmer2bits(b, k.size());  
	kmer_reverse_complement(kmerGWAS_kmer x, const uint8_t k_len)
	if (bt < b)
		return bt;
	else
		return b;
}