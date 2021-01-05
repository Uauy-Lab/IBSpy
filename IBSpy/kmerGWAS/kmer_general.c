#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include "kmer_general.h"
#include "nucleotide.h"


/**
 * @brief  Given a k-mer representation in bits convert it to string
 * @param  64-bit k-mer representation 
 * @return  string with bp representation
 */
void bits2kmer31(kmerGWAS_kmer w, const uint8_t kmer_size, char * res) {
	const static char dict_bp[] = {'A','C','G','T'};
	const static kmerGWAS_kmer mask2bits = KMERGWAS_MASK;
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
	kmerGWAS_kmer flag_forward = b & KMERGWAS_FORWARD;
	kmerGWAS_kmer flag_reverse = b & KMERGWAS_REVERSE;
	if(flag_forward)
		return true;
	if(flag_reverse)
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

kmerGWAS_kmer kmer_shift_and_insert(kmerGWAS_kmer kmer, char base, uint8_t kmer_size){
	kmerGWAS_kmer n = (kmerGWAS_kmer) char_to_binary_nucleotide(base);
	kmerGWAS_kmer kmer2 = kmer << 2;
	kmerGWAS_kmer wall  = (kmerGWAS_kmer)  KMERGWAS_MASK << kmer_size * 2;
	kmer2 |= (kmerGWAS_kmer) n;
	kmer2 &= ~wall;
	//fprintf(stderr, "shifting %llu %llu %c %llu %llu\n", kmer, kmer2, base, n, wall);
	return kmer2 ;
}

kmerGWAS_kmer kmer2bits(char * kmer, uint8_t kmer_size) {
	kmerGWAS_kmer b  = 0;
	kmerGWAS_kmer bt = 0;
	for(int i=0; i < kmer_size; i++) {
		b = kmer_shift_and_insert(b, kmer[i], kmer_size);
	}
	b &= KMERGWAS_ORIENTATION_MASK;
	// for(uint64_t i=0; i < kmer_size; i++) {
	// 	Nucleotide n = char_to_binary_nucleotide(kmer[kmer_size - i - 1]);
	// 	b |= ((uint64_t) n << (i*2));
	// }
	
	bt = kmer_reverse_complement(b,  kmer_size);
	if (bt < b)
		return bt | KMERGWAS_REVERSE;
	else
		return b  | KMERGWAS_FORWARD;
}

int kmer_compare_internal (const void * a, const void * b) {
	kmerGWAS_kmer ak = *(kmerGWAS_kmer*)a;
	kmerGWAS_kmer bk = *(kmerGWAS_kmer*)b; 
	ak &= KMERGWAS_ORIENTATION_MASK;
	bk &= KMERGWAS_ORIENTATION_MASK;
	return (ak > bk) - (ak < bk);
	//return ( ak - bk );
}

int kmer_compare (kmerGWAS_kmer  a, kmerGWAS_kmer  b){
	return kmer_compare_internal( &a, &b );
} 

