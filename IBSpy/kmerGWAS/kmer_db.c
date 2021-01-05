#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "kmer_general.h"
#include "nucleotide.h"
#include "kmer_db.h"


KmerGwasTable * kmer_gwas_table_new(uint8_t kmer_size){
	KmerGwasTable * kgt = malloc(sizeof(KmerGwasTable));
	assert( kgt != NULL);
	kgt->number_of_kmers = 0;
	kgt->capacity = 0;
	kgt->kmer_size = kmer_size;
	kgt->readonly = false;
	kgt->kmer = NULL;
	return kgt;
}


void kmer_gwas_table_free(KmerGwasTable ** kgt){
	if((*kgt)->kmer == NULL){
		free((*kgt)->kmer);
	}
	free(*kgt);
	*kgt = NULL;
}

void kmer_gwas_table_push_kmer(kmerGWAS_kmer kmer, KmerGwasTable * kgt){
	assert(kgt->number_of_kmers < kgt->capacity);
	kmerGWAS_kmer cannonical_kmer = kmer & KMERGWAS_ORIENTATION_MASK;;
	if(kmer_is_canonical(cannonical_kmer, kgt->kmer_size)){
		cannonical_kmer |= KMERGWAS_FORWARD;
	}else{
		cannonical_kmer = kmer_reverse_complement(cannonical_kmer, kgt->kmer_size);
		cannonical_kmer |= KMERGWAS_REVERSE;
	}
	//fprintf(stderr, "%llu: %llu %llu\n", kgt->number_of_kmers, kmer, cannonical_kmer );
	kgt->kmer[kgt->number_of_kmers] = cannonical_kmer; 
	kgt->number_of_kmers++;
}

void kmer_gwas_table_add_kmers_from_string(char * sequence, KmerGwasTable * kgt){
	uint8_t  kmer_size = kgt->kmer_size;
	uint64_t len = strlen(sequence);
	uint64_t total_kmers = len - kmer_size + 1;
	uint64_t current_stretch = 0;
	Nucleotide n;
	kgt->capacity    = kgt->number_of_kmers +  total_kmers;
	kgt->kmer = realloc(kgt->kmer, sizeof(kmerGWAS_kmer) * kgt->capacity );
	kmerGWAS_kmer tmp_kmer = 0L;
	for(uint64_t i = 0; i < len; i++){
		n = char_to_binary_nucleotide(sequence[i]);
		if(n == Undefined){
			current_stretch = 0;
			tmp_kmer = 0L;
		}else{
			current_stretch++;
		}
		fprintf(stderr, "%llu %llu\n",i, tmp_kmer);
		tmp_kmer = kmer_shift_and_insert(tmp_kmer, sequence[i], kmer_size);
		if(current_stretch >= kmer_size){
			kmer_gwas_table_push_kmer(tmp_kmer, kgt);
		}
	}
	kgt->kmer = realloc(kgt->kmer, sizeof(kmerGWAS_kmer) * kgt->number_of_kmers);
}

kmerGWAS_kmer kmer_gwas_table_get(uint64_t index, KmerGwasTable * kgt){
	assert(kgt->number_of_kmers > 0);
	return kgt->kmer[index];
}