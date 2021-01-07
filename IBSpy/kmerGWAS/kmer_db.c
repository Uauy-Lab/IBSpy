#include <assert.h>
#include <fcntl.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
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
	kgt->mem_table = 0;
	return kgt;
}


void kmer_gwas_table_free(KmerGwasTable ** kgt){
	
	if((*kgt)->mem_table != 0){
		kmer_gwas_table_mmap_close((*kgt));
	}
	if((*kgt)->kmer == NULL){
		free((*kgt)->kmer);
	}
	free(*kgt);
	*kgt = NULL;
}

void kmer_gwas_table_push_kmer(kmerGWAS_kmer kmer, KmerGwasTable * kgt){
	assert(kgt->number_of_kmers < kgt->capacity);
	kmerGWAS_kmer cannonical_kmer = kmer & KMERGWAS_ORIENTATION_MASK;

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

kmerGWAS_kmer * kmer_gwas_table_find(kmerGWAS_kmer kmer, KmerGwasTable * kgt){
	if(kgt->number_of_kmers == 0){
		return NULL;
	}
	//fprintf(stderr, "kmers_size %i %llu\n",sizeof(kmer), kgt->number_of_kmers);
	return bsearch(&kmer, kgt->kmer, kgt->number_of_kmers ,
		sizeof(kmer), kmer_compare_internal);
}


uint64_t kmer_gwas_sort_and_filter_unique(KmerGwasTable * kgt){
	uint64_t new_total_kmers = 0L;
	kmerGWAS_kmer last_kmer = 0L;
	kmerGWAS_kmer current_kmer = 0L;
	int cmp = 0;
	
	if(kgt->number_of_kmers < 2){
		return kgt->number_of_kmers;
	}
	qsort(kgt->kmer, kgt->number_of_kmers, sizeof(* kgt->kmer), kmer_compare_internal);
	last_kmer = kgt->kmer[0];
	uint64_t i;
	for(i=1; i < kgt->number_of_kmers; i++ ){
		current_kmer = kgt->kmer[i];
		cmp = kmer_compare_internal(&last_kmer, &current_kmer);
		assert(cmp >= 0);
		if(cmp == 0){
			last_kmer |= current_kmer;
		}else{
			kgt->kmer[new_total_kmers++] = last_kmer;
			last_kmer = current_kmer;
		}
	}
	kgt->kmer[new_total_kmers++] = last_kmer;
	kgt->number_of_kmers = new_total_kmers;
	kgt->kmer = realloc(kgt->kmer, sizeof(kmerGWAS_kmer) * kgt->number_of_kmers);
	kgt->capacity = new_total_kmers;
	return new_total_kmers;
}

void kmer_gwas_table_add_kmers_from_string(char * sequence, KmerGwasTable * kgt){
	uint8_t  kmer_size = kgt->kmer_size;
	uint64_t len = strlen(sequence);
	uint64_t total_kmers = len - kmer_size + 1;
	uint64_t current_stretch = 0;
	uint64_t i;
	assert(kgt->readonly == false);
	Nucleotide n;
	kgt->capacity    = kgt->number_of_kmers +  total_kmers;
	kgt->kmer = realloc(kgt->kmer, sizeof(kmerGWAS_kmer) * kgt->capacity );
	kmerGWAS_kmer tmp_kmer = 0L;

	for(i = 0; i < len; i++){
		n = char_to_binary_nucleotide(sequence[i]);
		if(n == Undefined){
			current_stretch = 0;
			tmp_kmer = 0L;
		}else{
			current_stretch++;
		}
		tmp_kmer = kmer_shift_and_insert(tmp_kmer, sequence[i], kmer_size);
		if(current_stretch >= kmer_size){
			kmer_gwas_table_push_kmer(tmp_kmer, kgt);
		}
	}
}

kmerGWAS_kmer kmer_gwas_table_get(uint64_t index, KmerGwasTable * kgt){
	assert(kgt->number_of_kmers > 0);
	return kgt->kmer[index];
}

void kmer_gwas_table_save(char * filename, KmerGwasTable * kgt){
	FILE *file = fopen(filename, "wb");
	fwrite(kgt->kmer, sizeof(kmerGWAS_kmer), kgt->number_of_kmers, file);
	fclose(file);
}

void kmer_gwas_table_mmap_read(char * file, KmerGwasTable * kgt ){
	struct stat s;
	uint64_t size;
	//FILE * fd = fopen (file, "rb");
	int fd = open (file, O_RDONLY);
	/* Get the size of the file. */
    int status = fstat (fd, & s);
    assert(status == 0);
    kmerGWAS_kmer kmer;
    size = s.st_size;

    //fprintf(stderr, "%llu, %lu, %i\n", size, sizeof(kmer), sizeof(s.st_size));

    kgt->kmer = (kmerGWAS_kmer *) mmap (NULL, size, PROT_READ, MAP_PRIVATE, fd, 0);
	kgt->number_of_kmers = size / sizeof(kmer);
	kgt->readonly = true;
	kgt->mem_table = fd;
	return;
}

void kmer_gwas_table_mmap_close(KmerGwasTable * kgt){
	struct stat s;
	int status = fstat (kgt->mem_table, & s);
	assert(status == 0);
	munmap(kgt->kmer, s.st_size);
	kgt->kmer = NULL;
	close(kgt->mem_table);
}

