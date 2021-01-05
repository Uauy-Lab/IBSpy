
#ifndef KMER_GWAS_DB_H
#define KMER_GWAS_DB_H
#include <stdint.h>
#include <stdbool.h>

typedef struct KmerGwasTable{
    uint64_t number_of_kmers;
    uint64_t capacity;
    kmerGWAS_kmer * kmer;
    uint8_t kmer_size;
    uint8_t readonly;
} KmerGwasTable;

KmerGwasTable * kmer_gwas_table_new(uint8_t kmer_size);
void kmer_gwas_table_free(KmerGwasTable * * kgt);
void kmer_gwas_table_add_kmers_from_string(char * sequence, KmerGwasTable * kgt);
kmerGWAS_kmer kmer_gwas_table_get(uint64_t index, KmerGwasTable * kgt);

#endif
