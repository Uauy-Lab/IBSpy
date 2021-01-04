
#ifndef NUCLEOTIDE_KMER_GWAS_H
#define NUCLEOTIDE_KMER_GWAS_H

typedef enum {
	Adenine = 0, Cytosine = 1, Guanine = 2, Thymine = 3, Undefined = 4
    //Never used for the binary representation of a kmer!.
} Nucleotide;

Nucleotide char_to_binary_nucleotide(char c);

#endif