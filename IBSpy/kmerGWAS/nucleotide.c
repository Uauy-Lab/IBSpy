#include "nucleotide.h"
Nucleotide char_to_binary_nucleotide(char c)
{
	switch (c) {
		case 'A':
			return Adenine;
		case 'C':
			return Cytosine;
		case 'G':
			return Guanine;
		case 'T':
			return Thymine;
		case 'a':
			return Adenine;
		case 'c':
			return Cytosine;
		case 'g':
			return Guanine;
		case 't':
			return Thymine;
		default:
			return Undefined;
	}
}