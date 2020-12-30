#include "nucleotide.h"

int char2bin[] = {
    Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined,
    Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined,
    Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined,
    Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined,
    Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined,
    Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined,
    Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined,
    Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined,
    Adenine, Undefined, Cytosine, Undefined, Undefined, Undefined, Guanine, Undefined,
    Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined,
    Thymine, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined,
    Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined,
    Adenine, Undefined, Cytosine, Undefined, Undefined, Undefined, Guanine, Undefined,
    Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined,
    Thymine, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined,
    Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined, Undefined,
};

Nucleotide char_to_binary_nucleotide(char c)
{
	return char2bin[(int) c];
}
