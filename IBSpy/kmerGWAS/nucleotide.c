#include "nucleotide.h"
#include <stdio.h>

const int nucelotide_char2bin_table[] = {
    Undefined, Undefined, Undefined, Undefined, Undefined, 
    Undefined, Undefined, Undefined, Undefined, Undefined, 
    Undefined, Undefined, Undefined, Undefined, Undefined, 
    Undefined, Undefined, Undefined, Undefined, Undefined, 
    Undefined, Undefined, Undefined, Undefined, Undefined,
    Undefined, Undefined, Undefined, Undefined, Undefined, 
    Undefined, Undefined, Undefined, Undefined, Undefined, 
    Undefined, Undefined, Undefined, Undefined, Undefined, 
    Undefined, Undefined, Undefined, Undefined, Undefined, 
    Undefined, Undefined, Undefined, Undefined, Undefined, 
    Undefined, Undefined, Undefined, Undefined, Undefined, 
    Undefined, Undefined, Undefined, Undefined, Undefined, 
    Undefined, Undefined, Undefined, Undefined, Undefined,
    Adenine,   Undefined, Cytosine,  Undefined, Undefined, 
    Undefined, Guanine,   Undefined, Undefined, Undefined, 
    Undefined, Undefined, Undefined, Undefined, Undefined, 
    Undefined, Undefined, Undefined, Undefined, Thymine, 
    Undefined, Undefined, Undefined, Undefined, Undefined, 
    Undefined, Undefined, Undefined, Undefined, Undefined, 
    Undefined, Undefined, Adenine,   Undefined, Cytosine,  
    Undefined, Undefined, Undefined, Guanine,   Undefined,
    Undefined, Undefined, Undefined, Undefined, Undefined, 
    Undefined, Undefined, Undefined, Undefined, Undefined, 
    Undefined, Thymine,   Undefined, Undefined, Undefined, 
    Undefined, Undefined, Undefined, Undefined, Undefined, 
    Undefined, Undefined, Undefined
};


int nucleotide_table_printed = 0;

void nucleotide_print_table(){
    int i;
    for(i = 0; i < 128; i++){
        fprintf(stderr, "%i: %c %i\n",i, i, nucelotide_char2bin_table[(int) i]) ; 
    }  
}

Nucleotide char_to_binary_nucleotide(char c)
{
    return nucelotide_char2bin_table[(int) c];
}
