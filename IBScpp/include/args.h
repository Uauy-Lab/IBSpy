#include <cxxopts/include/cxxopts.hpp>
#include<string>
#include<kmer_general.h>

using namespace std;

struct args_st
{
    string kmer_db_path; //data/to/path/kmer_with_strand
    string reference_file_path; // data/to/path/ref.fa
    uint32_t threads ; 
    uint32_t kmer_size;
    uint64_t window_size;
};

struct args_st deal_args(int argc, char* argv[]);

