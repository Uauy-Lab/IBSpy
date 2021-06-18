#include<pthread.h>
#include<fstream>
#include<iostream>
#include<string>
#include<reference.h>
#include<kmerdb.h>


using namespace std;

struct thread_data{
    string   reference;
    uint64_t joboffset;
    uint64_t jobcount;
    vector<window_info> *tasks;
    vector<chromosome_inf> *chr_inf;
    // KmerUint64Hash1 *kmers;
    vector<uint64_t> * kmers;
    // map<string,string> * ref;
    // KmersSet * kmers;
};

struct result{
    string seqname;
    uint64_t start;
    uint64_t end;
    uint64_t total_kmers;
    uint64_t observed_kmers;
    uint64_t variations;
    uint64_t kmer_distance;
};

void * thread_func(void *threadarg);