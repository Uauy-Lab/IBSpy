#include<string>
#include<vector>
#include<fstream>
#include<iostream>  //log
#include <stdlib.h>
#include<string.h>

#include <map>
using namespace std;
struct window_info
{
    string chromosome;
    uint64_t start;
    uint64_t end;
    // string sequence;
};

struct chromosome_inf
{
    string chromosome;
    uint64_t len;
    uint64_t start_offset;
    uint32_t line_bases;
    uint32_t line_width;
};

vector<chromosome_inf> get_reference_inf(string reference);
vector<window_info> produce_tasks_via_refernce(vector<chromosome_inf> ref_inf , uint64_t window_size, uint64_t kmer_size);
string get_seq_from_reference(string reference ,const vector<chromosome_inf> * chr_inf, string chromosome, uint64_t start , uint64_t end);

map<string,string> * get_reference_map(string reference);
string get_seq_from_reference2(map<string,string> * reference , string chromosome, uint64_t start , uint64_t end);