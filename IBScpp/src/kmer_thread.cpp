#include<kmer_thread.h>
#include<algorithm>
// struct thread_data{
//     string   reference;
//     uint64_t joboffset;
//     uint64_t jobcount;
//     vector<window_info> &tasks;
//     vector<chromosome_inf> &chr_inf;
// };

/*
2021-6-18 10:09:52
author:fengcong@caa.cn
*/

#define KMER_SIZE 31
vector<string> seq2kmers(string seq){
    uint64_t total_kmers = seq.length()-KMER_SIZE+1;
    vector<string> kmers;
    transform(seq.begin(),seq.end(),seq.begin(),::toupper);
    for (uint64_t i = 0 ; i< total_kmers;i++){
        // uint64_t end = i+KMER_SIZE;
        string kmer = seq.substr(i, KMER_SIZE); 
        if ( kmer.find_first_not_of("AGCT")  == string::npos){
            kmers.push_back(kmer);
        }
        
    }

    return kmers;
}
// struct result{
//     string seqname;
//     uint64_t start;
//     uint64_t end;
//     uint64_t total_kmers;
//     uint64_t observed_kmers;
//     uint64_t variations;
//     uint64_t kmer_distance;
// };
void add_variation(result *res,  uint32_t gap_size, uint32_t kmer_size){
    if (gap_size == 0){
        return;
    }else{
        res->variations ++;
        int32_t kmer_distance = gap_size - (kmer_size - 1);
        if( kmer_distance <= 0){
            kmer_distance = abs(kmer_distance + 1);
        }
        res->kmer_distance += kmer_distance;
    }

}
    

void * thread_func(void *threadarg){
    struct thread_data *my_data = (struct thread_data *) threadarg;

    struct result *ret_data = (struct result *) calloc(my_data->jobcount,sizeof(result));
    uint64_t count = 0;
    for(uint64_t offset = my_data->joboffset;count <  my_data->jobcount; offset++,count++){
        window_info task = (*my_data->tasks)[offset];
        string seq = get_seq_from_reference(my_data->reference,my_data->chr_inf,task.chromosome,task.start,task.end);
        // string seq = get_seq_from_reference2(my_data->ref,task.chromosome,task.start,task.end);
        (ret_data+count)->seqname = task.chromosome;
        (ret_data+count)->start = task.start;
        (ret_data+count)->end = task.end;

        vector<string> kmers = seq2kmers(seq);
        (ret_data+count)->total_kmers = kmers.size();
        
        // cerr << "thread go" <<endl;

        uint32_t gap_size = 0;
        for (uint64_t i =0 ; i < kmers.size() ; i++){
            // if ( find_kmer_db_hash(my_data->kmers, kmer2bits(kmers[i]),KMER_SIZE )){
            //     (ret_data+count)->observed_kmers ++;
            //     add_variation(ret_data+count,gap_size,KMER_SIZE);
            //     gap_size = 0;
            // }else{
            //     gap_size ++;
            // }
            if ( find_kmer_db_vector(my_data->kmers, kmer2bits(kmers[i]),KMER_SIZE )){
                (ret_data+count)->observed_kmers ++;
                add_variation(ret_data+count,gap_size,KMER_SIZE);
                gap_size = 0;
            }else{
                gap_size ++;
            }
        }
        add_variation(ret_data+count,gap_size,KMER_SIZE);


    }
    return ret_data ; 
}