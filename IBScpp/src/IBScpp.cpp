#include<IBScpp.h>
/*
2021-6-18 10:09:52
author:fengcong@caa.cn
This is the C / C ++ version of IBSpy, the original repository: https://github.com/Uauy-Lab/IBSpy
*/


int main(int argc, char* argv[]){
    //load and check parameters
    auto [kmer_db_path, reference_file_path , threads , kmer_size, window_size] = deal_args(argc,argv);

    //log
    cerr << "kmer_db_path: " << kmer_db_path << endl;
    cerr << "reference_file_path: " << reference_file_path << endl;
    cerr << "threads: " << threads << endl;
    cerr << "kmer_size: " << kmer_size << endl;
    cerr << "window_size: " << window_size << endl;

    //deal reference file
    vector<chromosome_inf> chr_inf =  get_reference_inf(reference_file_path);


    //get tasks list
    vector<window_info> tasks =  produce_tasks_via_refernce(chr_inf,window_size,kmer_size);


    //get reference map,this is an alternative.
    //No significant speed improvement, increased memory usage.so discard
    // map<string, string> * ref = get_reference_map(reference_file_path);


    //read kmers database
    //method 1, slow in loading process,more mem usage.
    // KmerUint64Hash1* kmers =  read_kmer_db_hash(  kmer_db_path);  
    //method 2, slow in loading process
    // KmersSet * kmers = load_kmer_raw_file2(kmer_db_path, 10000, false);
    //method 3, low mem usage,more faster 
    vector<uint64_t> *kmers = read_kmer_db_vector(  kmer_db_path);   

    //Assigning Tasks
    uint64_t tasks_num = tasks.size();
    cerr << "we total have " <<  tasks_num<< " tasks "<<endl;

    //calc tasknum for each thread
    cerr<<"calc tasknum"<<endl;
    vector<uint64_t> job_count(threads);
    for(size_t i=0;i<job_count.size()-1;i++){
        job_count[i]=(int)(tasks_num/threads);
        cerr << "\t" << job_count[i]; //
    }
    job_count[job_count.size()-1] = tasks_num-(job_count[0]*(job_count.size()-1));
    cerr << "\t" << job_count[job_count.size()-1] << "\n";// 

    

    //calc task offset for each thread
    cerr<<"calc offset"<<endl;
    vector<uint64_t> job_offset(threads);
    job_offset[0]=0;
    cerr << "\t" << job_offset[0] ; //
    for(size_t i=1;i<job_offset.size();i++){
        job_offset[i]=job_offset[i-1] + job_count[i-1];
        cerr << "\t" << job_offset[i]; //
    }
    cerr <<"\n"<<endl;

    // create the threads 
    cerr<<"create threads"<<endl;
    struct thread_data * td = (struct thread_data *)calloc(threads,sizeof(thread_data));
    pthread_t * tids=(pthread_t *)calloc(threads,sizeof(pthread_t));
    for(size_t i = 0; i < threads; ++i)
    {
        (td+i)->jobcount = job_count[i];
        (td+i)->reference =reference_file_path ;
        (td+i)->joboffset = job_offset[i];
        // vector<chromosome_inf> &tmp_chr_inf = chr_inf;
        (td+i)->chr_inf = &chr_inf;
        // vector<window_info> &tmp_tasks = tasks;
        (td+i)->tasks = &tasks;
        (td+i)->kmers = kmers;
        // (td+i)->ref = ref;
        //参数依次是：创建的线程id，线程参数，调用的函数，传入的函数参数
        int ret = pthread_create(tids+i, NULL, thread_func,(void *)(td+i));
        if (ret != 0)
        {
            cerr << "pthread_create error: error_code=" << ret << endl;
        }
    }


    //thread join
    cerr<<"waiting threads"<<endl;
    vector<void *> thread_res(threads);
    for(size_t i = 0; i < threads; ++i)
    {
        int ret = pthread_join ( *(tids+i), &thread_res[i] );
        if (ret != 0)
        {
            cerr << "pthread_join error: error_code=" << ret << endl;
        }

    }
    
    cout << "seqname\tstart\tend\ttotal_kmers\tobserved_kmers\tvariations\tkmer_distance" << endl;

    for(size_t i = 0; i < threads; ++i)
    {
        struct result * tmp_res = (struct result *)thread_res[i];
        for (size_t j=0;j < job_count[i];j++){
            cout << (tmp_res+j)->seqname << "\t";
            cout << (tmp_res+j)->start << "\t";
            cout << (tmp_res+j)->end << "\t";
            cout << (tmp_res+j)->total_kmers << "\t";
            cout << (tmp_res+j)->observed_kmers << "\t";
            cout << (tmp_res+j)->variations << "\t";
            cout << (tmp_res+j)->kmer_distance << endl;
        }

    }

    return 0;
}


