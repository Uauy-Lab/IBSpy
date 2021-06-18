#include<reference.h>
#include <algorithm>
/*
2021-6-18 10:09:52
author:fengcong@caa.cn
*/

vector<chromosome_inf> get_reference_inf(string reference){
    //2021-6-16 13:02:45
    vector<struct chromosome_inf> ref_inf;
    char buffer[1024]{0};
    // cout << reference+ string(".fai")<<endl;

    ifstream ref_index_handle (reference + string(".fai")); 
    if (ref_index_handle.is_open()){
        while (! ref_index_handle.eof() )
        {
            struct chromosome_inf tmpinf;

            ref_index_handle.getline (buffer,1024);  //读入每行
            if ( strlen(buffer) == 0 ){
                continue;
            }
            
            // C type string split 
            vector<string> res;
            char* temp = strtok(buffer, "\t");
            while(temp != NULL)
            {
                res.push_back(string(temp));
                temp = strtok(NULL, "\t");
            }

            tmpinf.chromosome = string(res[0]);
            tmpinf.len = strtoull( res[1].c_str(), NULL,10 );
            tmpinf.start_offset = strtoull( res[2].c_str(), NULL,10 );
            tmpinf.line_bases = atoi( res[3].c_str() );
            tmpinf.line_width = atoi( res[4].c_str() );

            ref_inf.push_back(tmpinf);

            // cout<< tmpinf.chromosome <<endl;
            // cout<< tmpinf.len <<endl;
            // cout<< tmpinf.start_offset <<endl;
            // cout<< tmpinf.line_bases <<endl;
            // cout<< tmpinf.line_width <<endl;
            
        }
    }
    
    ref_index_handle.close();


    return ref_inf;
}

vector<window_info> produce_tasks_via_refernce(vector<chromosome_inf> ref_inf , uint64_t window_size , uint64_t kmer_size){
    //2021-6-16 14:28:36
    vector<window_info> tasks;

    for (size_t i = 0 ; i < ref_inf.size() ; i++){
        chromosome_inf chr_inf = ref_inf[i];
        uint64_t pre_end = 0;
        for(size_t windstart = 0 ; windstart < chr_inf.len; windstart =  pre_end){
            window_info tmp_wd_inf ;
            tmp_wd_inf.chromosome = chr_inf.chromosome;

            if (pre_end == 0){
                tmp_wd_inf.start = 0;
            }else{
                tmp_wd_inf.start = pre_end-kmer_size;
            }

            tmp_wd_inf.end = tmp_wd_inf.start + window_size;
            if (tmp_wd_inf.end > chr_inf.len){
                tmp_wd_inf.end = chr_inf.len;
            }
            pre_end = tmp_wd_inf.end;
            tasks.push_back(tmp_wd_inf);
        }

        // window_info tmp_wd_inf ;
        // tmp_wd_inf.chromosome = chr_inf.chromosome;
        // tmp_wd_inf.start = pre_end-kmer_size;
        // tmp_wd_inf.end = chr_inf.len;
        // tasks.push_back(tmp_wd_inf);
        
    }

    return tasks;

}

string get_seq_from_reference(string reference ,const vector<chromosome_inf> * chr_inf, string chromosome, uint64_t start , uint64_t end){
    uint64_t offset_start = 0;
    uint64_t offset_end = 0;
    for (size_t i = 0 ; i < (*chr_inf).size() ;i++){
        // uint64_t start_offset;
        // uint32_t line_bases;
        // uint32_t line_width;
        if ( strcmp ( (*chr_inf)[i].chromosome.c_str(), chromosome.c_str()) == 0){
            offset_start = (*chr_inf)[i].start_offset + start + start/(*chr_inf)[i].line_bases;
            offset_end = (*chr_inf)[i].start_offset + end + end/(*chr_inf)[i].line_bases;
        }
    }
    
    ifstream ref_handle (reference );
    if (ref_handle.is_open()){
        char * buffer =  (char *)calloc(offset_end-offset_start+1,sizeof(char));
        ref_handle.seekg(offset_start,ios::beg);  
        ref_handle.read(buffer,offset_end-offset_start);

        // std:remove(buffer, buffer+strlen(buffer), '\n') == '\0';
        string buff_string(buffer);
        buff_string.erase(std::remove(buff_string.begin(), buff_string.end(), '\n'), buff_string.end());
        ref_handle.close();
        return buff_string;
    }else{
        cerr << "cant open refernce file"<< endl;
        exit(-1);
    }
}

int startsWith(string s, string sub){
        return s.find(sub)==0?1:0;
}

map<string,string> * get_reference_map(string reference){
    map<string,string> * ref = new map<string,string>();
    ifstream ref_handle (reference );
    string tmp_str="";
    char * buffer =  (char *)calloc(1024,sizeof(char));
    string current_chr = "";
    bool first = true;
    if (ref_handle.is_open()){
        while(1){
            ref_handle.getline(buffer,1024);
            if(!ref_handle.eof()){
                string buff_string(buffer);
                buff_string.erase(std::remove(buff_string.begin(), buff_string.end(), '\n'), buff_string.end());
                if (startsWith(buff_string,string(">"))) {
                    if (!first){
                        ref->insert ( std::pair<string,string>(current_chr,tmp_str) );
                    }else{
                        first = false;
                    }
                    current_chr =  buff_string.substr(1, buff_string.length() - 1);
                    tmp_str = "";
                }else{
                    tmp_str += buff_string;
                }
            }else{
                break;
            }
        }
        ref->insert ( std::pair<string,string>(current_chr,tmp_str) );
        ref_handle.close();

    }

    return ref;
}


string get_seq_from_reference2(map<string,string> * reference , string chromosome, uint64_t start , uint64_t end){
    

    return (*reference)[chromosome].substr(start,end-start);
     
}
