#include<args.h>
/*
2021-6-18 10:09:52
author:fengcong@caa.cn
*/
struct args_st deal_args(int argc, char* argv[]){
    //2021-6-16 12:09:54
    cxxopts::Options options("IBScpp", ",original repository:https://github.com/Uauy-Lab/IBSpy");
    try
	{
		options.add_options()
			("k,kmers_db", "kmers_with_strand path(required)", cxxopts::value<string>())
			("r,reference", "reference file path(required)", cxxopts::value<string>())
			("p,threads","threads(default=1)",cxxopts::value<uint32_t>()->default_value("1"))
			("s,kmer_size","kmer size(default=31)",cxxopts::value<uint32_t>()->default_value("31"))
            ("w,window_size","window size(default=1000)",cxxopts::value<uint64_t>()->default_value("1000"))
			("h,help", "print help")
			;
		auto result = options.parse(argc, argv);
		if (result.count("help"))
		{
			cerr << options.help() << endl;
			exit(0);
		}
        //check required parametrs
		vector<string> required_parametrs({"kmers_db",  "reference"});
		for(size_t i=0; i<required_parametrs.size(); i++) {
			if(result.count(required_parametrs[i]) == 0) {
				cerr <<"ERROR: " <<required_parametrs[i] << " is a required parameter" << endl;
				cerr << options.help() << endl;
				exit(1);
			}
		}
		// Load parameters
		string kmers_database(result["kmers_db"].as<string>());
		string reference_path(result["reference"].as<string>());
		uint32_t threads(result["threads"].as<uint32_t>());
		uint32_t kmer_size(result["kmer_size"].as<uint32_t>());
        uint64_t windsize(result["window_size"].as<uint64_t>());


		// Check if all input files exist
		vector<string> required_files({kmers_database, reference_path, reference_path+".fai"});
		for(size_t i=0; i<required_files.size(); i++) {
			if(!is_file_exist(required_files[i])) {
				cerr << "ERROR: Couldn't find file: " << required_files[i] << endl;
				exit(1);
			}
		}
		//check if window_size > kmer_size
		if (kmer_size >= windsize){
			cerr << "ERROR: window_size should bigger than kmer_size!"<< endl;
			exit(1);
		}

        return { kmers_database,reference_path, threads , kmer_size, windsize};
        
    }catch (const cxxopts::OptionException& e)
	{
		cerr << "error parsing options: " << e.what() << endl;
		cerr << options.help() << endl;
		exit(1);
	}


}