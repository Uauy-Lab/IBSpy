import argparse
from ast import arguments
import configparser
import logging
import os
import string
import psutil
from numpy import number
import pandas as pd


class IBSpyOptions:
    def __init__(self) -> None:
        self.window_size:int = 50000
        self.normalize:bool  = False
        self.filter_counts:float or number = None
        self.metadata:string = "metadata.tsv"
        self._samples:list=None 
        self._references:list=None 
        self.score:string="variations"
        self.stat:string = 'mean'
        self.affinity_blocks:int = 20
        self.affinity_window_size:int = 1000000
        self._out_folder:string = "./out/"
        self.no_cache_tables:bool = False
        self.chromosome_mapping:string = None
        self._chromosomes:pd.DataFrame = None
        self.block_mapping:string = None
        self.pool_size: int = 4
        self.chunks_in_pool: int = 1
        self.log_level =  logging.INFO
        self._mapping_seqnames: dict = None
        self._logger: logging.Logger = None
        self._dampings=[0.5, 0.6, 0.7, 0.8, 0.9] #TODO: Add this as options
        self.max_missing=5 #TODO: Add this as options (check if this is actually called. )
        self.iterations=100 #TODO: Add this as options
        self.min_iterations=10
        self._seed = 42
        self._chromosome_lengths = None
        self.name:string = "all"
        self._chromosome_suffix: dict = None
        self.chromosome_suffix_path : string = None
        self.chromosome: string = None 
        self.tsv_output = "out.tsv"
   

    @property
    def dampings(self):
        return self._dampings
    @dampings.setter
    def dampings(self, value):
        tmp = parse_str_to_list(value)
        def dampling_validate(d):
            d = float(d)
            if d < 0.5 or d >= 1:
                raise ValueError(f"Dampings must be between >=0.5, <1 ('{d}' in [{value}])")
            return d
        self._dampings = list(map(dampling_validate, tmp))
        
    @property
    def chromosomes(self):
        return self._chromosomes

    @chromosomes.setter
    def chromosomes(self, value):
        self._chromosomes = pd.read_csv(value, sep="\t")

    def chromosome_length(self, assembly, chromosome) -> int:
        #TODO: Validate or find a way to get this from the rest of the metadata. 
        df = self._chromosomes
        try:
            ret = df[ (df["assembly"] == assembly) &  (df["chr"] == chromosome) ]["end"].array[0]
        except Exception as e:
            ret = None 
        if ret is None:
            print(df)
            raise Exception(f"Unable to find length for {assembly}, {chromosome}")

        return ret

    @property
    def seed(self):
        return self._seed 

    @seed.setter
    def seed(self, value):
        self.seed = value

    @property
    def samples(self):
        return self._samples
    
    @samples.setter
    def samples(self, value):
        self._samples = parse_str_to_list(value)

    @property
    def references(self):
        return self._references

    @references.setter
    def references(self, value):
        self._references = parse_str_to_list(value)

    @property
    def metadata_filename(self):
        os.path.basename(self.metadata)

    @property
    def output_folder(self):
        if not os.path.isdir(self._out_folder): 
            os.mkdir(self._out_folder)
        return self._out_folder
    
    @output_folder.setter
    def output_folder(self, value):
        self._out_folder = value
        
    @property
    def preferences(self):
        return vars(self)
    
    @preferences.setter
    def preferences(self, value):
        ints = ["window_size", "affinity_blocks", "affinity_window_size", "pool_size", "chunks_in_pool", "max_missing", "iterations"]
        bools = ["normalize", "cache_tables"]
        floats = ["filter_counts"]
        config = configparser.ConfigParser()
        config.read(value)
        for key in config['DEFAULT']:
            if key in ints:
                value = config["DEFAULT"].getint(key)
            elif key in bools:
                value = config["DEFAULT"].getboolean(key)
            elif key in floats:
                value = config["DEFAULT"].getfloat(key)
            else:
                value = config['DEFAULT'][key]
            setattr(self, key, value)


    @property
    def cache_folder(self):
        cache_folder =  f'{self.output_folder}/{self.cache_file_prefix}'
        if not os.path.isdir(cache_folder): 
            os.mkdir(cache_folder)
        return cache_folder

    def folder_for_reference(self, reference: string):
        out_f = f'{self.cache_folder}/{reference}'
        if not os.path.isdir(out_f): 
            os.mkdir(out_f)
        return out_f
            
    @property
    def cache_tables(self):
        return not self.no_cache_tables

    @property
    def file_prefix(self) -> string:
        arr = ["matrix",
            str(self.window_size), "ws",
            str(self.score), "score",
            str(self.stat), "stat",
            ]
        if self.filter_counts is not None:
            arr.append(str(self.filter_counts), "filter")
        if self.normalize:
            arr.append("normalize")
        return "_".join(arr)

    @property
    def cache_file_prefix(self) -> string:
        arr = ["combined",
            str(self.window_size), "ws",
            ]
        return "_".join(arr)

    @property
    def mapping_seqnames(self):
        if self._mapping_seqnames is not None:
            return self._mapping_seqnames
        self._mapping_seqnames = {}
        if self.chromosome_mapping is None:
            return self._mapping_seqnames
        map_df = pd.read_csv(self.chromosome_mapping, delimiter='\t')
        for index, row in map_df.iterrows():
            self._mapping_seqnames[row['original']] = row['mapping']
        return self._mapping_seqnames

    @property
    def chromosome_suffix(self) -> dict:
        if self._chromosome_suffix is not None:
            return self._chromosome_suffix
        self._chromosome_suffix = {}
        if self.chromosome_suffix_path is None: 
            return self._chromosome_suffix
        map_df = pd.read_csv(self.chromosome_suffix_path, sep="\t")
        for index, row in map_df.iterrows():
            self._chromosome_suffix[row["reference"]] = row["suffix"]
        return self._chromosome_suffix

    def chromosome_lengths(self):
        if self._chromosome_lengths is not None:
            return self._chromosome_lengths
        self._chromosome_lengths = {}
        

    @property
    def logger(self) -> logging.Logger:
        if self._logger is not None:
            return self._logger
        logging.basicConfig(filename=f'{self.output_folder}/{self.file_prefix}.log', level=logging.DEBUG, format='%(asctime)s:%(name)s:%(levelname)s:%(message)s')
        self._logger = logging.getLogger("IBSpy")
        self._logger.setLevel(self.log_level)
        return self._logger

    def log(self, text) -> None:
        self.logger.info(text)
        print(text)
        gb_mem = psutil.Process().memory_info().rss / (1024 * 1024 * 1024)
        gb_mem = round(gb_mem, 3)
        self.logger.debug(f"Mem: {gb_mem} GB")

    def assembly_for_chromosome(self, chromosome):
        df = self._chromosomes
        try:
            ret = df[df["chr"] == chromosome ]["assembly"].array[0]
        except Exception as e:
            ret = None 
        if ret is None:
            print(df)
            raise ValueError(f"Unable to find chromosome {chromosome}")
        return ret


    @property
    def region(self):
        reg:string = self._region.replace("'", "")
        fields_1 = reg.split(":")
        if len(fields_1) != 2:
            raise ValueError(f"Invalid region (must contain ':') {self._region}")
        fields_2 = fields_1[1].split("-")
        if len(fields_2) != 2:
            raise ValueError(f"Invalid region (must contain '-' on the second field) {self._region}")
        chr   = fields_1[0]
        start = int(fields_2[0])
        end   = int(fields_2[1])
        return chr, start, end, self.assembly_for_chromosome(chr)
        

    @region.setter
    def region(self, value):
        self._region: string = value

def parse_IBSpyOptions_arguments():
    ret = IBSpyOptions()
    parser = argparse.ArgumentParser()
    parser.add_argument("-w", "--window_size", default=50000,
		help="Window size to analyze", type=int)
    parser.add_argument("-n", "--normalize", default=False, action="store_true", 
        help="If present, the values are nomalised in a range from 0 to 1")
    parser.add_argument("-f", "--filter_counts", type=number, default=None,
        help="Minimum value from the results of IBSpy counts per window to be considered")
    parser.add_argument("-m", "--metadata", default=None, 
        help="Tab separated file with the paths to the results from IBSpy_windows with the following columns: reference, query, file and path")
    parser.add_argument("-s", "--samples", default = None, 
        help="File with the samples to include in the analysis (one per line)")
    parser.add_argument("-r", "--references", default= None, 
        help="File with the references to include in the analysis (one per line)" )
    parser.add_argument("-S", "--score", default="variations", 
        help="Column with the score to use from the output of IBSpy_window_count")
    parser.add_argument("-T", "--stat",default="mean" , 
        help="To calculate on the grouped windows [mean, median, variance or std]")
    parser.add_argument("-a", "--affinity_blocks", default=20, type=int, 
        help="Number of windows windows to group [deprectiated]")
    parser.add_argument("-A", "--affinity_window_size", default=1000000, type=int, 
        help="Size of the groups to merge for the affinity propagatin, in basepairs")
    parser.add_argument("-o", "--output_folder", default="./out/", 
        help="Folder with the output, including the cached files")
    # parser.add_argument("-C", "--no_cache_tables", default=False, action="store_true", 
    #     help="By default, the intermediate tables are stored. This flag dissables the cache ")
    parser.add_argument("-c", "--chromosome_mapping", default=None, 
        help="Path to tab separated file to rename the chromosome names. Columns are [original, mapping]")
    parser.add_argument("-B", "--block_mapping", default=None,
        help="Path with the file containing the coordinate mapping across pangenome")
    parser.add_argument("-p", "--pool_size", default=1, type=int, 
        help="Number of threads to use")
    parser.add_argument("-P", "--chunks_in_pool", default=1, type=int, 
        help="Number of workers per pool. ")
    parser.add_argument("-F", "--preferences", default=None, 
        help="Preferences file ")
    parser.add_argument("-l", "--chromosomes", default=None, 
        help="Tab separated file with the following columns [assembly, chr, start, end]. The chromsosome lenght is used to determine the last position when using the tabix tables")
    parser.add_argument("-N", "--name", default="all", help="Name for the analysis, when subsetting samples")
    parser.add_argument("-t", "--chromosome_suffix_path", default=None, help="File with the chromosome suffixes for each reference. Columns [reference, suffix]")
    parser.add_argument("-C", "--chromosome", default=None, help="Run the analysis for a chromosome")
    parser.add_argument("-d", "--dampings", default="0.5,0.6,0.7,0.8,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98", help="Damping values for the affinity propagation separade by commas. Values >=0.5 and <1")
    parser.add_argument("-M", "--max_missing", default=6, type=int, help="Out of the  min_iterations, how many runs are allowed to fail in the clustering")
    parser.add_argument("-I", "--iterations", default=100,type=int, help="Number of monter carlo runs for the affinity prpagation")
    parser.add_argument("-i", "--min_iterations", default=10, type=int, help = "Initial numbber of montecarlo tests. If all of them converge to the same clustering, don't run all the iterations" )
    parser.add_argument("-R", "--region", help="Region to extract (chr:start-end)")
    parser.add_argument("-O", "--tsv_output", default="out.csv", help="path for the CSV with the extracted daata")
    parser.parse_args(namespace=ret)
    return ret

def get_options() -> IBSpyOptions:
    ret = parse_IBSpyOptions_arguments()
    return ret

def parse_str_to_list(value):
    if type(value) is list:
        return value
    if os.path.isfile(value):
        return file_to_list(value)
    if type(value) is str:
        return value.split(",")
    raise TypeError(f"Invalid sample {str(value)}")

def file_to_list(value):
    ret = []
    with open(value, 'r') as file:
        for line in file:
            ret.append(line.strip())
    return ret