import argparse
import os
import string

from numpy import number
import pandas as pd


class IBSpyOptions:
    def __init__(self) -> None:
        self.window_size:int = 50000
        self.normalize:bool  = False
        self.filter_counts:float or number = None
        self.metadata:string = "metadata.tsv"
        self.samples:list=None 
        self.references:list=None 
        self.score:string="variations"
        self.stat:string = 'mean'
        self.affinity_blocks:int = 20
        self.affinity_window_size:int = 1000000
        self.out_folder:string = "./out/"
        self.cache_tables:bool = True
        self.chromosome_mapping:string = None
        self.block_mapping:string = None
        self.pool_size: int = 1
        # self.chunks_in_pool = 100
        self._mapping_seqnames = None
        try:
            os.mkdir(self.out_folder)
        except OSError as error:
            pass
        self.metadata_filename = os.path.basename(self.metadata)
        

    @property
    def file_prefix(self) -> string:
        arr = [str(self.window_size), "ws",
            str(self.score), "score",
            str(self.stat), "stat",
            str(self.window_size), "ws"]
        if self.filter is not None:
            arr.append(str(self.filter_counts), "filter")
        if self.normalize:
            arr.append("normalize")
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



def parse_IBSpyOptions_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-w", "--window_size", default=50000,
		help="window size to analyze", type=int)
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
        help="to calculate on the grouped windows [mean, median, variance or std]")
    


    pass