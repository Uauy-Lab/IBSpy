import os
import string

from numpy import number


class IBSpyOptions:
    def __init__(self) -> None:
        self.window_size:int = 50000
        self.normalize:bool  = False
        self.filter_counts:float or number = None
        self.metadata:string = "metadata.tsv"
        self.samples:list=None 
        self.references:list=None 
        self.score:string="variations"
        self.window_size:int = 50000
        self.stat:string = 'mean'
        self.affinity_blocks:int = 20
        self.out_folder = "./out/"
        self.cache_tables = True
        try:
            os.mkdir(self.out_folder)
        except OSError as error:
            pass
        self.metadata_filename = os.path.basename(self.metadata)
        

    @property
    def file_prefix(self):
        arr = [str(self.window_size), "ws",
            str(self.score), "score",
            str(self.stat), "stat",
            str(self.window_size), "ws"]
        if self.filter is not None:
            arr.append(str(self.filter_counts), "filter")
        if self.normalize:
            arr.append("normalize")
        return "_".join(arr)
