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
        pass

