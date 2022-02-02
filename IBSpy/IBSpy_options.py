import string


class IBSpyOptions:
    def __init__(self) -> None:
        self.window_size = 50000
        self.normalize:bool  = False
        self.filter_counts = None
        self.metadata = "metadata.tsv"
        self.samples:list=None 
        self.references:list=None 
        self.score:string="variations"
        self.window_size:int = 50000
        self.stat:string = 'mean'
        pass

