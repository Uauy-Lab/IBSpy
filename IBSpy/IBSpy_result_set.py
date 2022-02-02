from locale import normalize
import string
import pandas as pd
from typing import List

from IBSpy import IBSpyResults

class IBSpyResultsSet:
    def __init__(self, filename="metadata.tsv", samples:list=None, 
        references:list=None, score:string="variations", normalize:bool=False,
        window_size:int = 50000, stat:string = 'mean'):
        
        self.samples_df = pd.read_csv(filename, delimiter='\t')
        if references is not None:
            self.samples_df =  self.samples_df[ self.samples_df['reference'].isin(references)]
        if samples is not None:
            self.samples_df = self.samples_df[self.samples_df['query'].isin(samples)]
        self.score:string  = score
        self.normalize:bool = normalize
        self.window_size:int = window_size
        self.stat:string = stat
        self.references = self.samples_df["reference"].unique()
        self.samples = self.samples_df["query"].unique()
        self._values_matrix  = None

    def values_matrix_for_reference(self, reference:string):
        samples = self.samples_df[self.samples_df['reference'] == reference]
        values_matrix = pd.DataFrame()
        for index, row in samples.iterrows():
            path=row['path'] + "/" + row['file']
            sample_name = row['query']
            ibr = IBSpyResults(path, self.window_size, score = self.score, normalize = self.normalize)
            df: pd.DataFrame = ibr.count_by_windows()
            values_matrix[sample_name] = df[self.stat]
        names = df[['seqname','start','end']]
        values_matrix = pd.concat([names,values_matrix], axis=1)
        return values_matrix


    @property
    def values_matrix(self):
        if self._values_matrix is not None:
            return self._values_matrix
        dfs = []
        for ref in self.references:
            df = self.values_matrix_for_reference(ref)
            dfs.append(df)
        self._values_matrix = pd.concat(dfs, join="inner")
        return self._values_matrix