from locale import normalize
import string
import pandas as pd
from typing import List

from IBSpy import IBSpyResults
from IBSpy.IBSpy_options import IBSpyOptions

class IBSpyResultsSet:
    def __init__(self, options: IBSpyOptions):
        
        self.samples_df = pd.read_csv(options.metadata, delimiter='\t')
        if options.references is not None:
            self.samples_df =  self.samples_df[ self.samples_df['reference'].isin(options.references)]
        if options.samples is not None:
            self.samples_df = self.samples_df[self.samples_df['query'].isin(options.samples)]
        self.options = options
        self.references = self.samples_df["reference"].unique()
        self.samples = self.samples_df["query"].unique()
        self._values_matrix  = None

    def values_matrix_for_reference(self, reference:string):
        samples = self.samples_df[self.samples_df['reference'] == reference]
        values_matrix = pd.DataFrame()
        for index, row in samples.iterrows():
            path=row['path'] + "/" + row['file']
            sample_name = row['query']
            ibr = IBSpyResults(path, self.options)
            df: pd.DataFrame = ibr.count_by_windows()
            values_matrix[sample_name] = df[self.options.stat]
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

    
    def values_matrix_seqname_iterator(self):
        mat = self.values_matrix
        seqnames = mat['seqname'].unique()
        for seq in seqnames:
            yield mat[mat['seqname'] == seq]

    def values_matrix_iterator(self):
        for seq in self.values_matrix_seqname_iterator():
            # total_windows = len(seq) / self.options.affinity_blocks
            for x in range(0, len(seq) ,self.options.affinity_blocks):
                yield seq.iloc[x:x+self.options.affinity_blocks]



