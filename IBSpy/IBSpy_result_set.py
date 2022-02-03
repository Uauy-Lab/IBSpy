from locale import normalize
import string
import pandas as pd
from typing import List
from .IBSpy_options import IBSpyOptions
from .IBSpy_values_matrix import IBSpyValuesMatrix
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

   

    @property
    def values_matrix(self):
        if self._values_matrix is not None:
            return self._values_matrix
        self._values_matrix = IBSpyValuesMatrix(self)
        return self._values_matrix
        
    def values_matrix_seqname_iterator(self):
        for x in self.values_matrix.seqname_iterator():
            yield x

    def values_matrix_iterator(self):
        for x in self.values_matrix.matrix_iterator():
            yield x

    def to_csv(self,  *args, **kwargs):
        self.values_matrix.to_csv( *args, **kwargs)
            


