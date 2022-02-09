from locale import normalize
import string
import pandas as pd
from typing import List
from pyranges import PyRanges
from multiprocess import Pool

from .BlockMapping import BlockMapping
from .IBSpy_options import IBSpyOptions
from .IBSpy_values_matrix import IBSpyValuesMatrix
import warnings

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
        self._block_mapping = None

    @property
    def values_matrix(self):
        if self._values_matrix is not None:
            return self._values_matrix
        self._values_matrix = IBSpyValuesMatrix(self)
        return self._values_matrix

    @property
    def block_mapping(self):
        if self._block_mapping is not None:
            return self._block_mapping
        if self.options.block_mapping is None:
            return None
        self._block_mapping = BlockMapping(self.options.block_mapping)
        return self._block_mapping
        
    def values_matrix_seqname_iterator(self):
        for x in self.values_matrix.seqname_iterator():
            yield x

    def values_matrix_iterator(self):
        for x in self.values_matrix.matrix_iterator():
            yield x

    def to_csv(self,  *args, **kwargs):
        self.values_matrix.to_csv( *args, **kwargs)
            
    def mapped_window(self, chromosome, start, end, assembly = None) -> PyRanges:
        if self.block_mapping is None:
            return self.values_matrix.values_matrix[chromosome, start:end]
        targets = self.block_mapping.all_regions_for(chromosome, start, end, assembly=assembly)
        return self.values_matrix.values_matrix.intersect(targets)

    def _function_window_wrapper(self, start, function, chromosome, assembly=None ):
        window = self.mapped_window(chromosome, start, start + self.options.affinity_window_size, assembly=assembly)
        return function(window)

    def map_window_iterator(self, chromosome= None, function= None):
        if function is None:
            function = lambda x: x
        
        chromosome_lengths = self.values_matrix.chromosome_lengths
        if chromosome is not None:
            chromosome_lengths = {chromosome : chromosome_lengths[chromosome]}

        ret = list()
        for chromosome, length in chromosome_lengths.items():
            assembly = None #TODO: have the hash of chromosomes to assembly
            with Pool(self.options.pool_size) as p:
                wrapped_function = lambda x: self._function_window_wrapper(x, function, chromosome, assembly = assembly)
                res = p.imap(wrapped_function, range(0, length , self.options.affinity_window_size),self.options.chunks_in_pool)
                ret.extend(res)
        return ret 

    
    