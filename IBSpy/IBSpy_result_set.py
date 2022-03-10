import gc
from locale import normalize
import os
import string
import pandas as pd
from typing import List
from pyranges import PyRanges
from multiprocess import Pool
import pysam
from IBSpy.window_affinity_propagation import AffinityRunResults, cluster_by_haplotype, select_best_cluster

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

    def mapped_window_tabix(self, chromosome, start, end, assembly = None) -> PyRanges:
        if self.block_mapping is None:
            ret = self.values_matrix.values_for_region(chromosome, start, end)
        else: 
            targets = self.block_mapping.all_regions_for(chromosome, start, end, assembly=assembly)
            ret = list()
            for index, row in targets.as_df().iterrows():
                # print(index)
                # print(row)
                tmp =  self.values_matrix.values_for_region(row["Chromosome"], row["Start"], row["End"])
                ret.append(tmp)
        return pd.concat(ret).drop_duplicates()

    def _function_window_wrapper(self, start, function, chromosome, assembly=None ):
        end = start + self.options.affinity_window_size
        self.options.log(f"Running: {chromosome}:{start}-{end}")
        window = self.mapped_window(chromosome, start, end, assembly=assembly)
        return function(window), chromosome, start, end 

    def _function_window_wrapper_tabix(self, start, function, chromosome, assembly=None ):
        end = start + self.options.affinity_window_size
        self.options.log(f"Running: {chromosome}:{start}-{end}")
        window = self.mapped_window_tabix(chromosome, start, end, assembly=assembly)
        return function(window), chromosome, start, end 

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

    
    def map_window_iterator_tabix(self, chromosome= None, function= None):
        if function is None:
            function = lambda x: x
        tabixes: dict[str, pysam.TabixFile] = self.values_matrix.merged_values
        ret = list()
        for assembly, tabix in tabixes.items():
            chromosomes = tabix.contigs
            if chromosome not in chromosomes:
                next
            
            length = self.options.chromosome_length(assembly, chromosome)
            with Pool(self.options.pool_size) as p:
                wrapped_function = lambda x: self._function_window_wrapper_tabix(x, function, chromosome, assembly = assembly)
                res = p.imap(wrapped_function, range(0, length , self.options.affinity_window_size),self.options.chunks_in_pool)
                ret.extend(res)
        return ret 

    
    def path_affinity(self,chr=None):
        prefix = self.options.output_folder
        file = self.options.file_prefix
        if chr is not None:
            chr = f".{chr}."
        else:
            chr = ""
        return f"{prefix}/{file}{chr}.affinity.csv.gz"

    def run_affinity_propagation(self) -> pd.DataFrame :
        self.options.log("Preparing affinity propagation")
        max_missing = self.options.max_missing
        dampings = self.options.dampings
        iterations = self.options.iterations
        seed = self.options.seed
        self.options.log(f"Searching for {self.path_affinity}")
        if os.path.isfile(self.path_affinity()): 
            return pd.read_csv(self.path_affinity(),sep="\t")
        self.options.log("Not found, building")
        def run_single_run(gr ):
            runs = cluster_by_haplotype(gr, seed=seed, iterations=iterations, dampings=dampings, max_missing=max_missing)
            best = select_best_cluster(runs)
            return  best 

        ret = list()
        chromosomes = self.values_matrix.values_matrix.chromosomes
        for chr in chromosomes:
            print(f"Affi for {chr}")
            gc.collect()
            for best, chromosome, start, end in self.map_window_iterator(function= run_single_run, chromosome=chr):
                
                if best is not None:
                    best.chromosome = chromosome
                    best.start = start 
                    best.end = end
                    ret.append(best.as_df()) 
            df = pd.concat(ret)
            df.to_csv(self.path_affinity(chr=f"{chr}"), sep="\t",index=False)
        return df

    

    
    