import gc
from locale import normalize
import os
import string
from numpy import dtype
import pandas as pd
from typing import List
from pyranges import PyRanges
from multiprocess import Pool
import pysam
from IBSpy.window_affinity_propagation import AffinityRunResults, cluster_by_haplotype, select_best_cluster
import dill
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
        ret = []
        if self.block_mapping is None:
            ret.append(self.values_matrix.values_for_region(chromosome, start, end))
        else: 
            targets = self.block_mapping.all_regions_for(chromosome, start, end, assembly=assembly)
            
            # print("[mapped_window_tabix] searching for mapped")
            for index, row in targets.as_df().iterrows():
                tmp =  self.values_matrix.values_for_region(row["Chromosome"], row["Start"], row["End"])
                if tmp is not None:
                    ret.append(tmp)
                    # print("++++++++++")
                    # print(tmp)
        ret = pd.concat(ret, ignore_index=True, axis=0, join="outer")
        # if ret is not None:
        # #     ret[['Start']] = ret[['Start']].apply(pd.to_numeric) 
        # #     ret[['End']] = ret[['End']].apply(pd.to_numeric) 
        # #     ret[[self.options.score]] = ret[[self.options.score]].apply(pd.to_numeric) 
        # #     ret[['mean']] = ret[['mean']].apply(pd.to_numeric) 
        # #     ret[['median']] = ret[['median']].apply(pd.to_numeric) 
        #     ret[['variance']] = ret[['variance']].apply(pd.to_numeric) 
        #     ret[['std']] = ret[['std']].apply(pd.to_numeric) 
        #     ret['sample'] = ret['sample'].astype('string') 
        #     ret['Chromosome'] = ret['Chromosome'].astype('string') 
        #     ret=ret.convert_dtypes()

        # print("~~~~~~")
        # print(ret.dtypes)
        # print(ret)
        return ret

    def _function_window_wrapper(self, start, function, chromosome, assembly=None ):
        end = start + self.options.affinity_window_size
        self.options.log(f"Running: {chromosome}:{start}-{end}")
        window = self.mapped_window(chromosome, start, end, assembly=assembly)
        return function(window), chromosome, start, end 

    def _function_window_wrapper_tabix(self, start, function, chromosome, assembly=None ):
        end = start + self.options.affinity_window_size
        self.options.log(f"[_function_window_wrapper_tabix] Running: {chromosome}:{start}-{end}")
        window = self.mapped_window_tabix(chromosome, start, end, assembly=assembly)
        # print(f"[_function_window_wrapper_tabix]About to run region: {chromosome}:{start}-{end}")
        # print(window)
        window = window.drop_duplicates(keep="last")
        m=window.pivot(index=["Chromosome", "Start", "End"], columns="sample", values="variations")
        m.reset_index(inplace=True)
        m.rename_axis(None, axis=1, inplace=True)
        m.columns
        # print(m)
        return function(m), chromosome, start, end 

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
        #print(tabixes)
        ret = list()
        for assembly, tabix in tabixes.items():
            chromosomes = tabix.contigs
            if chromosome not in chromosomes:
                continue
            # print(chromosomes)
            self.options.log(f"[map_window_iterator_tabix]About to get length for chromosome {chromosome}, {assembly}")
            length = self.options.chromosome_length(assembly, chromosome)
            # print(f"Lenght {length}")
            # print(f"Window {self.options.affinity_window_size}")
            # print(type(chromosome))
            # print(type(assembly))
            # print(f"chunks in pool: {self.options.chunks_in_pool}")
            
            with Pool(self.options.pool_size) as p:
                wrapped_function = lambda x: self._function_window_wrapper_tabix(x, function, chromosome, assembly = assembly)
                res = p.imap(wrapped_function, range(0, length , self.options.affinity_window_size),self.options.chunks_in_pool)
                #res = map(wrapped_function,  range(0, length , self.options.affinity_window_size))
                print("About to merge")
                newlist = [x for x in res]
                #print(list(res))
                ret.extend(newlist)
        #print(ret)
        return ret 

    
    def path_affinity(self,chr=None):
        prefix = self.options.output_folder
        file = self.options.file_prefix
        name = self.options.name
        if chr is not None:
            chr = f".{chr}."
        else:
            chr = ""
        return f"{prefix}/{file}{chr}.{name}.affinity.csv.gz"

    def run_affinity_propagation(self) -> pd.DataFrame :
        self.options.log("Preparing affinity propagation")
        max_missing = self.options.max_missing
        dampings = self.options.dampings
        iterations = self.options.iterations
        seed = self.options.seed
        # # self.options.log(f"Searching for {self.path_affinity(chr=chr)}")
        # # if os.path.isfile(self.path_affinity(chr=f"{chr}")): 
        # #     return pd.read_csv(self.path_affinity(chr=f"{chr}"),sep="\t")
        # self.options.log("Not found, building")
        def run_single_run(gr ):
            runs = cluster_by_haplotype(gr, seed=seed, iterations=iterations, dampings=dampings, max_missing=max_missing)
            best = select_best_cluster(runs)
            return  best 
            #return gr
        
        #chromosomes = self.values_matrix.values_matrix.chromosomes
        df = None
        chromosomes = self.options.chromosomes["chr"]
        ret_all = list()
        opt_chromosome = self.options.chromosome
        for chr in chromosomes:
            if opt_chromosome is not None and chr != opt_chromosome:
                continue
            ret = list()
            self.options.log(f"Affi for {chr}")
            gc.collect()
            affy_path = self.path_affinity(chr=f"{chr}")
            
            for best, chromosome, start, end in self.map_window_iterator_tabix(function= run_single_run, chromosome=chr):
                if best is not None:
                    best.chromosome = chromosome
                    best.start = start 
                    best.end = end
                    ret.append(best.as_df()) 
                else:
                    self.options.log(f"[run_affinity_propagation] failed to run {chromosome}:{start}-{end}")
            if(len(ret)) > 0:
                df = pd.concat(ret)
                df.to_csv(affy_path, sep="\t",index=False)
                ret_all.append(df)
            else:
                self.options.log(f"Unable to run affinity prpagation for {chr}")
        all = None
        if len(ret_all) > 0:
            all = pd.concat(ret_all)
            # affy_path = self.path_affinity(chr=f"all")
            # all.to_csv()
        return all
    

    
    