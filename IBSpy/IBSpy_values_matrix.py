import os
import string
import sys
from numpy import int64
import pandas as pd
import psutil
from pyranges import PyRanges
from multiprocess import Pool

 
from .IBSpy_options import IBSpyOptions
from .IBSpy_results import IBSpyResults

class IBSpyValuesMatrix:
    def __init__(self, result_set) -> None:
        self.options: IBSpyOptions = result_set.options
        self.result_set = result_set
        self.samples_df = pd.read_csv(self.options.metadata, delimiter='\t')
        self._samples = None
        self._references = None
        self.references = self.samples_df["reference"].unique()
        self.samples = self.samples_df["query"].unique()
        self._values_matrix  = None
        self._chromosome_lengths = None


    @property
    def references(self):
        return self._references
    @references.setter
    def references(self, value):
        if self.options.references is not None:
            value = filter(lambda r: r in self.options.references, value)
        self._references = value
        
    @property
    def samples(self):
        return self._samples
    @samples.setter
    def samples(self, value):
        if self.options.samples is not None:
            value = filter(lambda r: r in self.options.samples, value)
        self._samples = value

    def _summarize_single_line(self, i:int, samples: pd.DataFrame):
        row = samples.iloc[i]
        reference= row['reference']
        sample_name = row['query']
        score = self.options.score
        path_ref=f'{self.options.folder_for_reference(reference=reference)}'
        out_path=f"{path_ref}/{sample_name}.{score}.pkl.gz"
        if os.path.exists(out_path):
            self.options.log(f" File already exits: {out_path}")
            return i
        path=row['path'] + "/" + row['file']
        self.options.log(f" Reading {reference} {sample_name}: {path}")
        ibr = IBSpyResults(path, self.options)
        df: pd.DataFrame = ibr.count_by_windows()
        df.rename(columns = {'seqname':'Chromosome', 'chromosome':'Chromosome', 'start':'Start', 'end':"End", 'orientation':"Strand"}, inplace = True)
        df.to_pickle(out_path, compression={'method': 'gzip', 'compresslevel': 9} )
        return i        

    def _values_matrix_for_reference(self, reference:string):
        score = self.options.score
        samples = self.samples_df[self.samples_df['reference'] == reference]  
        path_ref=f'{self.options.folder_for_reference(reference=reference)}'      
        path_matrix=f'{path_ref}/{self.options.file_prefix}.merged.pickle.gz'
        self.options.log(f"Searching {path_matrix}")
        if os.path.isfile(path_matrix):
            return pd.read_pickle(path_matrix)
        nrows = samples.shape[0]
        with Pool(self.options.pool_size) as p:
            wrapped_function = lambda x: self._summarize_single_line(x, samples)
            p.map(wrapped_function,range(0, nrows),self.options.chunks_in_pool )

        values_matrix = pd.DataFrame()
        prep = {}
        for index, row in samples.iterrows():
            sample_name = row['query']
            self.options.log(f" Reading {path_ref}/{sample_name}.{score}.pkl.gz")
            df = pd.read_pickle(f"{path_ref}/{sample_name}.{score}.pkl.gz")
            prep[sample_name] = df[self.options.stat] 
        values_matrix = pd.DataFrame(prep)
        names = df[['Chromosome','Start','End']]
        values_matrix = pd.concat([names,values_matrix], axis=1)
        values_matrix.to_pickle(path_matrix, compression={'method': 'gzip', 'compresslevel': 9} )
        return values_matrix

    def _build_dataset(self) -> pd.DataFrame:
        dfs = []
        for ref in self.references:
            df = self._values_matrix_for_reference(ref)
            dfs.append(df)
        mem_usage = psutil.Process().memory_info().rss / (1024 * 1024 * 1024)
        self.options.log(f"Merging single values for matrix {mem_usage} Gb")

        return self.rename_sequnames( pd.concat(dfs, join="inner") )

    def rename_sequnames(self, df) -> pd.DataFrame:
        if self.options.chromosome_mapping is None:
            return df
        mapping = self.options.mapping_seqnames
        df['Chromosome'] = df['Chromosome'].apply(lambda x: mapping.get(x, x)) 
        return df

    @property
    def values_matrix(self) -> PyRanges:
        if self._values_matrix is not None:
            return self._values_matrix
        if os.path.isfile(self.path):
            self.options.log(f"Loading {self.path}")
            self._values_matrix = pd.read_csv(self.path, sep="\t")
        else:
            self.options.log("Building matrix")
            self._values_matrix = self._build_dataset()
        self.options.log("Converting to PyRanges")
        self._values_matrix = PyRanges(self._values_matrix, int64=True)
        return self._values_matrix

    @property
    def chromosome_lengths(self):
        if self._chromosome_lengths is not None:
            return self._chromosome_lengths
        chromosomes = self.values_matrix.chromosomes
        self._chromosome_lengths = {}
        for chr in chromosomes:
            self._chromosome_lengths[chr] = max(self.values_matrix[chr].End )
        return self._chromosome_lengths
    
    def seqname_iterator(self):
        mat = self.values_matrix
        for seq in mat.chromosomes:
            yield mat[seq], seq

    def matrix_iterator(self, window_size=1000000):
        chr_lens = self.chromosome_lengths
        for seq, chromosome in self.seqname_iterator():
            for x in range(0, chr_lens[chromosome] , window_size):
                yield seq[x:x+window_size]

    def to_csv(self,  *args, **kwargs):
        self.values_matrix.to_csv( *args, **kwargs)

    def save(self):
        self.values_matrix.to_csv(self.path,sep="\t")

    @property
    def path(self):
        prefix = self.options.output_folder
        file = self.options.file_prefix
        return f"{prefix}/{file}.csv.gz"
    
    @property
    def cache_path(self):
        prefix = self.options.cache_folder
        file = self.options.file_prefix
        return f"{prefix}/{file}.pickle.gz"

    