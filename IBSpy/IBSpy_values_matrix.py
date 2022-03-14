import os
import string
import sys
import threading
from typing import Dict
from numpy import int64
import pandas as pd
import psutil
import pysam
from pyranges import PyRanges
from multiprocess import Pool

 
from .IBSpy_options import IBSpyOptions
from .IBSpy_results import IBSpyResults

class IBSpyValuesMatrix:
    def __init__(self, result_set) -> None:
        self.options: IBSpyOptions = result_set.options
        self._locks : Dict[str, threading.Lock]  = {}
        self.result_set = result_set
        self.samples_df = pd.read_csv(self.options.metadata, delimiter='\t')
        self._samples = None
        self._references = None
        self.references = self.samples_df["reference"].unique()
        self.samples = self.samples_df["query"].unique()
        self._values_matrix  = None
        self._chromosome_lengths = None
        self._tabix = None 
        


    @property
    def references(self):
        return self._references
    @references.setter
    def references(self, value):
        if self.options.references is not None:
            value = filter(lambda r: r in self.options.references, value)
        self._references = value
        for ref in self._references:
            self._locks[ref] = threading.Lock()
        
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

    def _run_summirize_lines(self, reference:string) -> None:
        samples = self.samples_df[self.samples_df['reference'] == reference]  
        nrows = samples.shape[0]
        with Pool(self.options.pool_size) as p:
            wrapped_function = lambda x: self._summarize_single_line(x, samples)
            p.map(wrapped_function,range(0, nrows),self.options.chunks_in_pool )

    def _merge_values_matrix(self, samples:pd.DataFrame , path_ref: string) -> pd.DataFrame:
        score = self.options.score
        prep = {}
        for index, row in samples.iterrows():
            sample_name = row['query']
            self.options.log(f" Reading {path_ref}/{sample_name}.{score}.pkl.gz")
            df = pd.read_pickle(f"{path_ref}/{sample_name}.{score}.pkl.gz")
            prep[sample_name] = df[self.options.stat] 
        values_matrix = pd.DataFrame(prep)
        names = df[['Chromosome','Start','End']]
        values_matrix = pd.concat([names,values_matrix], axis=1)
        return values_matrix

    def _values_matrix_for_reference(self, reference:string) -> pd.DataFrame:
        samples = self.samples_df[self.samples_df['reference'] == reference]  
        path_ref=f'{self.options.folder_for_reference(reference=reference)}'      
        path_matrix=f'{path_ref}/{self.options.file_prefix}.merged.pickle.gz'
        self.options.log(f"Searching {path_matrix}")
        if os.path.isfile(path_matrix):
            return pd.read_pickle(path_matrix)
        self._run_summirize_lines(reference=reference)
        values_matrix = self._merge_values_matrix(samples, path_ref)
        values_matrix.to_pickle(path_matrix, compression={'method': 'gzip', 'compresslevel': 9} )
        return values_matrix

    def _read_single_line(self, i:int, samples: pd.DataFrame):
        row = samples.iloc[i]
        reference= row['reference']
        sample_name = row['query']
        score = self.options.score
        path_ref=f'{self.options.folder_for_reference(reference=reference)}'
        out_path=f"{path_ref}/{sample_name}.{score}.pkl.gz"
        df = pd.read_pickle(f"{path_ref}/{sample_name}.{score}.pkl.gz")
        df["sample"] = sample_name
        return df   

    def _merge_values_long(self, samples: pd.DataFrame, reference:string) -> pysam.TabixFile:
        path_ref   = f'{self.options.folder_for_reference(reference=reference)}' 
        path_df    = f'{path_ref}/{self.options.file_prefix}.merged.tsv'
        path_tabix = f'{path_df}.gz'
        if os.path.isfile(path_tabix):
            return pysam.TabixFile(filename=path_tabix,index=f"{path_tabix}.csi" )
        self._run_summirize_lines(reference=reference)
        samples = self.samples_df[self.samples_df['reference'] == reference]  
        nrows = samples.shape[0]
        with Pool(self.options.pool_size) as p:
            wrapped_function = lambda x: self._read_single_line(x, samples)
            dfs = p.map(wrapped_function,range(0, nrows),self.options.chunks_in_pool )
        
        df = pd.concat(dfs)
        df = self.rename_seqnames(df)
        df.sort_values(by=["Chromosome", "Start","End"], inplace=True)
        df.to_csv(path_df,sep="\t",index=False )
        pysam.tabix_index(path_df,seq_col=0, start_col=1, end_col=2, line_skip=1,csi=True)
        return pysam.TabixFile(filename=path_tabix, index=f"{path_tabix}.csi" )

    @property
    def merged_values(self) -> dict[str, pysam.TabixFile]:
        references = self.references
        # if self._tabix  is not None:
        #     return self._tabix
        ret = {}
        for ref in references:
            samples = self.samples_df[self.samples_df['reference'] == ref]
            ret[ref] = self._merge_values_long(samples=samples, reference=ref)
        # self._tabix = ret
        return ret

    def acquire(self, assembly):
        self._locks[assembly].acquire()

    def release(self, assembly):
        self._locks[assembly].release()     

    def values_for_region(self, chromosome, start, end):
        colnames = ["Chromosome",	"Start",	"End",	self.options.score,	"mean",	"median",	"variance",	"std",	"sample"]
        for assembly in self.merged_values.keys():
            tabix = self.merged_values[assembly]
            if chromosome not in tabix.contigs:
                continue
            # print(f'[values_for_region]{chromosome}:{start}-{end} is in {assembly}')
           # print(tabix.contigs)
            self.acquire(assembly)
            regions = pd.DataFrame(tabix.fetch(reference=chromosome, start=start, end=end, parser=pysam.asTuple()), columns=colnames)
            self.release(assembly)
            if regions is not None:
                regions[['Start']] = regions[['Start']].apply(pd.to_numeric) 
                regions[['End']] = regions[['End']].apply(pd.to_numeric) 
                regions[[self.options.score]] = regions[[self.options.score]].apply(pd.to_numeric) 
                regions[['mean']] = regions[['mean']].apply(pd.to_numeric) 
                regions[['median']] = regions[['median']].apply(pd.to_numeric) 
                regions[['variance']] = regions[['variance']].apply(pd.to_numeric) 
                regions[['std']] = regions[['std']].apply(pd.to_numeric) 
                # regions[['sample']] = regions[['sample']].apply(str) 
                #regions[['Chromosome']] = regions[['Chromosome']].apply(str) 
                regions=regions.convert_dtypes()
            # print(regions.dtypes)
            # print("[values_for_region]")
            # print(regions)
            # print("_____________")
            return regions
        raise f"Chromsome not found: {chromosome}"

    def _build_dataset(self) -> pd.DataFrame:
        dfs = []
        for ref in self.references:
            df = self._values_matrix_for_reference(ref)
            dfs.append(df)
        mem_usage = psutil.Process().memory_info().rss / (1024 * 1024 * 1024)
        self.options.log(f"Merging single values for matrix {mem_usage} Gb")

        return self.rename_seqnames( pd.concat(dfs, join="inner") )

    def rename_seqnames(self, df) -> pd.DataFrame:
        if self.options.chromosome_mapping is None:
            return df
        mapping = self.options.mapping_seqnames
        # print(df)
        df['Chromosome'] = df['Chromosome'].apply(lambda x: mapping.get(x, x)) 
        # print("...")
        # print(df)
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

    