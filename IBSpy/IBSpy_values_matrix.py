import string
import pandas as pd
from pyranges import PyRanges
 
from .IBSpy_options import IBSpyOptions
from .IBSpy_results import IBSpyResults

class IBSpyValuesMatrix:
    def __init__(self, result_set) -> None:
        self.options: IBSpyOptions = result_set.options
        self.result_set = result_set
        self.samples_df = pd.read_csv(self.options.metadata, delimiter='\t')
        self.references = self.samples_df["reference"].unique()
        self.samples = self.samples_df["query"].unique()
        self._values_matrix  = None
        self._chromosome_lengths = None

    def _values_matrix_for_reference(self, reference:string):
        samples = self.samples_df[self.samples_df['reference'] == reference]
        values_matrix = pd.DataFrame()
        for index, row in samples.iterrows():
            path=row['path'] + "/" + row['file']
            sample_name = row['query']
            ibr = IBSpyResults(path, self.options)
            df: pd.DataFrame = ibr.count_by_windows()
            values_matrix[sample_name] = df[self.options.stat]
        df.rename(columns = {'seqname':'Chromosome', 'chromosome':'Chromosome', 'start':'Start', 'end':"End", 'orientation':"Strand"}, inplace = True)
        names = df[['Chromosome','Start','End']]
        values_matrix = pd.concat([names,values_matrix], axis=1)
        return values_matrix

    def _build_dataset(self):
        dfs = []
        for ref in self.references:
            df = self._values_matrix_for_reference(ref)
            dfs.append(df)
        
        return self.rename_sequnames( pd.concat(dfs, join="inner") )

    def mapping_seqnames(self):
        map_df = pd.read_csv(self.options.chromosome_mapping, delimiter='\t')
        ret = {}
        for index, row in map_df.iterrows():
            ret[row['original']] = row['mapping']
        return ret
    
    def rename_sequnames(self, df):
        if self.options.chromosome_mapping is None:
            return df
        mapping = self.mapping_seqnames()
        df['Chromosome'] = df['Chromosome'].apply(lambda x: mapping[x]) 
        return df

    @property
    def values_matrix(self):
        if self._values_matrix is not None:
            return self._values_matrix
        self._values_matrix = self._build_dataset()
        self._values_matrix = PyRanges(self._values_matrix)
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

    