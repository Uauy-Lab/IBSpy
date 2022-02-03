import string
import pandas as pd

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

    def _build_dataset(self):
        dfs = []
        for ref in self.references:
            df = self.values_matrix_for_reference(ref)
            dfs.append(df)
        return pd.concat(dfs, join="inner")
    
    @property
    def values_matrix(self):
        if self._values_matrix is not None:
            return self._values_matrix
        self._values_matrix = self._build_dataset()
        return self._values_matrix
    
    def seqname_iterator(self):
        mat = self.values_matrix
        seqnames = mat['seqname'].unique()
        for seq in seqnames:
            yield mat[mat['seqname'] == seq]

    def matrix_iterator(self):
        for seq in self.seqname_iterator():
            for x in range(0, len(seq) ,self.options.affinity_blocks):
                yield seq.iloc[x:x+self.options.affinity_blocks]

    def to_csv(self,  *args, **kwargs):
        self.values_matrix.to_csv( *args, **kwargs)

    