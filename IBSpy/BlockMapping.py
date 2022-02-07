import pandas as pd
import pyranges as pr
from pyranges import PyRanges
 

class BlockMapping:
	def __init__(self, path) -> None:
		mapping = pd.read_csv(path, delimiter='\t')
		col_rename = {
			'chromosome':'Chromosome',
			'seqnames':'Chromosome', 
			'start':'Start', 'end':"End", 
			'orientation':"Strand"}
		mapping.rename(columns = col_rename, inplace = True)
		self.mapping = PyRanges(mapping)	

