import pandas as pd
import pyranges as pr
from pyranges import PyRanges
import warnings 

class BlockMapping:
	
	def __init__(self, path, int64=True) -> None:
		mapping = pd.read_csv(path, delimiter='\t')
		col_rename = {
			'chromosome':'Chromosome',
			'seqnames':'Chromosome', 
			'seqname':'Chromosome', 
			'start':'Start', 
			'end':"End", 
			'orientation':"Strand"}
		mapping.rename(columns = col_rename, inplace = True)
		self.int64 = int64
		self.mapping = PyRanges(mapping, int64=self.int64)	
		#TODO: remove once pyranges is updated to support numpy > 1.20
		warnings.filterwarnings("ignore") #,  module="numpy" , category=DeprecationWarning, 
		


	def mapping_for_range(self, chromosome, start, end, assembly = None):
		gr = self.mapping[chromosome, start:end]
		if assembly is not None:
			gr = gr[gr.assembly==assembly]
		blocks_nos = gr.block_no.unique()
		return self.mapping[self.mapping.block_no.isin(blocks_nos)]

	def all_regions_for(self, chromosome, start, end, assembly = None):
		region = pr.PyRanges(chromosomes=chromosome, starts=[start], ends=[end], strands=["+"], int64=self.int64)
		mapping = self.mapping_for_range(chromosome, start, end, assembly=assembly)
		mapping = mapping[mapping.Chromosome != chromosome ]
		mapping = pr.concat([region, mapping])
		merged = mapping.merge()
		return merged

