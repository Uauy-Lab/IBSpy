import pandas as pd
import pyranges as pr
from pyranges import PyRanges
 

class BlockMapping:
	def __init__(self, path) -> None:
		mapping = pd.read_csv(path, delimiter='\t')
		self.mapping = PyRanges(mapping)	

