from abc import ABC, abstractmethod
from pyfaidx import Fasta

class KmerDB(ABC):
	@abstractmethod
	def __contains__(self, key):
		raise NotImplementedError

	@abstractmethod
	def	__len__(self):
		raise NotImplementedError

class KmerBuilder(ABC):

	def __init__(self, kmer_size):
		self.kmer_size = kmer_size

	@abstractmethod
	def string_to_kmer(self,sequence):
		raise NotImplementedError

	@abstractmethod
	def kmer_to_string(self,sequence):
		raise NotImplementedError

	def sequence_to_kmers(self, sequence, filter_ambiguity = True):
		ret = [None] * (len(sequence) - self.kmer_size + 1)
		sequence = sequence.upper()
		for start in range(0, len(sequence) - self.kmer_size + 1):
			end = start + self.kmer_size
			kmer = sequence[start:end]
			count =  kmer.count('A')
			count += kmer.count('T')
			count += kmer.count('G')
			count += kmer.count('C') 
			if filter_ambiguity and count != self.kmer_size:
				continue 
			ret[start] = kmer
		return  [i for i in ret if i] 


class FastaChunkReader:
	
	def __init__(self, filename, chunk_size=10000 , kmer_size=31):
		self.fasta = Fasta(filename)
		self.current_ref   = 0
		self.current_start = 0
		self.chunk_size    = chunk_size
		self.kmer_size     = kmer_size
		self.seqnames      = list(self.fasta.keys())

	def __iter__(self):
		return self

	def __next__(self):
		if(len(self.seqnames) == self.current_ref):
			raise StopIteration
		#print(self.seqnames)
		#print(self.current_ref)
		seqname = self.seqnames[self.current_ref]
		start   = self.current_start
		end    = start + self.chunk_size
		self.current_start = end - self.kmer_size
		if(end >= len(self.fasta[seqname]) ):
			self.current_start  = 0
			self.current_ref   += 1
			end = len(self.fasta[seqname])
		return {
		'seqname': seqname,
		'start'  : start,
		'end'   : end,
		'seq': self.fasta[seqname][start:end].seq}

