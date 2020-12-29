from abc import ABC, abstractmethod

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