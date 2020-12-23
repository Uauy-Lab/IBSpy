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
	def str_to_kmer(sequence):
		raise NotImplementedError

	@abstractmethod
	def kmer_to_str(sequence):
		raise NotImplementedError