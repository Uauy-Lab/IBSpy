from .kmer_database import KmerDB, KmerBuilder, FastaChunkReader

# from ctypes import *
import os
import inspect
loaded=False
try:
    import dna_jellyfish as jf
    loaded=True
except ImportError:
    pass


class JellyfishSDB(KmerDB):
    def __init__(self, kmer_size):
        if not loaded:
            raise ImportError("Jellyfish not loaded")
        self.kmer_size = kmer_size
        self.db = None
        self.builder = JellyfishBuilder(kmer_size)
 
    def builder(self):
        return self.builder

    def __getitem__(self, index):
        return self.db[index]

    def __len__(self):
        raise NotImplementedError

    def __contains__(self, key):
        key.canonicalize()
        return self.db[key] > 0

    def load(self, filename):
        self.db = jf.QueryMerFile(filename)


class JellyfishBuilder(KmerBuilder):
    def __init__(self, kmer_size):
        if not loaded:
            raise ImportError("Jellyfish not loaded")
        self.kmer_size = kmer_size

    def string_to_kmer(self, sequence):
        binary_kmer = jf.MerDNA(sequence)
        return binary_kmer

    def kmer_to_string(self, binary_kmer):
        return str(binary_kmer)

