from .kmer_database import KmerDB, KmerBuilder, FastaChunkReader

# from ctypes import *
import os
import inspect
import kmerGWAS
import cython

# from Bio import SeqIO
# import pysam
from pyfaidx import Fasta


class KmerGWASDB(KmerDB):
    def __init__(self, kmer_size):
        self._builder = KmerGWASDBBuilder(kmer_size)
        self.kmer_size = kmer_size
        self.db  = kmerGWAS.KmerGWAS_database(kmer_size) 
        self.builder = KmerGWASDBBuilder(kmer_size)

    def builder(self):
        return self.builder

    def load_from_fasta(self, filename, buffer_size=10000):
        # fasta = pysam.FastaFile(filename)
        fasta_iter = FastaChunkReader(
            filename, chunk_size=buffer_size, kmer_size=self.kmer_size
        )        
        for chunk in fasta_iter:
            self.db.add_kmers(chunk['seq'])
        self.db.sort_unique()

    def __getitem__(self, index):
        return self.db[index]

    def __len__(self):
        return len(self.db)

    def __contains__(self, key):
        return key in self.db

    def save(self, filename):
        self.db.save(filename)

    def load(self, filename):
    	self.db.read_mmap(filename)


class KmerGWASDBBuilder(KmerBuilder):
    # cdef char* _builder
    def __init__(self, kmer_size):
        self.kmer_size = kmer_size
        self._builder = kmerGWAS.KmerGWAS_builder(kmer_size)

    def string_to_kmer(self, sequence):
        binary_kmer = self._builder.string_to_kmer(sequence)
        return binary_kmer

    def kmer_to_string(self, binary_kmer):
        return self._builder.kmer_to_string(binary_kmer).decode("UTF-8")

    def compare(self, a, b):
        ka = self.string_to_kmer(a)
        kb = self.string_to_kmer(b)
        return self._builder.compare(ka, kb)
