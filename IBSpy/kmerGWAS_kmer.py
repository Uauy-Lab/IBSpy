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

    def load_from_fasta(self, filename, buffer_size=10000):
        # fasta = pysam.FastaFile(filename)
        fasta_iter = FastaChunkReader(
            filename, chunk_size=buffer_size, kmer_size=self.kmer_size
        )
        for chunk in fasta_iter:
            print(chunk)
            # max_width = len(ref)

        # print(references)

    # def load_kmers_from_sequence(self, string):

    def __len__(self):
        raise NotImplementedError

    def __contains__(self, key):
        raise NotImplementedError


class KmerGWASDBBuilder(KmerBuilder):
    # cdef char* _builder
    def __init__(self, kmer_size):
        self.kmer_size = kmer_size
        self._builder = kmerGWAS.KmerGWAS_builder(kmer_size)
        # print(self._builder);

    def string_to_kmer(self, sequence):
        binary_kmer = self._builder.string_to_kmer(sequence)
        return binary_kmer

    def kmer_to_string(self, binary_kmer):
        return self._builder.kmer_to_string(binary_kmer).decode("UTF-8")
