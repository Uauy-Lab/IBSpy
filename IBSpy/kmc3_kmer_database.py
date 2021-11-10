from .kmer_database import KmerDB, KmerBuilder, FastaChunkReader

import os
import inspect
import sys

loaded=False
try:
    import py_kmc_api as pka
    loaded=True
except ImportError:
    pass

class KMC3DB(KmerDB):
    def __init__(self, kmer_size):
        if not loaded:
            raise ImportError("KMC3 not loaded")
        self.kmer_size = kmer_size
        self.db = None
        self.db_info = None
        self.kmer_data_base = pka.KMCFile()
        self.builder = KMCBuilder(kmer_size )

    def builder(self):
        return self.builder


    def load(self, filename):
        self.db = pka.KMCFile()
        opened = self.db.OpenForRA(filename)
        if not opened:
            raise Exception(f'Unable to open KMC3 database {filename}')
        self.db_info = self.db.Info()

    # def canonicalize(kmer):
    #     tmp = pka.CKmerAPI(kmer)
    #     tmp 
    # Finish this method if the unit test fails with 
    # kmers not found. 

    def __getitem__(self, kmer):
        if not self.db.IsKmer(kmer):
            kmer = None
        return kmer

    def __len__(self):
        return self.db_info.total_kmers

    def __contains__(self, kmer):
        ret = self.db.IsKmer(kmer)
        kmer.reverse()
        ret |= self.db.IsKmer(kmer)
        return ret


class KMCBuilder(KmerBuilder):

    def __init__(self, kmer_size):
        if not loaded:
            raise ImportError("KMC3 not loaded")
        self.kmer_size = kmer_size
        

    def string_to_kmer(self, sequence):
        kmer = pka.KmerAPI(self.kmer_size)
        converted = kmer.from_string(sequence)
        return kmer

    def kmer_to_string(self, binary_kmer):
        return str(binary_kmer)
