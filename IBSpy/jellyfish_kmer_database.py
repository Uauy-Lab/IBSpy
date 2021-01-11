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
            raise "Jellyfish not loaded"
        self.kmer_size = kmer_size
        #self.db  = kmerGWAS.KmerGWAS_database(kmer_size) 
        self.db = None
        self.builder = JellyfishBuilder(kmer_size)
 
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
        raise NotImplementedError

    def __contains__(self, key):
        key.canonicalize()
        return self.db[key] > 0

    def load(self, filename):
        self.db = jf.QueryMerFile(filename)


class JellyfishBuilder(KmerBuilder):
    def __init__(self, kmer_size):
        if not loaded:
            raise "Jellyfish not loaded"
        self.kmer_size = kmer_size

    def string_to_kmer(self, sequence):
        binary_kmer = jf.MerDNA(sequence)
        return binary_kmer

    def kmer_to_string(self, binary_kmer):
        return str(binary_kmer)

    def compare(self, a, b):
        ka = self.string_to_kmer(a)
        kb = self.string_to_kmer(b)
        return self._builder.compare(ka, kb)
