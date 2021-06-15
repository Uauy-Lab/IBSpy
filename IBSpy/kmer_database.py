from abc import ABC, abstractmethod, abstractproperty
from pyfaidx import Fasta
from functools import reduce 
from multiprocessing import Pool, TimeoutError, freeze_support

class WindowStats:
    __dict__ = ['seqname', 'start', 'end','total_kmers', 
    'observed_kmers', 'variations', 'kmer_distance']

    def header(self, sep="\t"):
        return sep.join(vars(self))

    def __init__(self):
        self.observed_kmers = 0
        self.variations     = 0
        self.kmer_distance  = 0

    def csv(self, sep="\t"):
        vs = map(lambda w: getattr(self, w), vars(self))
        return sep.join(map(str, vs))

    def add_variation(self, gap_size, kmer_size):
        if gap_size == 0:
            return
        self.variations += 1
        kmer_distance = gap_size - (kmer_size - 1)
        if kmer_distance <= 0:
            kmer_distance = abs(kmer_distance + 1)
        self.kmer_distance += kmer_distance

class KmerDB(ABC):
    
    @abstractproperty
    def builder(self):
        return None

    @abstractmethod
    def __contains__(self, key):
        raise NotImplementedError

    @abstractmethod
    def __len__(self):
        raise NotImplementedError

    def count_kmers_from_sequence(self, kmers): 
        f = filter(lambda k: k in self, kmers)
        return reduce(lambda x, y: x + 1 , f, 0)

    def kmers_stats_from_sequence(self, kmers): 
        stats = WindowStats()
        stats.total_kmers = len(kmers)
        last_present = True
        gap_size = 0
        i = 0
        for k in kmers: 
            i += 1
            present = k in self 
            if not present:
                gap_size += 1
            else:
                stats.observed_kmers += 1
                stats.add_variation(gap_size, self.kmer_size)
                gap_size = 0;
            last_present = present
        
        stats.add_variation(gap_size, self.kmer_size)
        return stats

    def kmers_in_windows(self, path, window_size=1000, processes=1, chunksize=1000):
        
        fasta_iter = FastaChunkReader(path, chunk_size = window_size, 
            kmer_size=self.kmer_size)
        fasta_iter_wrapper = FastaChunkReaderWrapper(fasta_iter, self)
        # with Pool(processes=processes) as pool:
        #     print("We are inside the context of the pool")
        #     ret = pool.map(window_summary, fasta_iter_wrapper,chunksize=chunksize)
        # return ret
        return map(window_summary, fasta_iter_wrapper)

def window_summary(it):
    # print(it)
    seq = it[0]
    kdb = it[1]
    # print(f'Running: {seq}')
    kmers = kdb.builder.sequence_to_kmers(seq['seq'], convert=True)
    stats = kdb.kmers_stats_from_sequence(kmers)
    stats.seqname = seq['seqname']
    stats.start   = seq['start']
    stats.end     = seq['end']
    return stats

class KmerBuilder(ABC):
    def __init__(self, kmer_size):
        self.kmer_size = kmer_size

    @abstractmethod
    def string_to_kmer(self, sequence):
        raise NotImplementedError

    @abstractmethod
    def kmer_to_string(self, sequence):
        raise NotImplementedError

    def compare(self, a, b):
        ka = self.string_to_kmer(a)
        kb = self.string_to_kmer(b)
        return self._builder.compare(ka, kb)

    def sequence_to_kmers(self, sequence, filter_ambiguity=True, convert=False):
        total_kmers = len(sequence) - self.kmer_size + 1
        ret = [None] * (total_kmers)
        sequence = sequence.upper()
        for start in range(0, total_kmers):
            end = start + self.kmer_size
            kmer = sequence[start:end]
            count = kmer.count("A")
            count += kmer.count("T")
            count += kmer.count("G")
            count += kmer.count("C")
            if filter_ambiguity and count != self.kmer_size:
                continue
            ret[start] = self.string_to_kmer(kmer) if convert else kmer 
        return [i for i in ret if i]


class FastaChunkReader:
    def __init__(self, filename, chunk_size=10000, kmer_size=31):
        self.fasta         = Fasta(filename)
        self.current_ref   = 0
        self.current_start = 0
        self.chunk_size    = chunk_size
        self.kmer_size     = kmer_size
        self.seqnames      = list(self.fasta.keys())
        #self.chunk         = chunk
        #self.total_chunks  = total_chunks

    def __iter__(self):
        return self

    def __next__(self):
        if len(self.seqnames) == self.current_ref:
            self.fasta.close()
            raise StopIteration
        seqname = self.seqnames[self.current_ref]
        start = self.current_start
        end = start + self.chunk_size
        self.current_start = end - self.kmer_size
        if end >= len(self.fasta[seqname]):
            self.current_start = 0
            self.current_ref += 1
            end = len(self.fasta[seqname])
        return {
            "seqname": seqname,
            "start": start,
            "end": end,
            "seq": self.fasta[seqname][start:end].seq,
        }

class FastaChunkReaderWrapper:
    def __init__(self, fasta_chunnk_reader, kmer_builder) -> None:
        self.fasta_chunnk_reader = fasta_chunnk_reader
        self.kmer_builder = kmer_builder

    def __iter__(self):
        return self
    
    def __next__(self):
        return [self.fasta_chunnk_reader.__next__(), self.kmer_builder]
