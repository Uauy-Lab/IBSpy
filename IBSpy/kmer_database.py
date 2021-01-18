from abc import ABC, abstractmethod, abstractproperty
from pyfaidx import Fasta
from functools import reduce 

class WindowStats:
    __slots__ = ['seqname', 'start', 'end','total_kmers', 'observed_kmers', 'variations', 'bubble_distance']

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
        stats.observed_kmers = 0
        stats.variations     = 0
        stats.bubble_distance  = 0
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
                if gap_size > 0:
                    stats.variations += 1
                    bubble_distance = gap_size - (self.kmer_size -1)
                    #print(str(i) + ": " + str(gap_size) + ": "+ str(bubble_distance))
                    if bubble_distance <= 0:
                        bubble_distance = abs(bubble_distance + 1)
                    #print(str(i) + ": " + str(gap_size) + ": "+ str(bubble_distance))
                    stats.bubble_distance += bubble_distance
                gap_size = 0;

            last_present = present
        #print(".....")
        if gap_size > 0:
            stats.variations += 1
            bubble_distance = gap_size - (self.kmer_size -1)
            if bubble_distance <= 0:
                bubble_distance = abs(bubble_distance + 1)
            stats.bubble_distance += bubble_distance
        return stats

    def kmers_in_windows(self, path, window_size=1000, chunk=0, total_chunks=1):
        def window_summary(seq):
            kmers = self.builder.sequence_to_kmers(seq['seq'], convert=True)
            stats = self.kmers_stats_from_sequence(kmers)
            stats.seqname = seq['seqname']
            stats.start   = seq['start']
            stats.end     = seq['end']
            return stats
        fasta_iter = FastaChunkReader(path, chunk_size = window_size, 
            kmer_size=self.kmer_size, chunk=chunk, total_chunks=total_chunks)
        return map(window_summary, fasta_iter )

class KmerBuilder(ABC):
    def __init__(self, kmer_size):
        self.kmer_size = kmer_size

    @abstractmethod
    def string_to_kmer(self, sequence):
        raise NotImplementedError

    @abstractmethod
    def kmer_to_string(self, sequence):
        raise NotImplementedError

    def sequence_to_kmers(self, sequence, filter_ambiguity=True, convert=False):
        ret = [None] * (len(sequence) - self.kmer_size + 1)
        sequence = sequence.upper()
        for start in range(0, len(sequence) - self.kmer_size + 1):
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
    def __init__(self, filename, chunk_size=10000, kmer_size=31, 
        chunk=0, total_chunks=1):
        self.fasta         = Fasta(filename)
        self.current_ref   = 0
        self.current_start = 0
        self.chunk_size    = chunk_size
        self.kmer_size     = kmer_size
        self.seqnames      = list(self.fasta.keys())
        self.chunk         = chunk
        self.total_chunks  = total_chunks

        

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
