from .kmer_database import KmerDB, KmerBuilder, FastaChunkReader
from .kmerGWAS_kmer import KmerGWASDB, KmerGWASDBBuilder
from .jellyfish_kmer_database import JellyfishSDB, JellyfishBuilder
from .IBSpy_window_count import window_count, parse_arguments
from .IBSpy_results import IBSpyResults
import argparse

def main(): 
	window_count(parse_arguments())

if __name__ == '__main__':
	main()

