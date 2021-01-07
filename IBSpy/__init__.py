from .kmer_database import KmerDB, KmerBuilder, FastaChunkReader
from .kmerGWAS_kmer import KmerGWASDB, KmerGWASDBBuilder
from .bin.IBSpy_window_count import window_count
import argparse


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-w", "--window_size", default=500000,
		help="window size to analyze", type=int)
	parser.add_argument("-k", "--kmer_size", default=31,
		help="Kmer size of the database", type=int)
	parser.add_argument("-d", "--database", default="./tests/data/test4B.jagger.kmerGWAS_k31",
		help="Kmer database")
	parser.add_argument("-r", "--reference", default="./tests/data/test4B.stanley.fa",
		help="The reference with the position of the kmers")
	parser.add_argument("-z", "--compress", default=False, action="store_true", 
		help="When an ouput file is present, it is compressed as .gz")
	parser.add_argument("-o", "--output", default=None, 
		help="Output file. If missing, the ouptut is sent to stdout")
	args = parser.parse_args()
	window_count(args)

if __name__ == '__main__':
	main()
