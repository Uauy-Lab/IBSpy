import argparse
import gzip
from sys import stdout
from IBSpy import KmerGWASDBBuilder, KmerGWASDB, FastaChunkReader, JellyfishSDB, JellyfishBuilder 


def print_line(out, line, compress):
	line += "\n"
	if compress:
		line=line.encode()
	out.write(line)

def open_out(args):
	out = stdout
	if(args.compress and args.output is not None) :
		out = gzip.open(args.output + ".gz", 'w') 
	if(not args.compress and args.output is not None) :
		out = open(args.output , 'w')
	return out

def close_out(out, args):
	if(args.output is not None):
		out.close()

def open_db(args):
	kmerdb = None
	if(args.database_format=="kmerGWAS", mmap=False):
		kmerdb = KmerGWASDB(args.kmer_size)
	if(args.database_format=="kmerGWAS_mmap"):
		kmerdb = KmerGWASDB(args.kmer_size, mmap=True)
	if(args.database_format=="jellyfish"):
		kmerdb = JellyfishSDB(args.kmer_size)
	kmerdb.load(args.database)
	return kmerdb


def window_count(args):
	window_size = 1000
	kmers_path  = "./tests/data/test4B.jagger.kmerGWAS_k31"
	reference  = "./tests/data/test4B.stanley.fa"
	kmer_size   = 31
	kmerdb = open_db(args)
	windows = kmerdb.kmers_in_windows(args.reference, window_size=args.window_size)
	printed = False
	out = open_out(args)

	for w in windows:
		if not printed:
			line = w.header()
			print_line(out, line, args.compress)
			printed = True
		line = w.csv()
		print_line(out, line, args.compress)
	close_out(out, args)
	

def parse_arguments():
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
	parser.add_argument("-f", "--database_format", default="kmerGWAS", choices=["kmerGWAS", "kmerGWAS_mmap","jellyfish"],
		help="Database format (kmerGWAS, jellyfish)")
	args = parser.parse_args()
	return args


def main():
	args = parse_arguments()
	window_count(args)

if __name__ == '__main__':
	main()
