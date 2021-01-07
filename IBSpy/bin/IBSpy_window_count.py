#!/usr/bin/env python

import argparse
import gzip
from sys import stdout
from IBSpy import KmerGWASDBBuilder, KmerGWASDB, FastaChunkReader 

def window_count(args):
	window_size = 1000
	kmers_path  = "./tests/data/test4B.jagger.kmerGWAS_k31"
	reference  = "./tests/data/test4B.stanley.fa"
	kmer_size   = 31

	kmerdb = KmerGWASDB(args.kmer_size)
	kmerdb.load(args.database)
	windows = kmerdb.kmers_in_windows(args.reference, window_size=args.window_size)
	printed = False
	out = stdout

	if(args.compress and args.output is not None) :
		out = gzip.open(args.output + ".gz", 'w') 

	if(not args.compress and args.output is not None) :
		out = open(args.output , 'w') 

	for w in windows:
		if not printed:
			line = "\t".join(w.keys())
			line += "\n"
			if args.compress:
				line=line.encode()
			out.write(line)
			printed = True
		line = "\t".join(map(str,w.values()))
		line += "\n"  
		if args.compress:
			line=line.encode()
		out.write(line)

	if(args.output is not None):
		out.close()


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

	# args = {
	# window_size : 1000000,
	# kmers_path  : "./tests/data/test4B.jagger.kmerGWAS_k31",
	# reference   : "./tests/data/test4B.stanley.fa",
	# kmer_size   : 31,
	# compress    : True,
	# output      : None
	# }

	#We need to add the different databases

	

	
	#parser = argparse.ArgumentParser()
	#parser.parse_args()

	