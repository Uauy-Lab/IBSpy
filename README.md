# IBSpy

![Python package](https://github.com/Uauy-Lab/IBSpy/workflows/Python%20package/badge.svg)
[![Maintainability](https://api.codeclimate.com/v1/badges/5a4b1b0e89f7f9f8c34c/maintainability)](https://codeclimate.com/github/Uauy-Lab/IBSpy/maintainability)

Python library to identify Identical By State regions

## Installyng IBSpy

There easiest way to install IBSpy is to use pip3. 

```sh
pip3 install IBSpy
```


If ```pip3``` fails, you can clone the project and compiling it with:

```sh
pip3 install cython biopython pyfaidx
python3 setup.py develop
```

Then you should have the  IBSpy command available. 


## Preparing the databases

IBSpy requires to have a kmer database from the sequencing files. Currently two formats are supported:

  1. Jellyfish: Follow the instructions in its [website](https://github.com/gmarcais/Jellyfish/blob/master/doc/Readme.md)
  2. kmerGWAS: Has an adhoc file format that contains only the kmers in a binary representation, sorted. This option is faster than the jellyfish version, but creating the kmer table is less straight forward. The manual is [here](https://github.com/voichek/kmersGWAS/blob/master/manual.pdf).

## Running IBSPy

IBSpy has relatively few options, you can look at them with the ```---help``` command. 

```sh
IBSPy --help
usage: IBSPy [-h] [-w WINDOW_SIZE] [-k KMER_SIZE] [-d DATABASE] [-r REFERENCE]
             [-z] [-o OUTPUT] [-f {kmerGWAS,jellyfish}]

optional arguments:
  -h, --help            show this help message and exit
  -w WINDOW_SIZE, --window_size WINDOW_SIZE
                        window size to analyze
  -k KMER_SIZE, --kmer_size KMER_SIZE
                        Kmer size of the database
  -d DATABASE, --database DATABASE
                        Kmer database
  -r REFERENCE, --reference REFERENCE
                        The reference with the position of the kmers
  -z, --compress        When an ouput file is present, it is compressed as .gz
  -o OUTPUT, --output OUTPUT
                        Output file. If missing, the ouptut is sent to stdout
  -f {kmerGWAS,jellyfish}, --database_format {kmerGWAS,jellyfish}
                        Database format (kmerGWAS, jellyfish)
```

To generate the table with the number of observed kmers and variants run the following command, using the kmer database from kmerGWAS use the following command:

```sh
 IBSpy --output "kmer_windows_LineXXX.tsv.gz" -z --database kmers_with_strand  --reference arinaLrFor.fa --window_size 50000 --compress --database_format kmerGWAS
```


