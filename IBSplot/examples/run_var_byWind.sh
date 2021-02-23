#!/bin/bash

scrip_dir='path/to/scripts'
db_dir='path/to/database'
chr_dir='path/to/chromosome_lengths_folder'
out_dir='path/to/output_folder'

IN_WINDOWS=100000
OUT_WINDOWS=500000
REF=symattis
QUERY=flame

python $scrip_dir/var_byWind.py \
	-db $db_dir/${QUERY}_${REF}_${IN_WINDOWS}.tsv.gz \
	-cl $chr_dir/chr_sizes_${REF}.genome.txt \
	-w ${OUT_WINDOWS} \
	-o $out_dir/${QUERY}_${REF}_${OUT_WINDOWS}.tsv