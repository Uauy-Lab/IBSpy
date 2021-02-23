#!/bin/bash

scrip_dir='path/to/scripts'
gmm_dir='path/to/gmm_file_folder'
db_dir='path/to/database'
chr_dir='path/to/chromosome_lengths_folder'
out_dir='path/to/output_folder'


WINDOWS=100000
REF=symattis
QUERY=flame
FILTER=400
CHROM=chr6A

python $scrip_dir/plot_HapGMM.py \
	-gm $gmm_dir/gmm_${REF}_vs_${QUERY}_${WINDOWS}_wd_${FILTER}_flt.tsv \
	-db $db_dir/${QUERY}_${REF}_${WINDOWS}.tsv.gz  \
	-cl $chr_dir/chr_sizes_${REF}.genome.txt \
	-rf ${REF} \
	-qr ${QUERY} \
	-chr ${CHROM} \
	-w ${WINDOWS} \
	-vf ${FILTER} \
	-st 3 \
	-o $out_dir/hap_blocks_${CHROM}_${REF}_vs_${QUERY}_${WINDOWS}_wd_${FILTER}_flt.jpg