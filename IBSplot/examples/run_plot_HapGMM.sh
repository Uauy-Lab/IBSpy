#!/bin/bash

REF=Jagger
QUERY=Flame
WINDOWS=500000
FILTER=1800
CHROM=chr1A

scrip_dir='path/to/scripts'
in_dir='path/to/database'
out_dir='path/to/output_folder'
chr_dir='path/to/genome_sizes_folder'

python $scrip_dir/plot_HapGMM.v2.py \
	-ho $in_dir/${CHROM}_${REF}_vs_${QUERY}_${WINDOWS}_wd_${FILTER}_flt.tsv.gz \
	-cl $chr_dir/chr_sizes_jagger.genome.txt \
	-rf ${REF} \
	-qr ${QUERY} \
	-chr ${CHROM} \
	-w ${WINDOWS} \
	-vf ${FILTER} \
	-o $out_dir/hap_blocks_${CHROM}_${REF}_vs_${QUERY}_${WINDOWS}_wd_${FILTER}_flt.jpg