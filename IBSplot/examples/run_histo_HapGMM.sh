#!/bin/bash


REF=Jagger
QUERY=Flame
WINDOWS=500000
FILTER=1800
CHROM=chr1A

scrip_dir='path/to/scripts'
in_dir='path/to/database'
out_dir='path/to/output_folder'

python $scrip_dir/histo_GMM.py \
	-db $in_dir/kmeribs-Wheat_${REF}-${QUERY}_${WINDOWS}.tsv.gz \
	-gm $in_dir/gmm_${REF}_vs_${QUERY}_${WINDOWS}_wd_${FILTER}_flt.tsv.gz \
	-rf ${REF} \
	-qr ${QUERY} \
	-w ${WINDOWS} \
	-vf ${FILTER} \
	-do $out_dir/hist_${REF}_vs_${QUERY}_${WINDOWS}_wd_${FILTER}_flt.pdf \
	-go $out_dir/hist_${REF}_vs_${QUERY}_${WINDOWS}_wd_whole_genome.pdf