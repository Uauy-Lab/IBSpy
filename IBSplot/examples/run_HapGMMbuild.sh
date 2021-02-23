#!/bin/bash

scrip_dir='path/to/scripts'
db_dir='path/to/database'
out_dir='path/to/output_folder'

WINDOWS=100000
REF=symattis
QUERY=flame
FILTER=400

python $scrip_dir/HapGMMbuild.py \
	-db $db_dir/${QUERY}_${REF}_${WINDOWS}.tsv.gz \
	-rf ${REF} \
	-qr ${QUERY} \
	-w ${WINDOWS} \
	-vf ${FILTER} \
	-nc 3 \
	-gmf $out_dir/gmm_${REF}_vs_${QUERY}_${WINDOWS}_wd_${FILTER}_flt.tsv \
	-hr $out_dir/whole_genome_variations_histogram_${REF}_vs_${QUERY}_${WINDOWS}_wd.pdf \
	-hg $out_dir/gmm_categorized_histogram_${REF}_vs_${QUERY}_${WINDOWS}_wd_${FILTER}_flt.pdf