#!/bin/bash

REF=Jagger
QUERY=Flame
WINDOWS=500000
FILTER=1800
CHROM=ch1A

scrip_dir='path/to/scripts'
in_dir='path/to/database'
out_dir='path/to/output_folder'


python $scrip_dir/HapGMMbuild.v2.py \
	-db $in_dir/kmeribs-Wheat_${REF}-${QUERY}_${WINDOWS}.tsv.gz \
	-rf ${REF} \
	-qr ${QUERY} \
	-w ${WINDOWS} \
	-vf ${FILTER} \
	-nc 3 \
	-chr $CHROM \
	-st 3 \
	-gm $out_dir/gmm_${REF}_vs_${QUERY}_${WINDOWS}_wd_${FILTER}_flt.tsv.gz \
	-ho $out_dir/${CHROM}_${REF}_vs_${QUERY}_${WINDOWS}_wd_${FILTER}_flt.tsv.gz 