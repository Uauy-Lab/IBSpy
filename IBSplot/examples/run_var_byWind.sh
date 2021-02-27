#!/bin/bash

REF=Jagger
QUERY=Flame
IN_WINDOWS=50000
OUT_WINDOWS=500000

scrip_dir='path/to/scripts'
in_dir='path/to/database'
out_dir='path/to/output_folder'
chr_dir='path/to/genome_sizes_folder'

python $scrip_dir/var_byWind.v2.py \
	-db $in_dir/kmeribs-Wheat_${REF}-${QUERY}.tsv.gz.gz \
	-cl $chr_dir/chr_sizes_Jagger.genome.txt \
	-w ${OUT_WINDOWS} \
	-o $out_dir/kmeribs-Wheat_${REF}-${QUERY}_${OUT_WINDOWS}.tsv.gz