# all arguments
IBSplot -i kmeribs-Wheat_Jagger-Flame.tsv.gz \
-w 500000 \
-f 1000 \
-n 3 \
-c full \
-s 3 \
-o gmm_ibs.tsv.gz \
-r Jagger \
-q Flame \
-p gmm_plots.pdf


# minimal arguments
IBSplot -i kmeribs-Wheat_Jagger-Flame.tsv.gz \
-w 500000 \
-o gmm_ibs.tsv.gz \
-r Jagger \
-q Flame \
-f 1000 \
-p gmm_plots.pdf
