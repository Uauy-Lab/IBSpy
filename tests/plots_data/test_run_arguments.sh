# minimal arguments
IBSplot --IBSpy_counts "kmeribs-Wheat_Jagger-Flame.tsv.gz" \
--window_size 400000 \
--output gmm_ibs.tsv.gz \
--reference Jagger \
--query Flame \
--plot_output gmm_plots.pdf

# all arguments
IBSplot --IBSpy_counts "kmeribs-Wheat_Jagger-Flame.tsv.gz" \
--window_size 400000 \
--output gmm_ibs.tsv.gz \
--reference Jagger \
--query Flame \
--plot_output gmm_plots.pdf \
--filter_counts 1000 \
--n_components 3 \
--covariance_type 'full' \
--stitch_number 3