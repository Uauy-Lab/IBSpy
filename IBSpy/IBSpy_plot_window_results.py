import argparse
from IBSpy import IBSpyResults, IBSpyPlots
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--IBSpy_counts', help='tvs file genetared by IBSpy output')
    parser.add_argument('-w', '--window_size', type=int, help='Windows size to count variations within')
    parser.add_argument('-f', '--filter_counts', default=None, type=int, help='Filter number of variaitons above this threshold to compute GMM model, default=None')
    parser.add_argument('-n', '--n_components', default=3, type=int, help='Number of componenets for the GMM model, default=3')
    parser.add_argument('-c', '--covariance_type', default='full', help='type of covariance used for GMM model, default="full"')
    parser.add_argument('-s', '--stitch_number', default=3, type=int, help='Consecutive "outliers" in windows to stitch, default=3')
    parser.add_argument('-o', '--output', help='tsv file with variations count by windows and summary statistics')

    # plot arguments
    parser.add_argument('-r', '--reference', help='genome reference name')
    parser.add_argument('-q', '--query', help='query sample')
    parser.add_argument('-p','--plot_output', help=' histograms and ascatter files in .PDF format')

    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()

    # count variations by windows from IBS_results
    IBS_results = IBSpyResults(args.IBSpy_counts, args.window_size, args.filter_counts)
    IBS_results = IBS_results.run_analysis(args.n_components, args.covariance_type, args.stitch_number)
    IBS_results.to_csv(args.output, index=False,  sep='\t', compression="gzip")

    plots =  IBSpyPlots(IBS_results, args.window_size, args.filter_counts, args.reference, args.query)
    with PdfPages(args.plot_output) as pdf:

        # histograms
        raw_hist = plots.raw_histogram()
        pdf.savefig(bbox_inches='tight')

        gmm_hist = plots.gmm_histogram()
        pdf.savefig(bbox_inches='tight')
        
        plt.close(raw_hist)
        plt.close(gmm_hist)

        # scatter plot
        chromosomes = IBS_results["seqname"].str.split('_').str[0].unique()
        for chromosome in chromosomes:
            scatter_by_chr = plots.scatterplot_by_chromosome(chromosome)
            pdf.savefig(bbox_inches='tight')
            plt.close(scatter_by_chr)

if __name__ == '__main__':
	main()