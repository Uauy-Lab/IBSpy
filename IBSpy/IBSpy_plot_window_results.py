import argparse
from IBSpy import IBSpyResults, IBSpyPlots
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--IBSpy_counts', help='tvs file genetared by IBSpy output')
    parser.add_argument('-w', '--window_size', type=int, help='Windows size to count variations within')
    parser.add_argument('-f', '--filter_counts', default=None, type=int, help='Filter number of variaitons above this threshold to compute GMM model')
    parser.add_argument('-n', '--n_components', default=3, type=int, help='Number of componenets for the GMM model')
    parser.add_argument('-c', '--covariance_type', default='full', help='type of covariance used for GMM model')
    parser.add_argument('-s', '--stitch_number', default=3, type=int, help='Consecutive "outliers" in windows to stitch')
    parser.add_argument('-o', '--output', help='tsv file with variations count by windows and summary statistics')

    # plot arguments
    parser.add_argument('-r', '--reference_name', help='genome reference name')
    parser.add_argument('-q', '--query_name', help='query sample')
    parser.add_argument('-p','--ouput_plot', help=' histograms and ascatter files in .PDF format')

    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()

    # count variations by windows from IBS_results
    IBS_results = IBSpyResults(args.IBSpy_counts, args.window_size, args.filter_counts)
    IBS_results = IBS_results.stitch_gmm_haplotypes(args.n_components, args.covariance_type, args.stitch_number)
    IBS_results.to_csv(args.output, index=False,  sep='\t', compression="gzip")


    plots =  IBSpyPlots(IBS_results, args.window_size, args.filter_counts, args.reference_name, args.query_name)
    with PdfPages(args.ouput_plot) as pdf:

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