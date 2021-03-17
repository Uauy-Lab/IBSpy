import argparse
from IBSpy import IBSpyResults, IBSpyHistoPlot

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--IBSpy_counts',  required=True, help='Tab separated file with variations genetared by IBSpy output')
    parser.add_argument('-c', '--chromosome_length',  required=True, help='Reference chromosome lenghts file')
    parser.add_argument('-w', '--window_size',  required=True, type=int, help='Windows size to count variations within')
    parser.add_argument('-f', '--filter_counts', type=int, help='Filter variations above this threshold to compute GMM model')
    parser.add_argument('-n', '--n_components', type=int, help='Number of componenets for the GMM model,default=3 ')
    parser.add_argument('-v', '--covariance_type', help='type of covariance used for GMM model, default, "full"')
    parser.add_argument('-s', '--stitch_number', type=int, help='Consecutive "outliers" in windows to stitch, default=3')
    parser.add_argument('-o', '--output', help='File with variations by windows only')

    # specific histo plot arguments
    # parser.add_argument('-g', '--g_filename', help='GMM model output file grouped by IBS (1) and non-IBS (0)')
    parser.add_argument('-r', '--reference_name', help='Name of the genome reference')
    parser.add_argument('-q', '--query_name', help='Name of the query sample')
    parser.add_argument('-ro','--raw_histogram', help='Output histogram of the variations from the whole genome in ".pdf" format')
    parser.add_argument('-go','--gmm_histogram', help='Output histogram of the variations from the whole genome by IBS and non-IBS in ".pdf" format')

    args = parser.parse_args()
    return args

# args = parse_arguments()
# plot_haploBlocks(args).savefig(args.outImg, format='jpg', dpi=150, bbox_inches='tight')


def main():
    args = parse_arguments()

    # File with variation count by windows, GMM and summary statistics
    results = IBSpyResults(args.IBSpy_counts, args.chromosome_length, args.window_size, args.filter_counts)
    results = results.stitch_haploBlocks(args.n_components, args.covariance_type, args.stitch_number)
    results.to_csv(args.output, index=False,  sep='\t', compression="gzip")

    plots =  IBSpyHistoPlot(results, args.window_size, args.filter_counts, args.reference_name, args.query_name)
    # Raw data histogram
    histo_ = plots.raw_data_histo()
    histo_.savefig(args.raw_histogram, format='jpg', dpi=150, bbox_inches='tight')
    # GMM data histogram
    histo_2 = plots.gmm_data_histo()
    histo_2.savefig(args.gmm_histogram, format='jpg', dpi=150, bbox_inches='tight')



    # plots = IBSpyChromosomeScatterPlot(results)
    # chromosomes=results["seqname"].unique()
    # plot_path=args.output_path + "_plot_"
    # for chr in chromosomes:
    #     plots.plot(chr, plot_path+chr+".pdf")



if __name__ == '__main__':
	main()