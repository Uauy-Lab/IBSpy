import argparse
from IBSpy import IBSpyResults
# import IBSpyResults
def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--IBSpy_counts',  required=True, help='Tab separated file with variations genetared by IBSpy output')
    parser.add_argument('-c', '--chromosome_length',  required=True, help='Reference chromosome lenghts file')
    parser.add_argument('-w', '--window_size',  required=True, help='Windows size to count variations within')
    parser.add_argument('-f', '--filter_counts', help='Filter variations above this threshold to compute GMM model')
    parser.add_argument('-n', '--n_components', help='Number of componenets for the GMM model,default=3 ')
    parser.add_argument('-v', '--covariance_type', help='type of covariance used for GMM model, default, "full"')
    parser.add_argument('-s', '--stitch_number', help='Consecutive "outliers" in windows to stitch, default=3')
    parser.add_argument('-o', '--output', help='File with variations by windows only')

    args = parser.parse_args()
    return args

def main():
	args = parse_arguments()
	results = IBSpyResults(args.IBSpy_counts, args.chromosome_length, args.window_size, args.filter_counts)
    # results.fit_gmm_model(args.n_components, args.covariance_type)
    results = results.stitch_haploBlocks(args.stitch_number)
	results_pd = results.stitch_haploBlocks()

	results_pd.to_csv(args.output, index=False,  sep='\t', compression="gzip")

if __name__ == '__main__':
	main()