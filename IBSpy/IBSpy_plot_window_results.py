import argparse
from IBSpy import IBSpyResults
def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--IBSpy_counts',  required=True, help='Tab separated file with variations genetared by IBSpy output')
    parser.add_argument('-c', '--chromosome_length',  required=True, help='Reference chromosome lenghts file')
    parser.add_argument('-w', '--window_size',  required=True, help='Windows size to count variations within')
    parser.add_argument('-o', '--output', help='File with variations by windows only')

    args = parser.parse_args()
    return args

def main():
	args = parse_arguments()
	results = IBSpyResults(args.IBSpy_counts, args.chromosome_length, args.window_size)
	results_pd = results.count_by_windows()
	results_pd.to_csv(args.output, index=False,  sep='\t', compression="gzip")

if __name__ == '__main__':
	main()