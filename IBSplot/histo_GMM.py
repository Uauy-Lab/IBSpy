import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, NullFormatter, ScalarFormatter)

def read_db():
	inDB = pd.read_csv(args.db_file, delimiter='\t')
	return inDB

def read_gmm_f():
    in_gmm_f = pd.read_csv(args.gmm_File, delimiter='\t')
    return in_gmm_f

# data
def plot_histo_db(args):
	inDB = read_db()
	X_a = inDB.loc[inDB['variations']]
	X_a = inDB['variations'].astype(int)

	#plot
	fig = plt.figure(figsize=(10, 2))
	plt.hist(X_a, log=False, alpha=1, bins=range(min(X_a), max(X_a)), color='darkred', label='All-regions')
	plt.title(f'All variations [{args.refId} vs {args.qryId}], {args.windSize}_windows', fontweight="bold", fontsize=12)
	plt.legend(prop={'size': 15})
	plt.grid(True)
	plt.ylabel('Frequency', fontweight="bold", fontsize=11)
	plt.xlabel('Variations count', fontweight="bold", fontsize=11)
	plt.xscale('symlog')
	plt.xlim(0,5000)
	fig.axes[0].xaxis.set_major_formatter(ScalarFormatter())
	plt.xticks([0,1,2,3,5,10,20,40,100,300,1000,2000,5000], fontsize = 10)
	return fig

#### GMM variations count histogram after classification ####
# data
def plot_histo_gmm(args):
	varDF = read_gmm_f()
	X_b = varDF.loc[varDF['gauss_mx_IBS'] == 1]
	X_b = X_b['variations'].astype(int)
	X_c = varDF.loc[varDF['gauss_mx_IBS'] == 0]
	X_c = X_c['variations'].astype(int)

	# plot
	fig = plt.figure(figsize=(10, 2))
	plt.hist(X_b, alpha=0.5, bins=range(min(X_b), max(X_b)), color='darkblue', label='IBS')
	plt.hist(X_c, alpha=0.5, bins=range(min(X_c), max(X_c)), color='darkgreen', label='non-IBS')
	plt.title(f'GMM: filtered variations <= {args.varfltr} [{args.refId} vs {args.qryId}], {args.windSize}_windows', fontweight="bold", fontsize=12)
	plt.legend(prop={'size': 15})
	plt.grid(True)
	plt.ylabel('Frequency', fontweight="bold", fontsize=11)
	plt.xlabel('Variations count', fontweight="bold", fontsize=11)
	plt.xscale('symlog')
	plt.xlim(0,5000)
	fig.axes[0].xaxis.set_major_formatter(ScalarFormatter())
	plt.xticks([0,1,2,3,5,10,20,40,100,300,1000,2000,5000], fontsize=10)
	return fig

def parse_arguments():
	parser = argparse.ArgumentParser()

	parser.add_argument('-db', '--db_file', help='Tab separated file with variations genetared by IBSpy output')
	parser.add_argument('-gm', '--gmm_File', help='GMM model output file grouped by IBS (1) and non-IBS (0)')
	parser.add_argument('-rf', '--refId', help='Name of the genome reference')
	parser.add_argument('-qr', '--qryId', help='Name of the query sample')
	parser.add_argument('-w', '--windSize', help='Windows size to count variations within')
	parser.add_argument('-vf', '--varfltr', help='Filter variations above this threshold to compute GMM model')
	parser.add_argument('-do','--db_output', help='Output histogram of the variations from the whole genome in ".pdf" format')
	parser.add_argument('-go','--gmm_output', help='Output histogram of the variations from the whole genome by IBS and non-IBS in ".pdf" format')

	args = parser.parse_args()
	return args

args = parse_arguments()

plot_histo_db(args).savefig(args.db_output, format='pdf', bbox_inches='tight')
plot_histo_gmm(args).savefig(args.gmm_output, format='pdf', bbox_inches='tight')
