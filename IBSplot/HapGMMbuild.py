import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from matplotlib.ticker import (MultipleLocator, NullFormatter, ScalarFormatter)


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description = "build GMM mode for IBS and non-IBS classification using variation count per windows")

	db_p = parser.add_argument_group('Sample name parameters')
	db_p.add_argument('-db', '--db_file',  required=True, help='Tab separated file with variations genetared by IBSpy output')
	db_p.add_argument('-rf', '--refId',  required=True, help='Name of the genome reference used')
	db_p.add_argument('-qr', '--qryId',  required=True, help='Name of the query sample')

	gmm_p = parser.add_argument_group('GMM model parameters')
	gmm_p.add_argument('-w', '--windSize',  required=True, help='Windows size to count variations within')
	gmm_p.add_argument('-vf', '--varfltr',  required=True, help='Filter variations above this threshold to compute GMM model')
	gmm_p.add_argument('-nc', '--comNumber',  required=True, help='Number of componenets for the GMM model')

	out_files = parser.add_argument_group('Output files')
	out_files.add_argument('-gmf','--output', help='GMM model output file grouped by IBS (1) and non-IBS (0)')
	out_files.add_argument('-hr','--out_hist_raw', help='Create histogram distribution of the whole genome using variations count raw data per windows')
	out_files.add_argument('-hg','--out_hist_gmm', help='Histogram of the GMM grouped by IBS and non-IBS')
	
#### command line arguments ####
	args = parser.parse_args()
	
	db_file = args.db_file
	refId = args.refId
	qryId = args.qryId
	windSize = args.windSize
	varfltr = int(args.varfltr)

#### prepare the input data ####
inFile = pd.read_csv(db_file, delimiter='\t')
filtered_inFile = inFile[inFile['variations'] <= varfltr]
varDF = pd.DataFrame(filtered_inFile[['seqname','start','variations']])
varDF.reset_index(drop=True, inplace=True)
varArray = np.array(varDF['variations'])
log_varArray = np.log(varArray, where=(varArray != 0))
log_varArray = np.nan_to_num(log_varArray, nan=0.0)
log_varArray = log_varArray.reshape((len(log_varArray),1))


#### build GMM model ####
def GMM_build(comNumber, output):
	model = GaussianMixture(n_components=comNumber, covariance_type='full')
	model.fit(log_varArray)
	varGMM_predict = model.predict(log_varArray)
	varGMM_predict = varGMM_predict.reshape((len(varGMM_predict),1)) 
	varGMM_predict = pd.DataFrame({'gauss_mx_IBS': varGMM_predict[:, 0]})         
	varDF['gauss_mx_IBS'] = varGMM_predict['gauss_mx_IBS']                                   
	get_IBS = varDF.groupby(by=['gauss_mx_IBS'])
	key_group = []
	value_group = []
	for key, group in get_IBS:
	    key_group.append(key)
	    value_group.append(group['variations'].median())
	dic_grupby = dict(zip(key_group, value_group))
	IBS_group = min(dic_grupby, key=dic_grupby.get)
	varDF['gauss_mx_IBS'] = np.where(varDF['gauss_mx_IBS'] == IBS_group, 1, 0)
	varDF.to_csv(output, index=False,  sep='\t')


################### plots #################
#### all variations count histogram before GMM model classification ####
# data
def plot_histo_raw(db_file, out_hist_raw):
	varCntFile = pd.read_csv(db_file, delimiter='\t')
	X_a = varCntFile.loc[varCntFile['variations']]
	X_a = varCntFile['variations'].astype(int)

	#plot
	fig = plt.figure(figsize=(10, 2))
	plt.hist(X_a, log=False, alpha=1, bins=range(min(X_a), max(X_a)), color='darkred', label='All-regions')
	plt.title(f'All variations [{refId} vs {qryId}], {windSize}_windows', fontweight="bold", fontsize=12)
	plt.legend(prop={'size': 15})
	plt.grid(True)
	plt.ylabel('Frequency', fontweight="bold", fontsize=11)
	plt.xlabel('Variations count', fontweight="bold", fontsize=11)
	plt.xscale('symlog')
	plt.xlim(0,5000)
	fig.axes[0].xaxis.set_major_formatter(ScalarFormatter())
	plt.xticks([0,1,2,3,5,10,20,40,100,300,1000,2000,5000], fontsize = 10)
	fig.savefig(out_hist_raw, format='pdf', bbox_inches='tight')

#### GMM variations count histogram after classification ####
# data
def plot_histo_gmm(varDF, varfltr, out_hist_gmm):
	X_b = varDF.loc[varDF['gauss_mx_IBS'] == 1]
	X_b = X_b['variations'].astype(int)
	X_c = varDF.loc[varDF['gauss_mx_IBS'] == 0]
	X_c = X_c['variations'].astype(int)

	# plot
	fig = plt.figure(figsize=(10, 2))
	plt.hist(X_b, alpha=0.5, bins=range(min(X_b), max(X_b)), color='darkblue', label='IBS')
	plt.hist(X_c, alpha=0.5, bins=range(min(X_c), max(X_c)), color='darkgreen', label='non-IBS')
	plt.title(f'GMM: filtered variations count <= {varfltr} [{refId} vs {qryId}], {windSize}_windows', fontweight="bold", fontsize=12)
	plt.legend(prop={'size': 15})
	plt.grid(True)
	plt.ylabel('Frequency', fontweight="bold", fontsize=11)
	plt.xlabel('Variations count', fontweight="bold", fontsize=11)
	plt.xscale('symlog')
	plt.xlim(0,5000)
	fig.axes[0].xaxis.set_major_formatter(ScalarFormatter())
	plt.xticks([0,1,2,3,5,10,20,40,100,300,1000,2000,5000], fontsize=10)
	fig.savefig(out_hist_gmm, format='pdf', bbox_inches='tight')

#### Write out files ####
if __name__ == '__main__':
	args = parser.parse_args()
	GMM_build(int(args.comNumber), args.output)
	plot_histo_raw(args.db_file, args.out_hist_raw)
	plot_histo_gmm(varDF, int(args.varfltr), args.out_hist_gmm)
