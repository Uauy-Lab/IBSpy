import argparse
import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture


def read_db():
	inDB = pd.read_csv(args.db_file, delimiter='\t')
	return inDB

#### prepare input data ####
def prepare_data(args):
	inDB = read_db()
	filtered_inDB = inDB[inDB['variations'] <= int(args.varfltr)]
	varDF = pd.DataFrame(filtered_inDB[['seqname','start','variations']])
	varDF.reset_index(drop=True, inplace=True)
	varArray = np.array(varDF['variations'])
	log_varArray = np.log(varArray, where=(varArray != 0))
	log_varArray = np.nan_to_num(log_varArray, nan=0.0)
	log_varArray = log_varArray.reshape((len(log_varArray),1))
	return log_varArray, varDF
	
#### build GMM model ####
def GMM_build(args):
	comNumber = int(args.comNumber)
	model = GaussianMixture(n_components=comNumber, covariance_type='full')
	log_varArray, varDF = prepare_data(args)
	model.fit(log_varArray)
	varGMM_predict = model.predict(log_varArray)
	varGMM_predict = varGMM_predict.reshape((len(varGMM_predict),1)) 
	varGMM_predict = pd.DataFrame({'gauss_mx_IBS': varGMM_predict[:, 0]})     
	varDF['gauss_mx_IBS'] = varGMM_predict['gauss_mx_IBS']                                   
	gmmIBS_DF = varDF.groupby(by=['gauss_mx_IBS'])
	return gmmIBS_DF, varDF

def transform_gmm(args):
	gmmIBS_DF, varDF = GMM_build(args)
	key_group = []
	value_group = []
	for key, group in gmmIBS_DF:
	    key_group.append(key)
	    value_group.append(group['variations'].median())
	dic_grupby = dict(zip(key_group, value_group))
	IBS_group = min(dic_grupby, key=dic_grupby.get)
	varDF['gauss_mx_IBS'] = np.where(varDF['gauss_mx_IBS'] == IBS_group, 1, 0)
	return varDF

def gmm_by_chromosome():
	inDB = transform_gmm(args)
	byChrDf = inDB[inDB['seqname'].str.contains(args.chrNme)].copy()
	byChrDf.reset_index(drop=True, inplace=True)
	return byChrDf

#### stitching haplotypes: number of non-IBS "outliers" that must a appear consecutively in a windows to be called non-IBS ####
def stitch_haploBlocks(args):
    byChrDf = gmm_by_chromosome()
    gmmDta = byChrDf['gauss_mx_IBS'].copy()
    StitchVarNum = int(args.StitchVarNum)
    for i in range(len(gmmDta)):
        if gmmDta[i] == 1:
            if i < (len(gmmDta) - StitchVarNum):
                count = 0
                for n in range(1, (StitchVarNum+1)):
                    if gmmDta[i+n] == 0:
                        count += 1
                if count == StitchVarNum:
                    continue
                else:
                    gmmDta[i+1] = 1
    hapBlock = gmmDta
    hapBlock = pd.Series(hapBlock)
    byChrDf['hap_block'] = hapBlock.values
    # put back full dataset of data for final haplo blocks
    hapCntFile = read_db()
    hapCntFile = hapCntFile[['seqname', 'start', 'variations']].copy()
    hapCntFile = hapCntFile[hapCntFile['seqname'].str.contains(args.chrNme)].copy()
    byChrDf = pd.merge(hapCntFile, byChrDf, left_on='start', right_on='start', how='left')
    byChrDf.loc[:,'gauss_mx_IBS':'hap_block'] = np.where(byChrDf.loc[:, 'gauss_mx_IBS':'hap_block'] == 1, 1, 0)
    byChrDf.rename(columns={'seqname_x':'seqname', 'variations_x':'variations'}, inplace=True)
    byChrDf = byChrDf[['seqname', 'start', 'variations', 'gauss_mx_IBS', 'hap_block']].copy()
    return byChrDf

def parse_arguments():
	parser = argparse.ArgumentParser()

	parser.add_argument('-db', '--db_file', help='Tab separated file with variations genetared by IBSpy output')
	parser.add_argument('-rf', '--refId', help='Name of the genome reference used')
	parser.add_argument('-qr', '--qryId', help='Name of the query sample')
	parser.add_argument('-chr', '--chrNme', help='Chromosome name to be plotted')
	parser.add_argument('-w', '--windSize', help='Windows size to count variations within')
	parser.add_argument('-vf', '--varfltr', help='Filter variations above this threshold to compute GMM model')
	parser.add_argument('-nc', '--comNumber', help='Number of componenets for the GMM model')
	parser.add_argument('-st', '--StitchVarNum', help='Consecutive "outliers" in windows to stitch')
	parser.add_argument('-gm','--output', help='GMM model output file grouped by IBS (1) and non-IBS (0)')
	parser.add_argument('-ho','--hap_output', help='GMM model output file grouped by IBS (1) and non-IBS (0)')
	
	args = parser.parse_args()
	return args

args = parse_arguments()

transform_gmm(args).to_csv(args.output, index=False,  sep='\t', compression="gzip")
stitch_haploBlocks(args).to_csv(args.hap_output, index=False,  sep='\t', compression="gzip")

