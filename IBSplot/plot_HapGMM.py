import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, NullFormatter, ScalarFormatter)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Build haplotypes and make scatter plot for vizualization")

    db_p = parser.add_argument_group('Sample name parameters')
    db_p.add_argument('-gm', '--gmm_file',  required=True, help='Tab separated file with IBS and non_IBS categorized by GMM model')
    db_p.add_argument('-db', '--db_file',  required=True, help='Tab separated file with variations genetared by IBSpy output')
    db_p.add_argument('-rf', '--refId',  required=True, help='Name of the genome reference used')
    db_p.add_argument('-qr', '--qryId',  required=True, help='Name of the query sample')
    db_p.add_argument('-chr', '--chrNme',  required=True, help='Chromosome name to be plotted')
    db_p.add_argument('-cl', '--chr_length_file',  required=True, help='Reference chromosome lenghts file')

    hap_block_p = parser.add_argument_group('Haplotype blocks parameters')
    hap_block_p.add_argument('-w', '--windSize',  required=True, help='Windows size to count variations within')
    hap_block_p.add_argument('-vf', '--varfltr',  required=True, help='Filter variations above this threshold to compute GMM model')
    hap_block_p.add_argument('-st', '--StitchVarNum',  required=True, help='Stitching haplotypes: number of non-IBS "outliers" that must a appear consecutively in a windows to be called non-IBS')

    out_files = parser.add_argument_group('Output files')
    out_files.add_argument('-o','--out_img_file', help='Output scatter plot and bars with haplotypes in ".jpg" format ')

    args = parser.parse_args()

    gmm_file = args.gmm_file
    db_file = args.db_file
    refId = args.refId
    qryId = args.qryId
    chrNme = args.chrNme
    chr_length_file = args.chr_length_file
    windSize = int(args.windSize)
    varfltr = int(args.varfltr)
    StitchVarNum = int(args.StitchVarNum)
    out_img_file = args.out_img_file


#### GMM input file ####
'''
It is a temporary file where IBS and non-IBS were computed (to 1 as IBS, and 0 non-IBS)
using the number of variations found in a windows of a size employing the GMM model.
'''
inFile = pd.read_csv(gmm_file, delimiter='\t')
byChrDf = inFile[inFile['seqname'] == chrNme].copy()
byChrDf.reset_index(drop=True, inplace=True)

#### stitching haplotypes: number of non-IBS "outliers" that must a appear consecutively in a windows to be called non-IBS ####
def stitch_haplotypes(byChrDf):
    gmmDta = byChrDf['gauss_mx_IBS'].copy()
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
    hapBlock = pd.Series(hapBlock) # final haplotype block after stitching
    byChrDf['hap_block'] = hapBlock.values # new column for stitched haplotypes
    return byChrDf

#### put back the k-mer variations Data Base file with all variations  count per windows position ####
hapCntFile = pd.read_csv(db_file, delimiter='\t')
hapCntFile = hapCntFile[['seqname', 'start', 'variations']].copy()
hapCntFile = hapCntFile[hapCntFile['seqname'] == chrNme].copy()

byChrDf = stitch_haplotypes(byChrDf)
byChrDf = pd.merge(hapCntFile, byChrDf, left_on='start', right_on='start', how='left')
byChrDf.loc[:,'gauss_mx_IBS':'hap_block'] = np.where(byChrDf.loc[:, 'gauss_mx_IBS':'hap_block'] == 1, 1, 0)
byChrDf.rename(columns={'seqname_x':'seqname', 'variations_x':'variations'}, inplace=True)
byChrDf = byChrDf[['seqname', 'start', 'variations', 'gauss_mx_IBS', 'hap_block']].copy()

######## plots ########
#### GMM stitched #### 
IBS = byChrDf[byChrDf['hap_block'] == 1].copy() 
strtIbsCln = np.array(IBS['start'])
hapBlkIbsCln = np.array(IBS['hap_block'])
hapArrIbsCln = np.vstack((strtIbsCln, hapBlkIbsCln))

non_IBS = byChrDf[byChrDf['hap_block'] == 0].copy()
chr_st_non_IBS = np.array(non_IBS['start'])
hp_blk_non_IBS = np.array(non_IBS['hap_block'])
chr_hp_arr_non_IBS = np.vstack((chr_st_non_IBS, hp_blk_non_IBS))

#### GMM before stitching ####
IBS_org = byChrDf[byChrDf['gauss_mx_IBS'] == 1].copy()
strtIbsCln_org = np.array(IBS_org['start'])
IBS_org = np.array(IBS_org['gauss_mx_IBS'])
hapArrIbsCln_org = np.vstack((strtIbsCln_org, IBS_org))

non_IBS_org = byChrDf[byChrDf['gauss_mx_IBS'] == 0].copy()
chr_st_non_IBS_org = np.array(non_IBS_org['start'])
non_IBS_org = np.array(non_IBS_org['gauss_mx_IBS'])
chr_hp_arr_non_IBS_org = np.vstack((chr_st_non_IBS_org, non_IBS_org))


#### base subplots ####
fig, axs = plt.subplots(3, figsize=(20, 7), sharex=True, gridspec_kw={'height_ratios': [1, 1,3]})
fig.suptitle(f'{chrNme}: {refId} vs {qryId}, {windSize} bp windows, stitched var: {StitchVarNum}, filter: {varfltr}', fontweight="bold", fontsize=25)

#### horizontal bar plot ####
### GMM_stitched
axs[0].eventplot(hapArrIbsCln, colors='darkcyan', alpha=1, lineoffsets=[1,1], linelengths=1.5)
axs[0].eventplot(chr_hp_arr_non_IBS, colors='goldenrod', lineoffsets=[1,1], linelengths=1)

### GMM_original (before stitching) ###
axs[1].eventplot(hapArrIbsCln_org, colors='darkcyan', alpha=1, lineoffsets=[1,1], linelengths=1.5)
axs[1].eventplot(chr_hp_arr_non_IBS_org, colors='goldenrod', lineoffsets=[1,1], linelengths=1)

### scatter plot ###
IBS_block = byChrDf[(byChrDf['hap_block'] == 1)]
x = IBS_block['start']
y = IBS_block['variations']
axs[2].scatter(x,y, c="darkcyan", alpha=0.5, s=4, label='IBS')

non_IBS_block = byChrDf[(byChrDf['hap_block'] == 0)]
a = non_IBS_block['start']
b = non_IBS_block['variations']
axs[2].scatter(a,b, c="goldenrod", alpha=0.6, s=4, marker='>', label='non-IBS')
axs[2].set_yscale('symlog')
axs[2].set_ylabel('Variations count', fontweight="bold", fontsize=17)
axs[2].set_xlabel('Chromosome position [Mbp]', fontweight="bold", fontsize=17)
axs[2].yaxis.set_major_formatter(ScalarFormatter())

### get chromosome length to plot X_axis limits ###
chr_file = pd.read_csv(chr_length_file, delimiter='\t', header=None)
chr_length = int(chr_file[chr_file[0] == chrNme][1].values) + 10000000

#### labels and axis limits ####
names = ['gmm_stitched', 'gmm_origin']
for axs, names in zip(axs, names):
        axs.text(-0.01, 0.5, names, va='center', ha='right', fontweight= "bold", fontsize=17, transform=axs.transAxes)
        plt.grid(True, color='gray', linestyle='--', linewidth=0.5)
        
        axs.set_xlim(0,chr_length)
        x_labels = list(range(0,chr_length,20000000)) # one label every 20 Mb
        list_labels = int(chr_length/1000000) # transform to Mb to reduce number size
        x_ticklabels = list(range(0,list_labels,20)) # use tranformed labels every 20 Mb
        plt.xticks(x_labels, x_ticklabels, rotation=45, fontsize=15)
        axs.yaxis.set_major_formatter(ScalarFormatter())
        plt.yticks([1,2,3,10,20,40,100,300,1000,2000,5000], fontsize = 15)
        plt.legend(prop={'size': 10}, markerscale=4, borderaxespad=0.3, handletextpad=0.1)
        axs.set_axis_off()
fig.savefig(out_img_file, format='jpg', dpi=150, bbox_inches='tight')