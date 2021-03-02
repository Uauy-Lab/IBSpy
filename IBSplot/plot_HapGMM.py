import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, NullFormatter, ScalarFormatter)

def read_hapBlocks_f():
    hoFile = pd.read_csv(args.hoFile, delimiter='\t')
    return hoFile
def read_chrLen_f():
    in_chrLen = pd.read_csv(args.chrLen_f, delimiter='\t', header=None)
    return in_chrLen

def plot_haploBlocks(args):
    hoFile = read_hapBlocks_f()
    IBS = hoFile[hoFile['hap_block'] == 1]
    strtIbsCln = np.array(IBS['start'])
    hapBlkIbsCln = np.array(IBS['hap_block'])
    hapArrIbsCln = np.vstack((strtIbsCln, hapBlkIbsCln))

    non_IBS = hoFile[hoFile['hap_block'] == 0].copy()
    chr_st_non_IBS = np.array(non_IBS['start'])
    hp_blk_non_IBS = np.array(non_IBS['hap_block'])
    chr_hp_arr_non_IBS = np.vstack((chr_st_non_IBS, hp_blk_non_IBS))

    #### GMM before stitching ####
    IBS_org = hoFile[hoFile['gauss_mx_IBS'] == 1].copy()
    strtIbsCln_org = np.array(IBS_org['start'])
    IBS_org = np.array(IBS_org['gauss_mx_IBS'])
    hapArrIbsCln_org = np.vstack((strtIbsCln_org, IBS_org))

    non_IBS_org = hoFile[hoFile['gauss_mx_IBS'] == 0].copy()
    chr_st_non_IBS_org = np.array(non_IBS_org['start'])
    non_IBS_org = np.array(non_IBS_org['gauss_mx_IBS'])
    chr_hp_arr_non_IBS_org = np.vstack((chr_st_non_IBS_org, non_IBS_org))

    #### base subplots ####
    fig, axs = plt.subplots(3, figsize=(20, 7), sharex=True, gridspec_kw={'height_ratios': [1, 1,3]})
    fig.suptitle(f'{args.chrNme}: {args.refId} vs {args.qryId}, {args.windSize} bp windows, filter: {args.varfltr}', fontweight="bold", fontsize=25)

    #### horizontal bar plot ####
    ### GMM_stitched
    axs[0].eventplot(hapArrIbsCln, colors='darkcyan', alpha=1, lineoffsets=[1,1], linelengths=1.5)
    axs[0].eventplot(chr_hp_arr_non_IBS, colors='goldenrod', lineoffsets=[1,1], linelengths=1)

    ### GMM_original (before stitching) ###
    axs[1].eventplot(hapArrIbsCln_org, colors='darkcyan', alpha=1, lineoffsets=[1,1], linelengths=1.5)
    axs[1].eventplot(chr_hp_arr_non_IBS_org, colors='goldenrod', lineoffsets=[1,1], linelengths=1)

    ### scatter plot ###
    IBS_block = hoFile[(hoFile['hap_block'] == 1)]
    x = IBS_block['start']
    y = IBS_block['variations']
    axs[2].scatter(x,y, c="darkcyan", alpha=0.5, s=4, label='IBS')

    non_IBS_block = hoFile[(hoFile['hap_block'] == 0)]
    a = non_IBS_block['start']
    b = non_IBS_block['variations']
    axs[2].scatter(a,b, c="goldenrod", alpha=0.6, s=4, marker='>', label='non-IBS')
    axs[2].set_yscale('symlog')
    axs[2].set_ylabel('Variations count', fontweight="bold", fontsize=17)
    axs[2].set_xlabel('Chromosome position [Mbp]', fontweight="bold", fontsize=17)
    axs[2].yaxis.set_major_formatter(ScalarFormatter())

    ### get chromosome length to plot X_axis limits ###
    chrLen_f = read_chrLen_f()
    chr_length = int(chrLen_f[chrLen_f[0] == args.chrNme][1].values) + 10000000

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
    return fig

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-ho', '--hoFile', help='Tab separated file with IBS and non_IBS categorized by GMM model')
    parser.add_argument('-rf', '--refId', help='Name of the genome reference used')
    parser.add_argument('-qr', '--qryId', help='Name of the query sample')
    parser.add_argument('-chr', '--chrNme', help='Chromosome name to be plotted')
    parser.add_argument('-cl', '--chrLen_f', help='Reference chromosome lenghts file')
    parser.add_argument('-w', '--windSize', help='Windows size to count variations within')
    parser.add_argument('-vf', '--varfltr', help='Filter variations above this threshold to compute GMM model')
    parser.add_argument('-o','--outImg', help='Output scatter plot and bars with haplotypes in ".jpg" format')
    
    args = parser.parse_args()
    return args

args = parse_arguments()
plot_haploBlocks(args).savefig(args.outImg, format='jpg', dpi=150, bbox_inches='tight')
