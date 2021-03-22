import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

class IBSpyPlots:
	# class variables go here
	figure_size = (20, 5)
	X_hist_lim = (0,10000)
	v_count_ticks = [0,1,2,5,10,20,50,100,200,400,1000,2000,5000,10000]
	ax_font_size = 17
	axx_tick_size = 10
	axy_tick_size = 15

	def __init__(self, gmm_results, window_size, filter_counts, r_name, q_name):
		
		self.gmm_results = gmm_results
		self.window_size = window_size
		self.filter_counts = filter_counts
		self.r_name = r_name
		self.q_name = q_name

	# raw_data
	def raw_histogram(self):
		inDB = self.gmm_results
		X_a = inDB['variations']

		fig = plt.figure(figsize= (20, 4))
		plt.hist(X_a, log=False, alpha=1, bins=range(min(X_a), max(X_a)), color='darkred', label='All-counts')
		plt.title(f'Raw variations count, whole genome: {self.r_name} vs {self.q_name}, {self.window_size}_windows', fontweight="bold", fontsize=22)
		plt.legend(prop={'size': 15})
		plt.grid(True)
		plt.ylabel('Frequency', fontweight="bold", fontsize=self.ax_font_size)
		plt.xlabel('Variations count', fontweight="bold", fontsize=self.ax_font_size)
		plt.xscale('symlog')
		plt.xlim(self.X_hist_lim)
		fig.axes[0].xaxis.set_major_formatter(ScalarFormatter())
		plt.xticks(self.v_count_ticks, fontsize = self.axx_tick_size)
		plt.yticks(fontsize = self.axy_tick_size)
		return fig
	
	# GMM
	def gmm_histogram(self):
		inDB = self.gmm_results
		X_b = inDB['variations'][inDB['v_gmm'] == 1]
		X_c = inDB['variations'][inDB['v_gmm'] == 0]

		fig = plt.figure(figsize=(20, 4))
		plt.hist(X_b, alpha=0.5, bins=range(min(X_b), max(X_b)), color='darkblue', label='IBS')
		plt.hist(X_c, alpha=0.5, bins=range(min(X_c), max(X_c)), color='darkgreen', label='non-IBS')
		plt.title(f'GMM variations count <= {self.filter_counts}, whole genome: {self.r_name} vs {self.q_name}, {self.window_size}_windows', fontweight="bold", fontsize=22)
		plt.legend(prop={'size': 15})
		plt.grid(True)
		plt.ylabel('Frequency', fontweight="bold", fontsize=self.ax_font_size)
		plt.xlabel('Variations count', fontweight="bold", fontsize=self.ax_font_size)
		plt.xscale('symlog')
		plt.xlim(self.X_hist_lim)
		fig.axes[0].xaxis.set_major_formatter(ScalarFormatter())
		plt.xticks(self.v_count_ticks, fontsize=self.axx_tick_size)
		plt.yticks(fontsize = self.axy_tick_size)
		return fig

	def scatterplot_by_chromosome(self, chromosome_id):
		variations_db = self.gmm_results
		variations_db = variations_db[variations_db['seqname'].str.contains(chromosome_id)]
		
		IBS = variations_db[variations_db['vh_block'] == 1]
		IBS_window_position = np.array(IBS['window'])
		IBS_haplo_block = np.array(IBS['vh_block'])
		IBS_window_block_array = np.vstack((IBS_window_position, IBS_haplo_block))

		non_IBS = variations_db[variations_db['vh_block'] == 0]
		non_IBS_window_position = np.array(non_IBS['window'])
		non_IBS_haplo_block = np.array(non_IBS['vh_block'])
		non_IBS_window_block_array = np.vstack((non_IBS_window_position, non_IBS_haplo_block))

		# base subplots
		fig, axs = plt.subplots(2, figsize=(20, 4), sharex=True, gridspec_kw={'height_ratios': [1, 3]})
		plt.subplots_adjust(hspace=False)
		fig.suptitle(f'{chromosome_id}: {self.r_name} vs {self.q_name}, {self.window_size} bp windows, filter: {self.filter_counts}', fontweight="bold", fontsize=22)

		# horizontal bar plot
		axs[0].eventplot(IBS_window_block_array, colors='darkcyan', alpha=1, lineoffsets=[1,1], linelengths=1.5)
		axs[0].eventplot(non_IBS_window_block_array, colors='goldenrod', lineoffsets=[1,1], linelengths=1)

		# scatter plot
		x = IBS['window']
		y = IBS['variations']
		axs[1].scatter(x,y, c="darkcyan", alpha=0.5, s=4, label='IBS')
		a = non_IBS['window']
		b = non_IBS['variations']
		axs[1].scatter(a,b, c="goldenrod", alpha=0.6, s=4, marker='>', label='non-IBS')

		axs[1].set_yscale('symlog')
		axs[1].set_ylabel('Variations count', fontweight="bold", fontsize = self.ax_font_size)
		axs[1].set_xlabel('Chromosome position [Mbp]', fontweight="bold", fontsize = self.ax_font_size)
		axs[1].yaxis.set_major_formatter(ScalarFormatter())
		axs[1].set_ylim(1,10000)

		# get chromosome length to plot X_axis limits
		chr_length = variations_db['window'].max() + 2000000
		names = ['GMM']
		for axs, names in zip(axs, names):
			axs.text(-0.01, 0.5, names, va='center', ha='right', fontweight= "bold", fontsize = self.ax_font_size, transform=axs.transAxes)
			plt.grid(True, color='gray', linestyle='--', linewidth=0.5)
			axs.set_xlim(0,chr_length)
			x_labels = list(range(0,chr_length,10000000)) # one label every 20 Mb
			list_labels = int(chr_length/1000000) # transform to Mb to reduce number size
			x_ticklabels = list(range(0,list_labels,10)) 
			plt.xticks(x_labels, x_ticklabels, rotation=45, fontsize = self.axx_tick_size)
			plt.yticks(self.v_count_ticks, fontsize = self.axy_tick_size)
			plt.legend(prop={'size': 10}, markerscale=4, borderaxespad=0.3, handletextpad=0.1)
			axs.set_axis_off()
		return fig
