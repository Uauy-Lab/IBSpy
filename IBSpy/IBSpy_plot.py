import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, NullFormatter, ScalarFormatter)

class IBSpyHistoPlot:
# class variables go here
	figure_size = (10, 2)

	def __init__(self, filename, window_size, filter_counts, r_name, q_name):
		#instance variables
		self.db = filename
		self.window_size = window_size
		self.filter_counts = filter_counts
		self.r_name = r_name
		self.q_name = q_name

	# Raw_data
	def raw_data_histo(self):
		inDB = self.db
		X_a = inDB['variations']

		#Raw variation counts plot
		fig = plt.figure(figsize=self.figure_size)
		plt.hist(X_a, log=False, alpha=1, bins=range(min(X_a), max(X_a)), color='darkred', label='All-counts')
		plt.title(f'Raw variations count [{self.r_name} vs {self.q_name}], {self.window_size}_windows', fontweight="bold", fontsize=12)
		plt.legend(prop={'size': 15})
		plt.grid(True)
		plt.ylabel('Frequency', fontweight="bold", fontsize=11)
		plt.xlabel('Variations count', fontweight="bold", fontsize=11)
		plt.xscale('symlog')
		plt.xlim(0,5000)
		fig.axes[0].xaxis.set_major_formatter(ScalarFormatter())
		plt.xticks([0,1,2,3,5,10,20,40,100,300,1000,2000,5000], fontsize = 10)
		return fig
	
	# GMM
	def gmm_data_histo(self):
		inDB = self.db
		X_b = inDB['variations'][inDB['v_gmm'] == 1]
		X_c = inDB['variations'][inDB['v_gmm'] == 0]

		fig = plt.figure(figsize=self.figure_size)
		plt.hist(X_b, alpha=0.5, bins=range(min(X_b), max(X_b)), color='darkblue', label='IBS')
		plt.hist(X_c, alpha=0.5, bins=range(min(X_c), max(X_c)), color='darkgreen', label='non-IBS')
		plt.title(f'GMM variations count <= {self.filter_counts} [{self.r_name} vs {self.q_name}], {self.window_size}_windows', fontweight="bold", fontsize=12)
		plt.legend(prop={'size': 15})
		plt.grid(True)
		plt.ylabel('Frequency', fontweight="bold", fontsize=11)
		plt.xlabel('Variations count', fontweight="bold", fontsize=11)
		plt.xscale('symlog')
		plt.xlim(0,5000)
		fig.axes[0].xaxis.set_major_formatter(ScalarFormatter())
		plt.xticks([0,1,2,3,5,10,20,40,100,300,1000,2000,5000], fontsize=10)
		return fig

