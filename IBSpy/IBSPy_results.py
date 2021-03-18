import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture

class IBSpyResults:
    # class variables go here
    # filter_counts = 2000

    def __init__(self, filename, window_size, filter_counts):
        #instance variables
        # self.filename = filename
        self.db = pd.read_csv(filename, delimiter='\t')
        # self.lengths =  pd.read_csv(chromosome_lengths, delimiter='\t', header=None)
        self.window_size = window_size
        self.filter_counts = filter_counts

    def count_by_windows(self):
        windSize = self.window_size
        db_DF = self.db[['seqname','start','end','variations']]
        # get the longest chromosome
        chrLen = db_DF['end'].max()
        # windows to iterate
        w_pos = 0
        db_byChr = pd.DataFrame()
        while w_pos <= chrLen:
            db_DF_ByWind = db_DF[(db_DF['start'] >= w_pos) & (db_DF['start'] < w_pos + windSize)]
            db_DF_ByWind = db_DF_ByWind.loc[:, ('seqname','start','variations')]
            db_DF_ByWind['window'] = w_pos + windSize
            w_pos += windSize
            db_byChr = db_byChr.append(db_DF_ByWind)

        # calculate statistics by chromosome, window, and variations
        db_ByGenom = db_byChr.groupby(['seqname','window'])['variations'].sum().reset_index()
        # mean
        db_byChr_mean = db_byChr.groupby(['seqname','window'])['variations'].mean().reset_index()
        db_ByGenom['v_mean'] = db_byChr_mean['variations'].astype(int)
        # median
        db_byChr_median = db_byChr.groupby(['seqname','window'])['variations'].median().reset_index()
        db_ByGenom['v_median'] = db_byChr_median['variations'].astype(int)
        # variance
        db_byChr_var = db_byChr.groupby(['seqname','window'])['variations'].var().reset_index()
        db_ByGenom['v_var'] = db_byChr_var['variations'].fillna(0).astype(int)
        # standar dev
        db_byChr_std = db_byChr.groupby(['seqname','window'])['variations'].std().reset_index()
        db_ByGenom['v_std'] = db_byChr_std['variations'].fillna(0).astype(int)

        self.tmp_table = db_ByGenom
        return db_ByGenom
    
    def normalize_data(self):
        by_windows_db = self.count_by_windows()
        filtered_by_windows_db = by_windows_db[by_windows_db['variations'] <= self.filter_counts]
        varDF = pd.DataFrame(filtered_by_windows_db[['seqname','window','variations']])
        varDF.reset_index(drop=True, inplace=True)
        varArray = np.array(varDF['variations'])
        log_varArray = np.log(varArray, where=(varArray != 0))
        log_varArray = np.nan_to_num(log_varArray, nan=0.0)
        log_varArray = log_varArray.reshape((len(log_varArray),1))
        return log_varArray, varDF

    def fit_gmm_model(self, n_components, covariance_type='full'):
        log_varArray, varDF = self.normalize_data() # this DF its filtered, above
        model = GaussianMixture(n_components=n_components, covariance_type=(covariance_type))
        model.fit(log_varArray)
        varGMM_predict = model.predict(log_varArray)
        varGMM_predict = varGMM_predict.reshape((len(varGMM_predict),1)) 
        varGMM_predict = pd.DataFrame({'v_gmm': varGMM_predict[:, 0]})     
        varDF['v_gmm'] = varGMM_predict['v_gmm']                                   
        gmmIBS_DF = varDF.groupby(by=['v_gmm'])
        # build binary
        key_group = []
        value_group = []
        for key, group in gmmIBS_DF:
            key_group.append(key)
            value_group.append(group['variations'].median())
        dic_grupby = dict(zip(key_group, value_group))
        IBS_group = min(dic_grupby, key=dic_grupby.get)
        varDF['v_gmm'] = np.where(varDF['v_gmm'] == IBS_group, 1, 0)
        
        self.tmp_varDF = varDF
        return varDF

    def stitch_haploBlocks(self, n_components,covariance_type,stitch_number):
    # change here, instead of chromosome, full genome
        stitch_db = self.fit_gmm_model(n_components, covariance_type='full')
        gmmDta = stitch_db['v_gmm'].copy()

        StitchVarNum = stitch_number
        for i in range(len(gmmDta)):
            if gmmDta[i] == 1:
                if i < (len(gmmDta) - StitchVarNum):
                    count = 0
                    for n in range(1, (StitchVarNum + 1)):
                        if gmmDta[i+n] == 0:
                            count += 1
                    if count == StitchVarNum:
                        continue
                    else:
                        gmmDta[i+1] = 1
        hapBlock = gmmDta
        hapBlock = pd.Series(hapBlock)
        stitch_db['vh_block'] = hapBlock.values
        # put back filtered data with haploblocks
        hapCntFile = self.count_by_windows()
        stitch_db = pd.merge(hapCntFile, stitch_db, left_on=['seqname','window'], right_on=['seqname','window'], how='left')
        stitch_db.loc[:,'v_gmm':'vh_block'] = np.where(stitch_db.loc[:, 'v_gmm':'vh_block'] == 1, 1, 0)
        stitch_db.rename(columns={'seqname_x':'seqname', 'variations_x':'variations'}, inplace=True)
        stitch_db = stitch_db.drop(['variations_y'], axis=1)
        return stitch_db
    
    def run_analysis(self):
        counts     = self.count_by_windows()
        # by_windows_db = self.count_by_windows()
        # normalised = self.normalize_data(by_windows_db)
#         self.fit_gmm_model()
#         self.stitch_haploBlocks()