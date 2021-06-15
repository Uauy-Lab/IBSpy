import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture

class IBSpyResults:
    # class variables go here
    def __init__(self, filename, window_size, filter_counts):

        self.db = pd.read_csv(filename, delimiter='\t')
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

        # calculate statistics by chromosome, windows, and variations
        out_db = db_byChr.groupby(['seqname','window'])['variations'].sum().reset_index()
        out_db['mean'] = db_byChr.groupby(['seqname','window'])['variations'].mean().values.round(2)
        out_db['median'] = db_byChr.groupby(['seqname','window'])['variations'].median().values.round(2)
        out_db['variance'] = db_byChr.groupby(['seqname','window'])['variations'].var().values.round(2)
        out_db['std'] = db_byChr.groupby(['seqname','window'])['variations'].std().values.round(2)
        return out_db
    
    def transform_counts_to_log(self, counts):
        by_windows_db = counts 
        #by_windows_db = self.count_by_windows()
        if self.filter_counts is not None:
            applied_filter = by_windows_db[by_windows_db['variations'] <= self.filter_counts]
        else:
            applied_filter = by_windows_db
        varDF = pd.DataFrame(applied_filter[['seqname','window','variations']])
        varDF.reset_index(drop=True, inplace=True)
        varArray = varDF['variations'].to_numpy()
        log_varArray = np.log(varArray, where=varArray>0)
        log_varArray = log_varArray.reshape((len(log_varArray),1))
        return log_varArray, varDF

    def build_gmm_model(self, log_varArray, varDF, n_components, covariance_type):
        #log_varArray, varDF = self.transform_counts_to_log()
        model = GaussianMixture(n_components=n_components, covariance_type=covariance_type)
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

    def stitch_gmm_haplotypes(self, model, stitch_number):
        #model = self.build_gmm_model(n_components, covariance_type)
        gmmDta = model['v_gmm'].copy()

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
        model['vh_block'] = hapBlock.values
        # put back filtered data with haploblocks
        hapCntFile = self.count_by_windows()
        model = pd.merge(hapCntFile, model, left_on=['seqname','window'], right_on=['seqname','window'], how='left')
        model.loc[:,'v_gmm':'vh_block'] = np.where(model.loc[:, 'v_gmm':'vh_block'] == 1, 1, 0)
        model.rename(columns={'seqname_x':'seqname', 'variations_x':'variations'}, inplace=True)
        model = model.drop(['variations_y'], axis=1)
        return model
    
    def run_analysis(self, n_components, covariance_type, stitch_number):
        counts     = self.count_by_windows()
        log_test, pd = self.transform_counts_to_log(counts)

        model = self.build_gmm_model(log_test, pd, n_components, covariance_type) 
        hap_pd = self.stitch_gmm_haplotypes(model, stitch_number)
        return hap_pd
        # normalised = self.transform_counts_to_log(by_windows_db)
#         self.build_gmm_model()
#         self.stitch_gmm_haplotypes()