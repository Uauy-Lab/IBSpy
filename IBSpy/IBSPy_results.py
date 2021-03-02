import pandas as pd

class IBSpyResults:
    def __init__(self, filename, chromosome_lengths, window_size):
        self.filename = filename
        self.db = pd.read_csv(filename, delimiter='\t')
        self.lengths =  pd.read_csv(chromosome_lengths, delimiter='\t', header=None)
        self.window_size = window_size


#### file with chromosome lengths ####
    def count_by_windows(self):
        windSize = int(self.window_size)
        db_ByGenom = pd.DataFrame()
        for index, row in self.lengths.iterrows():
            ## get individual chromosome and size
            chrID = self.lengths[self.lengths[0] == row[0]]
            chrLen = chrID.iloc[0][1]
            ## access db data 
            db_DF = self.db[self.db['seqname'].str.contains(row[0])] # made a change to fit new files from IBSpy
            db_DF = db_DF[['seqname','start','variations']]

            ## empty DF to merge by chromosome
            w_pos = 0
            db_byChr = pd.DataFrame()
            while w_pos <= chrLen:
                db_DF_ByWind = db_DF[(db_DF['start'] > w_pos) & (db_DF['start'] <= w_pos + windSize)]
                db_DF_ByWind = db_DF_ByWind.loc[:, ('seqname','start','variations')]
                db_DF_ByWind['w_num'] = w_pos + windSize
                w_pos += windSize
                db_byChr = db_byChr.append(db_DF_ByWind)

            db_byChr = db_byChr.groupby(['w_num']).sum().reset_index()
            db_byChr = db_byChr[['w_num', 'variations']]
            db_byChr['seqname'] = row[0]
            db_byChr = db_byChr[['seqname', 'w_num', 'variations']]
            db_byChr.rename(columns={'w_num':'start'}, inplace=True)
            db_ByGenom = db_ByGenom.append(db_byChr)
        self.tmp_table = db_ByGenom
        return db_ByGenom

    def transform_gmm(self):
        return None


    def run_analysis(self):
        self.count_by_windows()
        self.transform_gmm()