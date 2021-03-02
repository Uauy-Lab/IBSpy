import sys
import gzip
import argparse
import pandas as pd

def read_db():
    inDB = pd.read_csv(args.db_file, delimiter='\t')
    return inDB

def read_chrLen_f():
    in_chrLen = pd.read_csv(args.chrLen_f, delimiter='\t', header=None)
    return in_chrLen

#### file with chromosome lengths ####
def count_by_windows(args):
    chrLen_DF = read_chrLen_f()
    inDB = read_db()
    windSize = int(args.windSize)

    db_ByGenom = pd.DataFrame()
    for index, row in chrLen_DF.iterrows():
        ## get individual chromosome and size
        chrID = chrLen_DF[chrLen_DF[0] == row[0]]
        chrLen = chrID.iloc[0][1]
        ## access db data 
        db_DF = inDB[inDB['seqname'].str.contains(row[0])] # made a change to fit new files from IBSpy
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
    return db_ByGenom
    
def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-db', '--db_file',  required=True, help='Tab separated file with variations genetared by IBSpy output')
    parser.add_argument('-cl', '--chrLen_f',  required=True, help='Reference chromosome lenghts file')
    parser.add_argument('-w', '--windSize',  required=True, help='Windows size to count variations within')
    parser.add_argument('-o','--output', help='File with variations by windows only')

    args = parser.parse_args()
    return args

args = parse_arguments()
count_by_windows(args).to_csv(args.output, index=False,  sep='\t', compression="gzip")