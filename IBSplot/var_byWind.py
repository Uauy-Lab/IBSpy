import sys
import gzip
import argparse
import pandas as pd

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Count variations by windows sizes")

    db_p = parser.add_argument_group('Sample name parameters')
    db_p.add_argument('-db', '--db_file',  required=True, help='Tab separated file with variations genetared by IBSpy output')
    db_p.add_argument('-cl', '--chr_length_file',  required=True, help='Reference chromosome lenghts file')
    db_p.add_argument('-w', '--windSize',  required=True, help='Windows size to count variations within')

    out_files = parser.add_argument_group('Output files')
    out_files.add_argument('-o','--output', help='File with variations count by windows only')

    args = parser.parse_args()

    #### command line input arguments ####
    db_file = args.db_file
    chr_length_file = args.chr_length_file
    windSize = int(args.windSize)
    output = args.output

#### file with chromosome lengths ####
chrLen_DF = pd.read_csv(chr_length_file, delimiter='\t', header=None)
in_db_file = pd.read_csv(db_file, delimiter='\t')

db_ByGenom = pd.DataFrame()
for index, row in chrLen_DF.iterrows():
    ## get individual chrm and size
    chrID = chrLen_DF[chrLen_DF[0] == row[0]]
    chrLengh = chrID.iloc[0][1]
    ## access db data 
    db_DF = in_db_file[in_db_file['seqname'] == row[0]]
    db_DF = db_DF[['seqname','start','variations']]

    ## empty DF to merge per chromosome
    w_pos = 0
    db_byChr = pd.DataFrame()
    while w_pos <= chrLengh:
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

db_ByGenom.to_csv(output, index=False,  sep='\t')