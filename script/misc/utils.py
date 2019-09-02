import pandas as pd
import numpy as np
import os

def ensure_dir(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

def average_replicates(mrna):
    '''
    Average time points
    ''' 
    #mrna = mrna.iloc[:, :-3]
    avg = mrna.T.groupby(level=[0]).mean()

    return avg.T

def load_encode(path):
    '''
    Read in wgEncodeRegTfbsClusteredV3.bed from path
    usage: 
        path = os.path.join('..', 'data', 'raw', 'wgEncodeRegTfbsClusteredV3.bed')
        ENCODE_data = load_encode(path)
    '''
    schema = ['chrom', 'chromStart', 'chromEnd', 'symbol', 'score', 'expCount', 'expNums', 'expScores']
    data = pd.read_csv(path, sep='\t', header=None, names=schema)

    return data

if __name__ == "__main__":
    celltype = 'HES3_GFP_ESC'
    se = pd.read_csv(os.path.join("../data", celltype, "super-enhancer_exp.tsv"), sep="\t", index_col=0, header=list(range(2)))
    se = average_replicates(se) 
    print(se)
    

