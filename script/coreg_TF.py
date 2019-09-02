"""
Compute co-regulation score of co-expressed seRNA-mRNA modules
"""
import pickle, os, sys
import numpy as np
import pandas as pd
from misc import utils as misc_u
import pyximport; pyximport.install(pyimport=True)
from coreglib import *

def main(celltype, num_state=2):

    # Load in co-express seRNA & mRNA location data
    comrna_lst = []
    for state_id in range(num_state):
        path = os.path.join('..', 'data', celltype, "coexpress_locus_state{}.pkl".format(state_id))

        # From Mutual_Rank.ipynb
        with open(path, "rb") as f:
            comrna = pickle.load(f)
        comrna_lst.append(comrna)

    #=========================================
    # Load in tf data
    path = os.path.join('..', 'data', 'raw', 'wgEncodeRegTfbsClusteredV3.bed')
    ENCODE_data = misc_u.load_encode(path)

    #==========================================
    # random select the same size of mRNAs
    # Load in mRNA data
    gene = pd.read_csv(os.path.join('..', 'data', celltype, 'Super_Enhancers',
                                    'genes_{}_avg.tsv'.format(celltype)), 
                       sep='\t', index_col=0)
    
    coreg = [0 for i in range(num_state)]
    for i in range(num_state):
        path = os.path.join('..', 'data', celltype, 'state{}_tf_bmat'.format(i))
        coreg[i] = compute_coreg_and_test(comrna_lst[i], gene, ENCODE_data, path, r=500, k=100)
        
        print('Start writing coreg into pkl file...')
        pkl_path = os.path.join('..', 'data', celltype, 'coreg{}.pkl'.format(i))
        with open(pkl_path, 'wb') as f:
            pickle.dump(coreg[i], f)
        
            
if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(prog='TF Enrich on associated genes')
    parser.add_argument("--celltype", type=str, help="Which celltype", required=True)
    parser.add_argument("--num_state", type=int, help="how many state from NMF")
    args = parser.parse_args()
    try:
        main(args.celltype, args.num_state)
    except Exception as e:
        print(e)
        print('Usage: python3.5 coreg_TF.py [celltype] [num_state]')
