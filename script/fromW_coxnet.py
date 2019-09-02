import pandas as pd
import numpy as np
import os, sys, json
import subprocess
import rpy2.robjects as robjects
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
import matplotlib.pyplot as plt
import seaborn as sns
from misc import utils
from network import bi_coxnet
from nmf import select
from preprocess import data_preprocess as pre

# Suppress R script warnings
import warnings
from rpy2.rinterface import RRuntimeWarning

warnings.filterwarnings("ignore", category=RRuntimeWarning)

def main(celltype, state_id=0):
    
    # Load in mRNA data
    mrna = pd.read_csv(os.path.join('../data', celltype, 'Super_Enhancers',
                                    'genes_'+celltype+'_avg.tsv'), sep='\t', index_col=0, header=0)

    # Load in seRNA data
    se = pd.read_csv(os.path.join('../data', celltype, 'Super_Enhancers',
                                  'super-enhancer_avg_exp.tsv'), sep='\t', index_col=0)
    se['locus'] = se.index
    seRNA_w = os.path.join('../data', celltype, 'NMF', 'state2seRNA.csv')
    selected_outpath = os.path.join('..', 'data', celltype,
                                    'selected_seRNA_state{}.csv'.format(state_id))

    # Select state-specific seRNA
    selected_serna = select.select_from_W(se, w_path=seRNA_w, outfilepath=selected_outpath, 
                                         topn='all', state_id=state_id)

    print('Selection of state {}-specfic seRNA done'.format(state_id))
    
    # Construct seRNA to mRNA bipartite co-expression network
    net = bi_coxnet.BiCoxNet(mrna, se, 'fromW-'+celltype, state_id) 
    net.calculate_corr('pearson', p_value=0, save=False)
    #net.corr.to_csv('net_corr_mat{}.csv'.format(state_id))

    # Select mRNA based on mutual rank
    filtered = filter_mRNA(net, celltype, state_id)
    ss_mrna = net.corr.loc[filtered.loc[:, selected_serna['locus']].astype(bool).sum(axis=1) > 0, selected_serna['locus']]
    print('State-specific mRNA: {}'.format(ss_mrna.shape))
    ss_mrna.to_csv(os.path.join('..', 'data', celltype,
                                'state{}-specific_mRNA.csv'.format(state_id)), index=True)

    # Convert mRNA loci to gene symbol
    # And get gene symbols and locus for every given seRNA
    save_sym_path = os.path.join('..', 'data', celltype,
                                 'coexpress_mrna_state{}.json'.format(state_id))
    save_loci_path = os.path.join('..', 'data', celltype,
                                  'coexpress_locus_state{}.json'.format(state_id))
    interaction = filtered.loc[:, selected_serna['locus']].astype(bool)
    get_coexpress_symbol_and_loci(interaction, celltype, save_sym_path, save_loci_path)

    
def filter_mRNA(net, celltype, state_id, threshold=5):
    '''
    Filter mRNA based on mutual rank 
    @param: net, BiCoxNet class object
    @param: celltype, state_id, for saveing file.
    '''
    seRNA_rank = net.corr.rank(axis=0, ascending=False)
    mRNA_rank = net.corr.rank(axis=1, ascending=False)
    # Compress rank to equal range
    old_min = seRNA_rank.min().min()
    old_max = seRNA_rank.max().max()
    new_min = mRNA_rank.min().min()
    new_max = mRNA_rank.max().max()

    seRNA_rank_compress = (((new_max-1)/(old_max-1)) * (seRNA_rank-1)) + 1 
    
    mutual_rank = np.sqrt(mRNA_rank.values * seRNA_rank_compress.values)
    mutual_rank = pd.DataFrame(mutual_rank)
    mutual_rank.columns = net.serna_locus
    mutual_rank.index = net.mrna_locus

    filtered = mutual_rank[mutual_rank <= threshold].fillna(0)
    return filtered


def get_coexpress_symbol_and_loci(filtered, celltype, save_sym_path=None, save_loci_path=None):
    '''
    Convert mRNA loci to gene symbol
    And get gene symbols and locus for every given seRNA
    Save as json in save_path
    '''
    symbol = pd.read_csv(os.path.join('..', 'data', celltype, "Super_Enhancers",
                                      'symbol_{}.tsv'.format(celltype)), 
                        sep='\t', header=None, index_col=0)
    ss_mrna = filtered
    coexpress_mrna = {}
    coexpress_locus = {}
    for ser in ss_mrna:
        mrna_lst = []
        for i in range(ss_mrna[ser].shape[0]):
            if not ss_mrna[ser].iloc[i]:
                continue
            else:
                mrna_lst.append(ss_mrna.index[i])
        sym = symbol.loc[mrna_lst].dropna().values.tolist()

        coexpress_mrna[ser] = list(set([i[0] if len(i) == 1 else [j for j in i] for i in sym]))
        coexpress_locus[ser] = list(set([i for i in mrna_lst]))
    if save_loci_path is not None:
        with open(save_loci_path, 'w') as f:
            json.dump(coexpress_locus, f)
    
    if save_sym_path is not None:
        with open(save_sym_path, 'w') as f:
            json.dump(coexpress_mrna, f)

    return coexpress_mrna, coexpress_locus


if __name__ == "__main__":
    try:
        nstates = int(sys.argv[2])
        for state_id in range(nstates):
            main(sys.argv[1], state_id)
    except IndexError:
        print('Usage: python3 fromW_coxnet.py [celltype] [# of states from NMF]')
