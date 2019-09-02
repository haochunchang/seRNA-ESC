"""
Automate functional enrichment analysis through R

Produce functional homogeneity score by GO term semantic similarity
Call R script from fea/
"""
import json, os, sys, argparse
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from misc import utils as misc_u

# Suppress R script warnings
import warnings
from rpy2.rinterface import RRuntimeWarning

warnings.filterwarnings("ignore", category=RRuntimeWarning)
pandas2ri.activate()


def auto_go(coexpress_mrna, celltype, state_id):
    '''
    Wrapper of FEA.R FEA
    '''
    for se, v in coexpress_mrna.items():
        se = se.replace('..', '_').replace(':', '_')
        r_sess = robjects.r
        r_sess.source('./fea/FEA.R')
        r_sess.TopGO(v, celltype, state_id, se)
        print("TopGO of {} done.".format(se))

        
def auto_gosim(comrna, pool_path):
    
    from rpy2.robjects.conversion import localconverter
    
    # drop duplicate gene symbols
    pool = pd.read_csv(pool_path, sep='\t', index_col=0, header=None, names=['symbol'])
    pool = pool.drop_duplicates()
    
    with localconverter(robjects.default_converter + pandas2ri.converter):
        sym_pool = robjects.conversion.py2rpy(pool)   
    #sym_pool = pandas2ri.py2ri(pool)

    r_sess = robjects.r
    r_sess.source('./fea/GOSIM.R')
    print("Getting GOSemSim data\n")
    hsGO = r_sess.godata('org.Hs.eg.db', keytype="SYMBOL", ont="BP")
    print("Clean symbol pool to only hgnc symbols\n")
    sym_pool = r_sess.only_hgnc(sym_pool)
    
    fcon = []
    p_val = []
    for se, mrna in comrna.items():
        result = r_sess.gosim(mrna, hsGO, sym_pool) 
        consist = result[0]
        serna_p = result[1]
        fcon.append(consist)
        p_val.append(serna_p)
        print("{}: homogeneity = {}, p-value: {}".format(se, consist, serna_p))
        
    pv = pd.DataFrame()
    pv['se'] = comrna.keys()
    pv['consistency'] = fcon
    pv['p_value'] = p_val
    return pv


def auto_GOSIM_main(celltype, num_state=2):
    for state_id in range(num_state):
        coexpress_path = os.path.join('..', 'data', celltype, 'coexpress_mrna_state{}.json'.format(state_id))
        with open(coexpress_path, "r") as f:
            comrna = json.load(f)
            
        pool_path = os.path.join('..', 'data', celltype, 'Super_Enhancers', 'symbol_{}.tsv'.format(celltype))
        pv = auto_gosim(comrna, pool_path)
        pv.to_csv(os.path.join('..', 'data', celltype, 'state{}_functional_homogeneity.csv'.format(state_id)))

        
def auto_FEA_main(celltype, num_state=2):
    for state_id in range(num_state):
        coexpress_path = os.path.join('..', 'data', celltype,
                                      'coexpress_mrna_state{}.json'.format(state_id))
        with open(coexpress_path, "r") as f:
            comrna = json.load(f)
            
        auto_go(comrna, celltype, state_id)
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--type", required=True, type=str, help="Which type of automatic Gene set analysis? (FEA / GOSIM)")
    parser.add_argument("--celltype", required=True, type=str, help="Which celltype?")
    parser.add_argument("--num_state", required=True, type=int, help="How many hidden state? (refer to NMF.py)")
    args = parser.parse_args()
    if args.type == 'GOSIM':
        auto_GOSIM_main(args.celltype, args.num_state)
    elif args.type == "FEA":
        auto_FEA_main(args.celltype, args.num_state)