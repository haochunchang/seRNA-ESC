import pandas as pd
import sys, os
import argparse
from nmf import time_course_NMF
from nmf import utils
from misc import utils as misc_u
from rpy2.robjects import pandas2ri

pandas2ri.activate()

def main(celltype, avg_flag=False):
    
    # Load in seRNA and gene expression profile  
    gene = pd.read_csv(os.path.join("../data", celltype, "Super_Enhancers",
                                    "genes_"+celltype+".tsv"), sep='\t', index_col=0, header=list(range(2))) 
    
    gene_avg = misc_u.average_replicates(gene)
    gene_avg.to_csv(os.path.join('../data', celltype, "Super_Enhancers",
                                 'genes_'+celltype+'_avg.tsv'), sep='\t', index=True)
    
    # Average Time points if avg_flag == True
    if avg_flag:
        se = pd.read_csv(os.path.join("../data", celltype, "Super_Enhancers",
                                      "super-enhancer_exp.tsv"), sep="\t", index_col=0, header=list(range(2)))
        
        se = misc_u.average_replicates(se) 
        se.to_csv(os.path.join("../data", celltype, "Super_Enhancers",
                               "super-enhancer_avg_exp.tsv"), sep="\t", index=True) 
    else:
        se = pd.read_csv(os.path.join("../data", celltype, "Super_Enhancers",
                                      "super-enhancer_exp.tsv"), sep="\t", index_col=0, header=list(range(2)))
    
    # NMF, select best k
    k_lst =  [x for x in range(2, se.shape[1])] 
    graph_path = os.path.join('..', 'graph', celltype, 'seRNA')
    misc_u.ensure_dir(graph_path)
    
    sample2state, state2seRNA, best_k = time_course_NMF.model_select(se, k_lst, graph_path)
    misc_u.ensure_dir(os.path.join('../data', celltype, 'NMF'))
    sample2state.to_csv(os.path.join('../data', celltype, 'NMF', 'sample2state.csv'), index=True)
    state2seRNA.to_csv(os.path.join('../data', celltype, 'NMF', 'state2seRNA.csv'), index=True)
    print("H and W stored.")
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='NMF')
    parser.add_argument("-a", "--average", help="Flag indicating whether averaging time points.", action="store_true")
    parser.add_argument('celltype', action = 'store') 
    args = parser.parse_args()
    main(args.celltype, args.average)
