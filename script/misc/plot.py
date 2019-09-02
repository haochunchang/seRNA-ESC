import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys, os

def plot_expression_by_time():
    '''
    Plot boxplot of expression levels.
    X: Time point
    Y: log2 expression level of each locus.
    '''
    #super_exp = pd.read_csv("../DATA/super_enhancer_log2_exp.tsv", sep = "\t", index_col = 0, header = list(range(3)))
    super_exp = pd.read_csv("../DATA/super_enhancer_log2_normalized.tsv", sep = "\t", index_col = 0)
    print(super_exp)
    labels = [day[0] for day in super_exp.columns]
    
    super_exp = np.array(super_exp)
    plt.boxplot(super_exp, labels = labels)
    locs, labels = plt.xticks()
    plt.setp(labels, rotation = 90)
    plt.title("Super-enhancer")
    plt.show()
    plt.close()


def plot_corr(celltype, state_id, flag='reconstructed'):
    '''
    Input: 
        string celltype: e.g. HES3_GFP_ESC
        int state_id: e.g. 0, 1, ...
    '''

    # Load in re_mrna and super_enhancer
    sample2state = pd.read_csv(os.path.join('../../data', celltype, 'sample2state.csv'))
    state2seRNA = pd.read_csv(os.path.join('../../data', celltype, 'state2seRNA.csv'))
    n_serna = state2seRNA.shape[0]

    if flag == 'fitted':
        # mRNA
        X = np.loadtxt(os.path.join('../data', celltype, 'mRNA_fitted.txt'))
        

        # seRNA
        seRNA = state2seRNA.iloc[:, state_id].values.reshape(n_serna, 1)
        sample = sample2state.iloc[state_id, :].values.reshape(1, 13)
        re_serna = pd.DataFrame(np.dot(seRNA, sample))
        
        # Combine re_mrna and re_serna
        merged = pd.concat([re_mrna, re_serna], axis=0, join='inner')
    
    elif flag == 'reconstructed':
        exp_profile_path = os.path.join('../../data', celltype, 'Reconstructed_state{}_exp.csv'.format(state_id))
        merged = pd.read_csv(exp_profile_path, sep='\t')         
        merged = merged.drop(['locus'], axis=1) 
   
    else:
        print('Usage: flag= reconstructed(default) / fitted') 
    #merged = merged.set_index(locus_list[0][:1000])
    cor = merged.iloc[:2000, :].T.corr('pearson')
    cor = cor.sort_values(cor.columns[1], axis=0)
    fig, ax = plt.subplots(figsize=(20, 10))
    corr = ax.matshow(cor)
    fig.colorbar(corr)
    plt.show()

if __name__ == "__main__":
    celltype = sys.argv[1]
    state_id = int(sys.argv[2])

    plot_corr(celltype, state_id, flag=sys.argv[3])
