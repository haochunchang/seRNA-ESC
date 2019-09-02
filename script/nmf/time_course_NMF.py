import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from preprocess import data_preprocess as pre
import sys, os.path
import joblib
from sklearn.decomposition import NMF
from sklearn.metrics import silhouette_score
import fastcluster

sys.setrecursionlimit(100000)

def nmf(data, k, random_seed = 0):
    '''Use NMF to factorize super-enhancer expression matrix
    
    data: Expression matrix.
    k: Factorize into k states.
    random_seed: Initialize matrix configuration.
    
    '''
    model = NMF(n_components=k, init="nndsvd", random_state=random_seed)
    model.fit(data.values)
    H = pd.DataFrame(model.components_, columns=data.columns.get_level_values(0))
    W = pd.DataFrame(model.transform(data.values), index=data.index)
    return (H, W)

def show_nmf(mat, k, outfilepath, Wflag = False): 
    '''Save NMF clustermap to outfilepath
    
    mat: NMF matrix component
    Wflag: flag of indicating H is (locus x state) when it is True, else if H is (sample x state), it is False.
    
    '''
    if not Wflag:
        H_colors = serial_palette(mat.columns)
    elif Wflag:
        W_colors = serial_palette(mat.index)
    cmap = sns.light_palette((260, 100, 60), input="husl", as_cmap=True)

    with sns.plotting_context("talk"):
        #sns.set(font_scale = 3)
        if not Wflag:
            cm = sns.clustermap(mat, row_cluster=False, col_cluster=False,\
                                cmap=cmap, figsize=(24, 16))
            plt.setp(cm.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
            fig = cm.fig 
            fig.savefig(os.path.join(outfilepath, "nmf_"+str(k)+"_H.png"), dpi=800) 
            plt.close()
        elif Wflag:
            cm = sns.clustermap(mat.T, row_cluster=False, col_cluster=True, cmap=cmap, figsize=(15, 12))
            plt.setp(cm.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
            plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels())
            fig = cm.fig
            fig.savefig(os.path.join(outfilepath, "nmf_"+str(k)+"_W.png"), dpi=800) 
            plt.close()        
    

def serial_palette(objs):
    from functools import reduce
    '''Map objs to color series, then return the mapping relationship'''
    
    odered_set = reduce(lambda x, y: x if y in x else [*x, y], objs, [])
    pal = sns.diverging_palette(240, 10, n=len(odered_set), center="dark")
    lut = dict(zip(map(str, odered_set), pal))
    return pd.Series(lut)


def model_select(data, k_lst, outfilepath): 
    '''Return factorplot
    
    x-axis: a list of given k 
    y-axis: different sihouette scores
    
    '''
    W_err = 0
    Hscores = pd.DataFrame(columns = ['k', 'score'])
    Wscores = pd.DataFrame(columns = ['k', 'score'])
    for k in k_lst:
        for random_seed in range(200):
            H, W = nmf(data, k, random_seed = random_seed)
            score = silhouette_score(H.T, sample2state(H), random_state = random_seed)
            Hscores = Hscores.append(pd.DataFrame([(k, score)], columns = ['k', 'score']))
            try:
                score = silhouette_score(W, locus2state(W), random_state = random_seed)
            except MemoryError:
                W_err = 1
                print("Memory Error: k = {}".format(k))
                continue
            Wscores = Wscores.append(pd.DataFrame([(k, score)], columns = ['k', 'score']))
    # Save matrix which has largest silhouette score
    best_k = int(Hscores.loc[Hscores['score'] == Hscores['score'].max()]['k']) 
    Hscores.to_csv(os.path.join(outfilepath, 'silhouette_score_H.csv'))
    # Plot silhouette score plot and Factorized matrix.
    sns.set(font_scale=1.5) 
    fig = sns.catplot(x='k', y='score', data=Hscores, ci='sd',
                      kind='point', height=7, aspect=1.5) 
    fig.savefig(os.path.join(outfilepath, "silhouette_score_H.png"), dpi=800) 
    plt.close()        
        
    best_H, best_W = nmf(data, best_k)
    show_nmf(best_H, best_k, outfilepath)
    
    if not W_err:
        sns.set(font_scale=1.5) 
        fig = sns.catplot(x='k', y='score', data=Wscores, 
                          kind='point', height=7, aspect=1.5) 
        show_nmf(best_W, best_k, outfilepath, Wflag = True)
        

    return best_H, best_W, best_k
 
def sample2state(mat):
    '''Determine which sample belongs to which state by selecting the maximum value in its column'''
    return np.array([col.idxmax() for _, col in mat.iteritems()])

def locus2state(mat):
    '''Determine which locus belongs to which state by selecting the maximum value in its row'''
    return np.array([row.idxmax() for _, row in mat.iterrows()])

if __name__ == "__main__":
    pass
    
   
