import pandas as pd
import numpy as np
import os, sys, pickle

## Select top n mRNA or seRNA from NMF-factorized W matrix
# @param exp_profile: exp profile of mRNA or seRNA
# @param w_path: filepath of state x (mRNA or seRNA) .csv file
# @param outfilepath: output file 
# @param topn: top n mRNA or seRNA, 'all' means all positive weight
# @param state_id: int of which state
# @param select: indicate selecting target: mRNA or seRNA
# return: selected top n target
def select_from_W(exp_profile, outfilepath, w_path=None, w=None, topn='all', state_id=0):
    
    # scale to mean = 1, and subtract max(other states)
    if w_path != None:
        state2target = pd.read_csv(w_path)
        state2target = state2target.loc[:, state2target.dtypes == np.float64].values
    elif w != None:
        state2target = w
    else:
        raise ValueError
    state2target = state2target / state2target.mean(axis=0)
    a = np.ma.array(state2target, mask=False)
    a.mask[:, state_id] = True
    state2target = state2target - a.max(axis=1).reshape((a.shape[0], 1))
    
    # selection based on scaled + compared W
    if topn != 'all':
        portion = state2target.shape[0] // topn
        top_per_idx = state2target[:, state_id].argsort()[::-1][:portion]
        selected = pd.DataFrame(exp_profile.iloc[top_per_idx, :])
        selected['locus'] = exp_profile.index[top_per_idx]
    else:
        # Select >= 2*std of difference distribution of eRNA as 'important' to state_id
        tmp = state2target[:, state_id]
        std = state2target[:, state_id].std()
        selected = pd.DataFrame(exp_profile[tmp>=2*std])
        selected['locus'] = selected.index
        
    selected.to_csv(outfilepath, index=False)

    return selected
