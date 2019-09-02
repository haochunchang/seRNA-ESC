# Brwose many GO results
import pandas as pd
import glob, os, pickle
#import auto_GO as auto

def browseGO(celltype, state_id, s):
    
    s = {0: 'A', 1: 'B'}
    paths = glob.glob(os.path.join('..', 'data', celltype, 'FEA', 'stage{}'.format(state_id), 
                             '*', 'top20GO_bp.csv'))
    
    # Rank the seRNA with the most # of mRNA interaction
    #coexpress = auto.get_coexpress_mrna(celltype, state_id)
    with open('coexpress_mrna_state{}.pkl'.format(s[state_id]), 'rb') as f:
        coexpress = pickle.load(f)

    mrna_count_rank = []
    
    for se, mrna in coexpress.items():
        mrna_count_rank.append((se, len(mrna)))
    mrna_count_rank = sorted(mrna_count_rank, key=lambda x: x[1], reverse=True)
    serna_rank = [i[0] for i in mrna_count_rank]
    print(len(mrna_count_rank))

    go = {}
    for p in paths:
        locus = p.split('/')[-2]
        gobp = pd.read_csv(p).iloc[:,1]
        go[locus] = gobp
    
    GOtable = pd.DataFrame()
    for se, term in go.items():
        single = pd.Series(term, name=se)
        GOtable = GOtable.append(single, ignore_index=True)
    GOtable.index = go.keys()
    print(GOtable)
    GOtable.to_csv('./firstGO_state{}.csv'.format(state_id), index=True)

    
def getGO_by_fcon(celltype, state_id, alpha=0.05):

    copath = os.path.join('..', 'data', celltype, 'coexpress_mrna_state{}.pkl'.format(state_id))
    with open(copath, 'rb') as f:
        coexpress = pickle.load(f)

    fcon = pd.read_csv(os.path.join('..', 'data', celltype, 
                        'state{}_functional_consistency.csv'.format(state_id)), usecols=['se', 'consistency'], index_col=0)
    fcon = fcon['consistency'].apply(lambda x: float(x.split(' ')[-1].strip()))
    fcon = fcon[fcon < alpha]
    homo = [os.path.join('..', 'data', celltype, 'FEA', 'stage{}'.format(state_id), 
                             '{}'.format(i), 'top20GO_bp.csv') for i in list(fcon.index)]
    # get GO terms
    l = []
    g = []
    top = 3
    homotable = pd.DataFrame()
    for h in homo:
        locus = h.split('/')[-2]
        go = pd.read_csv(h).iloc[:top, 1]
        l.append(locus)
        g.append(go.values.tolist())

    homotable['se'] = l
    for i in range(top):
        homotable['go_{}'.format(i)] = [x[i] for x in g]
    homotable.to_csv(os.path.join('..', 'data', celltype, 
                        'homogeneous_go_term{}.csv'.format(state_id)), sep=',', index=False)

if __name__ == "__main__":
    
    celltype = "HES3_GFP_ESC"
    state_id = 1
    getGO_by_fcon(celltype, state_id)
