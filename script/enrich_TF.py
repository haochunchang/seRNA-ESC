# Enrichment of ENCODE TFs in state-specific seRNAs
#### Perform Fisher's exact test on state-specific seRNA and ENCODE TFs

from scipy import stats
import os, sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from misc import utils as misc_u

def calculate_tf_enrich(data, allse, se, state_id):
    '''
    data: ENCODE TF data
    allse: All seRNA location
    se: [state-specific seRNA location]
    return: {TF name: p-value in state_id}
    '''
    tf_enrich = {}
    xs = {}
    Ms = {}
    otherse = allse.drop(se.index, axis=0)
    for name, tfdata in data.groupby('symbol'): 
        x = 0
        M = 0
        tfdata_group = tfdata.groupby('chrom')
        for (other_chrom, other_data), (se_chrom, se_data) in zip(otherse.groupby('chrom'), se.groupby('chrom')):
            try:
                tf_chrom = tfdata_group.get_group(se_chrom)
                other_chrom = tfdata_group.get_group(other_chrom)
            except KeyError:
                continue
            contain_start = other_chrom['chromStart'].apply(lambda x: x >= other_data['start'])
            contain_end = other_chrom['chromEnd'].apply(lambda x: x <= other_data['end'])
            contain = contain_start & contain_end
            M += contain.sum().sum()
            
            contain_start = tf_chrom['chromStart'].apply(lambda x: x >= se_data['start'])
            contain_end = tf_chrom['chromEnd'].apply(lambda x: x <= se_data['end'])
            contain = contain_start & contain_end
            x += contain.sum().sum()
            
        xs[name] = x
        Ms[name] = M
        
    for name, x in xs.items():
        M = Ms[name]
        K = sum(xs.values()) - x
        N = sum(Ms.values()) - M
        table = [[x, M], [K, N]]
        try:
            o, p = stats.fisher_exact(table)
        except:
            print(table)
        tf_enrich[name] = -np.log10(p)
    return pd.Series(tf_enrich)

def plot_enrichment(tf_enrich, outfilename, alpha=0.05):
    
    sns.set(font_scale = 1.2)
    jp = sns.jointplot(x='state0 -log10(p-value)', y='state1 -log10(p-value)', data=tf_enrich, stat_func=None, height=12)
    
    for line in range(tf_enrich.shape[0]):
        #if tf_enrich['stateA -log10(p-value)'][line] > alpha or tf_enrich['stateB -log10(p-value)'][line] > alpha:
        jp.ax_joint.text(tf_enrich['state0 -log10(p-value)'][line]+0.01, tf_enrich['state1 -log10(p-value)'][line], 
                            tf_enrich.index[line], horizontalalignment='left', va='top', size='small', 
                             color='black', weight='semibold')

    jp.ax_joint.axhline(-np.log10(alpha))
    jp.ax_joint.axvline(-np.log10(alpha))
    jp.savefig(outfilename, dpi=600)
    print('Plot saved to {}'.format(outfilename))

def main(celltype, num_state=2):
    # Load in tf data
    path = os.path.join('..', 'data', 'raw', 'wgEncodeRegTfbsClusteredV3.bed')
    ENCODE_data = misc_u.load_encode(path)
    r = 500

    # Load in all seRNA location
    se = pd.read_csv(os.path.join('..', 'data', celltype, 'Super_Enhancers',
                                  'super-enhancer_avg_exp.tsv'), sep='\t', index_col=0)
    se['locus'] = se.index
    chrom = se['locus'].str.split(':').apply(lambda x: x[0])
    start = se['locus'].apply(lambda x: int(x.split(':')[1].split('..')[0]))
    end = se['locus'].apply(lambda x: int(x.split(':')[1].split('..')[1]))

    allse = pd.DataFrame()
    allse['chrom'] = chrom
    allse['start'] = start - r
    allse['end'] = end + r
    allse.index = se.index

    state_se_list = []
    for state_id in range(num_state):
        # Load in state-specific seRNA location
        state_se = pd.read_csv(os.path.join('..', 'data', celltype, 'selected_seRNA_state{}.csv'.format(state_id)))
        state_se_list.append(allse.loc[state_se['locus'], :])

        print(state_se.shape)

    # enrichment TF: [p-value for stateA, p-value for stateB]
    alpha = -np.log10(0.05)
    tf_enrich_inA = calculate_tf_enrich(ENCODE_data, allse, state_se_list[0], 0)
    tf_enrich_inB = calculate_tf_enrich(ENCODE_data, allse, state_se_list[1], 1)
    tf_enrich = pd.concat([tf_enrich_inA, tf_enrich_inB], axis=1)
    tf_enrich.columns = ['state0 -log10(p-value)', 'state1 -log10(p-value)']
    tf_enrich.to_csv(os.path.join('..', 'data', celltype, 'TF_enrichment.csv'), index=True)

    # Plot and save figure
    outfilename = os.path.join('..', 'graph', celltype, 'TF_enrichment.png')
    plot_enrichment(tf_enrich, outfilename, alpha)
    
if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(prog='TF Enrich on stage-specific seRNAs')
    parser.add_argument("--celltype", type=str, help="Which celltype", required=True)
    args = parser.parse_args()
    
    try:
        main(args.celltype)
    except Exception as e:
        print(e)
        print('Usage: python3.5 enrich_TF.py [celltype]')
