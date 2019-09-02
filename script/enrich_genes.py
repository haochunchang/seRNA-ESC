from scipy import stats
import random, json, argparse, os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sns

def plot_enrich_heatmap(pval, es, coexpress, tfbmat, output_path):
    sns.set(font_scale = 1.2)
    
    pval = pval.loc[:, (pval < 0.05).any(axis=0)]
    pval = pval.loc[(pval < 0.05).any(axis=1), :]
    es = es.loc[pval.index, pval.columns]
    es = round(es, 2).astype(str)
    es = es.applymap(lambda x: "")
    es[pval < 0.05] = "*"
    print(es)
    # Check significant & large es associated genes also have their se bounded by TFs
    significant = es[pval < 0.05][es > 1].fillna(0)
    
    pval = round(pval.applymap(lambda x: -np.log10(x)), 2)
    sig = {}
    for seRNA, sig_series in significant.iterrows():
        sig[seRNA.replace('..', '-')] = [x for x in sig_series.index if tfbmat.loc[seRNA, x] > 0]
    
    # Plot
    print(pval.shape, es.shape)
    pval.index = pval.index.str.replace('..', '-', regex=False)
    cg = sns.clustermap(pval, method='average', metric='euclidean', row_cluster=False, col_cluster=False, 
                        linewidths=.05, square=False, 
                        annot=es.values, annot_kws={"size": 16}, fmt='', 
                        figsize=(12,8), center=-np.log10(0.05), 
                        cmap="Oranges", cbar_kws={'label': '-log10(p-value)'})
    
    for sig_seRNA, sig_tf in sig.items():
        for tf in sig_tf:
            if tf not in pval.columns: continue
            cg.ax_heatmap.add_patch(Rectangle((pval.columns.get_loc(tf), pval.index.get_loc(sig_seRNA)), 
                                              1, 1, fill=False, edgecolor='darkred', lw=2))
    #cg.ax_heatmap.set_xticks([])
    #cg.ax_heatmap.set_yticks([])
    cg.ax_heatmap.set_ylabel('State-specific seRNAs')
    cg.ax_heatmap.set_xlabel('Transcription Factors')
    
    cg.savefig(output_path, dpi=600, figsize=(12,8))
    #plt.show()
    #plt.close()


def calculate_gene_enrich(tfbmat, coexpress, num_genes, k=1000):
    '''
    Calculate TF enrichment score of each associated genes. 
    By Fisher's exact test:
    --------------------------------------------------------
    |     Count      | Associated genes | Randomized genes |
    --------------------------------------------------------
    | Bs of given TF |        X         |       M          | 
    --------------------------------------------------------
    | Bs of other TF |        K         |       N          |
    --------------------------------------------------------
    
    tfbmat: index=se+genes, columns=TFs, values=count of overlap bs
    coexpress: {seRNA: associated genes}
    
    return: DataFrame, index=seRNA, column=ENCODE TF, 
            value1 = p-value
            value2 = enrichment score
    '''
    tfbmat = tfbmat.astype(bool)
    
    random_candid = [i for i in range(num_genes)]
    num_se = tfbmat.shape[0] - num_genes
    gene_tfbmat = tfbmat.iloc[num_se:, :]
    se_tfbmat = tfbmat.iloc[:num_se, :]
    
    result_pval = pd.DataFrame(index=coexpress.keys(), columns=tfbmat.columns)
    result_es = pd.DataFrame(index=coexpress.keys(), columns=tfbmat.columns)
    # test each group
    for se, genes in coexpress.items():
        
        gene_bind = gene_tfbmat.loc[genes, :].sum(axis=0)
        if gene_bind.sum() == 0: print(gene_bind)
            
        random_gene_bind = []
        for i in range(k):
            random_index = random.sample(random_candid, len(genes))
            random_gene_bind.append(gene_tfbmat.iloc[random_index, :].sum(axis=0))
        
        # Ms: value=mean(overlapping_bs_count), columns=TF
        Ms = pd.DataFrame(random_gene_bind).mean(axis=0)
        for tf in tfbmat.columns:
            X = gene_bind[tf]
            K = gene_bind.sum() - X
            M = Ms[tf]
            N = Ms.sum() - M
            
            table = [[X, M], [K, N]]
            odds, pvalue = stats.fisher_exact(table, "greater")
            es = (X/(M+1e-10))/(K/N)
            if es > 1 and pvalue < 0.05 : 
                print(se, odds, pvalue, es)
                print(table)
            result_pval.loc[se, tf] = pvalue
            result_es.loc[se, tf] = es
            
    return result_pval, result_es

if __name__ == "__main__":
    
    from misc import utils as misc_u
    
    parser = argparse.ArgumentParser(prog='TF Enrich on associated genes')
    parser.add_argument("--celltype", type=str, help="Which celltype", required=True)
    parser.add_argument("--num_state", type=int, help="how many state from NMF")
    args = parser.parse_args()
    
    path = os.path.join('..', 'data', args.celltype)
    output_dir = os.path.join(path, 'TF_enrichment')
    
    tfbmat = pd.read_csv(os.path.join(path, 'tf_binding_matrix.csv'), index_col=0)
    mrna_pool = pd.read_csv(os.path.join(path, 'Super_Enhancers', 'genes_{}_avg.tsv'.format(args.celltype)), 
                            sep='\t', index_col=0)
    num_genes = mrna_pool.shape[0]
    
    # Load in co-express seRNA & mRNA location data
    for state_id in range(args.num_state):

        # From Mutual_Rank.ipynb
        with open(os.path.join(path, "coexpress_locus_state{}.json".format(state_id)), "r") as f:
            comrna = json.load(f)
        
        pval, es = calculate_gene_enrich(tfbmat, comrna, num_genes)
        
        misc_u.ensure_dir(output_dir)
        pval.to_csv(os.path.join(output_dir, "TFenrich_on_gene_state{}_pval.csv".format(state_id)), index=True)
        es.to_csv(os.path.join(output_dir, "TFenrich_on_gene_state{}_es.csv".format(state_id)), index=True)
        
        pval = pd.read_csv(os.path.join(output_dir, "TFenrich_on_gene_state{}_pval.csv".format(state_id)), index_col=0)
        es = pd.read_csv(os.path.join(output_dir, "TFenrich_on_gene_state{}_es.csv".format(state_id)), index_col=0)
        output_path = os.path.join('..', 'graph', args.celltype, 'TF_enrich_on_genes_state{}_new.jpg'.format(state_id))
        plot_enrich_heatmap(pval, es, comrna, tfbmat, output_path)
        
        
 