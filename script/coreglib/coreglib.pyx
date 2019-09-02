## Co-regulatory score of seRNA & associated genes
"""
## Co-regulatory score of seRNA & associated genes
### Main: coreg_TF.py

* Binding to super-enhancer: binding sites within +- 'r' bp of super-enhancer
* Binding to genes: binding sites within +- 'r' bp of START site of genes
* Compute (Jaccard + Overlap)/2 between genes and se

"""
import pickle, os, sys
import numpy as np
import pandas as pd

def single_coreg_score(x, y):
    '''
    x, y: 1d bool array
    jaccard: intersection(A, B) / Union(A, B)
    overlap_coeff: intersection(A, B) / min(A, B)
    '''
    if x.sum() == 0 and y.sum() == 0:
        return np.NaN
    inter = (x & y).sum()
    min_elem = min(x.sum(), y.sum())
    J = inter / (x.sum() + y.sum() - inter)
    
    if x.sum() == 0 or y.sum() == 0:
        O = 0
    else:
        O = inter / min_elem
    return (J+O)/2


def populate_tf_binding_matrix(coexpress, tf, bmat):
    '''
    coexpress: DataFrame: index=seRNA, columns=[chrom, start, end]
    tf: DataFrame: ENCODE_data
    bmat: DataFrame: empty binding matrix
    '''
    for sym_name, tfdata in tf.groupby('symbol'):
        tfdata_group = tfdata.groupby('chrom')
        for co_chrom, co_data in coexpress.groupby('chrom'):
            if co_chrom in tfdata_group.groups.keys():
                tf_chrom_data = tfdata_group.get_group(co_chrom)
            else:
                continue
            contain_start = co_data['start'].apply(lambda x: x <= tf_chrom_data['chromStart'])
            contain_end = co_data['end'].apply(lambda x: x >= tf_chrom_data['chromEnd'])
            is_binding = (contain_end & contain_start).any(axis=1)
            bmat.loc[is_binding.index, sym_name] = is_binding
    
    return bmat.astype(bool)


def build_tf_binding_matrix(serna_pool, mrna_pool, ENCODE_data, path, r=500):
    '''
    For genes, use TSS +- r as range
    For serna, use Loci +- r as range
    
    mrna_pool: DataFrame, all genes loci
    serna_pool: DataFrame, all serna loci
    ENCODE_data: DataFrame, TF binding site
    r: extended base pairs
    Return: DataFrame, index:all loci, columns: TF symbols, value: bool(is_binding)
    '''
    # Preprocess mrna_pool, for genes, use TSS +- r as range 
    chrom = [x.split(':')[0] for x in mrna_pool.index]
    start = [int(x.split(':')[-1].split('..')[0])-r for x in mrna_pool.index]
    end = [x + r for x in start]
    mrna_df = pd.DataFrame({'chrom':chrom, 'start':start, 'end':end}, index=mrna_pool.index)
    
    # Preprocess serna_pool, use range +- r as range 
    chrom = [x.split(':')[0] for x in serna_pool.index]
    start = [int(x.split(':')[-1].split('..')[0])-r for x in serna_pool.index]
    end = [int(x.split(':')[-1].split('..')[-1])+r for x in serna_pool.index]
    serna_df = pd.DataFrame({'chrom':chrom, 'start':start, 'end':end}, index=serna_pool.index)
    
    all_pool = serna_df.append(mrna_df)
    
    # Preprocess ENCODE_data and build binding matrix
    colnames = ENCODE_data.symbol.unique()
    rownames = list(all_pool.index)
    bools = np.zeros((len(rownames), len(colnames)))
    empty_bmat = pd.DataFrame(bools, index=rownames, columns=colnames)
    
    tf_binding_mat = populate_tf_binding_matrix(all_pool, ENCODE_data, empty_bmat)
    tf_binding_mat.to_csv(os.path.join(path, 'tf_binding_matrix.csv'), index=True)
    print('Binding matrix of saved')
    return tf_binding_mat


def compute_coreg_and_test(comrna, mrna_pool, ENCODE_data, path, r=500, k=100):
    '''
    comrna: dict: {seRNA locus: [co-expressed gene loci]}
    mrna_pool: DataFrame, all genes locus
    ENCODE_data: DataFrame, TF binding site
    r: extended base pairs
    Return: dict: {seRNA locus: [co-regulatory score, p-value]}
    '''
    # Preprocess mrna_pool 
    chrom = [x.split(':')[0] for x in mrna_pool.index]
    start = [int(x.split(':')[-1].split('..')[0])-r for x in mrna_pool.index]
    end = [x + r for x in start]
    pool = pd.DataFrame({'chrom':chrom, 'start':start, 'end':end}, index=mrna_pool.index)
    
    # Preprocess ENCODE_data and build empty matrix
    colnames = ENCODE_data.symbol.unique()
    rownames = list(pool.index) + list(comrna.keys())
    bools = np.zeros((len(rownames), len(colnames)))
    empty_bmat = pd.DataFrame(bools, index=rownames, columns=colnames, dtype=bool)
    
    coreg = {}
    for se, coexpress in comrna.items():
        # get coexpress location: for genes, use TSS +- r as range
        loc = pool.loc[coexpress,:]
        seloc = pd.Series({'chrom': se.split(':')[0], 
                           'start': int(se.split(':')[-1].split('..')[0]) - r, 
                           'end': int(se.split(':')[-1].split('..')[-1]) + r}, name=se)
        loc = loc.append(seloc)

        subset_empty_bmat = empty_bmat.loc[loc.index]
        tf_binding_mat = populate_tf_binding_matrix(loc, ENCODE_data, subset_empty_bmat)
        tf_binding_mat.to_csv(os.path.join(path, 'tf_bmat_{}.csv'.format(se)), index=True)
        print('Binding matrix of {} saved'.format(se))
        
        # calculate (J + O)/2 for each mRNA to seRNA and use mean(ri) as score of co-regulatory
        try:
            serna = tf_binding_mat.loc[se]
        except KeyError:
            coreg[se] = 0
            continue
        mrna = tf_binding_mat.drop(se, axis=0)

        ri = mrna.apply(lambda x: single_coreg_score(x, y=serna), axis=1)
        score = ri.mean(skipna=True)
        print('Coregulatory score of {}: {}'.format(se, score))
        
        #===================================
        # Testing by random sample genes
        
        count_of_bigger = 0
        for j in range(k):
            random_gene = pool.sample(n=loc.shape[0], replace=False)
            random_gene = random_gene.append(seloc)
            
            subset_empty_bmat = empty_bmat.loc[random_gene.index]
            random_tf_binding_mat = populate_tf_binding_matrix(random_gene, ENCODE_data, subset_empty_bmat)

            try:
                serna = random_tf_binding_mat.loc[se]
            except KeyError:
                random_score = 0
                if random_score > score: count_of_bigger += 1
                continue
            mrna = tf_binding_mat.drop(se, axis=0)
            random_score = mrna.apply(lambda x: single_coreg_score(x, y=serna), axis=1).mean(skipna=True)
            if random_score > score: count_of_bigger += 1
        
        p_val = count_of_bigger / k
        coreg[se] = [score, p_val]
        
    return coreg