import joblib
import matplotlib.pyplot as plt
import pandas as pd
from stitch import stitch_utility as sutil
import networkx as nx
import os, glob

def stitch_main(data, gene, num, path):
    '''
    stitching criteria: 
    1. locus centers within 12500(bp)
    2. on the same chromosome
    3. do not overlap with gene loci
    
    data: list of enhancer expression grouped by chorm.
    gene_loci: for checking overlapping with genes when stitching
    num: ith number of chromosome
    path: output filepath of stitched enhancer expression profiles.
    '''
    se = data[0].copy().iloc[:1, :]  # copy the first row in order to obtain the same header
    
    # Data for grouping
    chrom = data[num]
    chrom = chrom.sort_values('locus')
    
    # Build graph for recording groups
    net = nx.Graph()
    net.add_nodes_from(chrom.index)
    
    # Sorted locus. From the first, find the next point if available.
    indexes = chrom.index.tolist()
    for i in range(0, len(chrom) - 1):
        if chrom.loc[indexes[i+1], 'locus'].sum() - chrom.loc[indexes[i], 'locus'].sum() <= 12500:
            series = chrom.iloc[i, :]
            next_series = chrom.iloc[i + 1, :]
            idx = chrom.iloc[i].name
            next_idx = chrom.iloc[i + 1].name
            new_region = sutil.new_label(pd.DataFrame([series, next_series], index = [idx, next_idx])) 
            if not sutil.overlapped(new_region, gene):
                net.add_edge(idx, next_idx)
    
    # sum up expression levels of each connected components
    for cluster in nx.connected_components(net):
        stitch = chrom.loc[cluster, :]

        # rename label and append it to se 
        se = se.append(pd.Series(stitch.sum(), name = sutil.new_label(stitch)))

    se = se.drop(se.index[:1]) # drop the first row from the beginning copy
    joblib.dump(se, path)
    return (num, path)

def main(data, gene, n_chroms):
    if not os.path.exists('./stitch/stitch_temp/'):
        os.makedirs('./stitch/stitch_temp/')  

    # Run through all dataset
    paths = []
    for num in range(n_chroms):
        paths.append("./stitch/stitch_temp/stitched_" + str(num) + "_.pkl")
    
    for num in range(len(paths)):
        print("Now stitching number %d chrom." % (num))   
        stitch_main(data, gene, num, paths[num]) 
        print("Result: stitched number %d chrom. pickled in %s." % (num, paths[num]))

    merged = sutil.merge(paths)
    for i in glob.glob('./stitch/stitch_temp/*.pkl'):
        os.remove(i)
    return merged

if __name__ == "__main__":
    main()
