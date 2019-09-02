from sklearn.cluster import AgglomerativeClustering as AggCluster
import pandas as pd
import os

def determine_stage(mat, exp_profile, n_clusters, outfilepath):
    '''
    Determine stage of differentiation by NMF H matrix.
    Input: H matrix, gene expression profile
    Return: list of gene expression profiles of each stage
    '''
    # Perform Average linkage Clustering
    clustering = AggCluster(n_clusters=n_clusters, linkage='ward')
    clustering.fit(mat.T)
    labels = clustering.fit_predict(mat.T) 
    
    # Create a list of dataframes to store separated exp. of each stage.
    grouped = [pd.Series(exp_profile.index, index=exp_profile.index)] * len(set(labels))
    
    # Record which days is which stages.
    stage_index = []
    for stage_id in set(labels):
        stage_index.append([i for i, x in enumerate(labels) if x == stage_id])
    
    # Seperate gene expression profiles based on clustered labels. 
    for stage in range(len(grouped)):
        for day in stage_index[stage]:
            if day >= 10:
                grouped[stage] = pd.concat([grouped[stage], exp_profile.filter(regex=('day'+str(day)+'.*'))], axis=1)      
            else:
                grouped[stage] = pd.concat([grouped[stage], exp_profile.filter(regex=('day0'+str(day)+'.*'))], axis=1)      

    # Remove redundant column.
    for i in range(len(grouped)):
        grouped[i] = grouped[i].drop(0, axis=1)
    
    # Save expression profile of each stage
    for stage in range(len(grouped)):
        grouped[stage].to_csv(os.path.join(outfilepath, "stage"+str(stage)+"_genes.tsv"), sep="\t", index=True)    

    return grouped

def determine_state(mat, outfilepath):
    '''
    Determine state of each seRNA by NMF W matrix.
    Input: W matrix
    Return: seRNA to state mapping table.
    '''    
    seRNA = []
    states = []
    for rowname, row in mat.iterrows():
        seRNA.append(rowname)
        states.append(row.idxmax())
    map_table = pd.DataFrame([seRNA, states])
    
    map_table.to_csv(outfilepath, sep="\t", index=False)

    return map_table
