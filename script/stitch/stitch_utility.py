import joblib
import pandas as pd

def overlapped(enhancer, gene):
    '''
    Input: string of enhancer locus (ex: chr1:100055-100080)
    gene_loci: hg19 whole gene loci data
    Return True if enhancer locus is overlapped with genes, else False
    '''
    
    chrom = enhancer.split(':')[0]
    start = enhancer.split(':')[1].split('-')[0]
    end = enhancer.split('-')[1]

    # Select rows which locates within enhancer locus
    # gene = pd.read_csv('../data/gene_loci.bed', sep = "\t")
    gene = gene[gene['#CHROM'] == chrom]
    gene_start = gene[gene['START'] >= int(start)]
    gene_end = gene_start[gene_start['STOP'] <= int(end)]
    
    if len(gene_end) > 0:
        return True
    else:
        return False
               
def new_label(proximal):   
    '''
    Return new label name of the combined row.
    Example: chr10:000000-000001 + chr10:000003-000010 = chr10:000000-000010
    '''
    upper = proximal[proximal['locus'] == max(proximal['locus'])].index[0].split("..")[1]
    lower = proximal[proximal['locus'] == min(proximal['locus'])].index[0].split("..")[0]
    return lower + "-" + upper

def get_center(coor):
    '''Return given locus midpoint'''
    pos = coor.split(":")[1].split("..")
    return (int(pos[0]) + int(pos[1])) / 2

 
def group_by_chrom(data):    
    '''
    Group by chromosome
    Return: a list containing expression profiles of each chromosomes.
    '''
    data['separation'] = [False for i in data.index]
    data['locus'] = [get_center(index) for index in data.index]
    data['chrom'] = [i.split(':')[0] for i in data.index]

    chroms = []
    data.sort_values('chrom', inplace = True)
    groups = data.groupby('chrom')
    for name in groups.groups.keys():
        chroms.append(groups.get_group(name))     
    
    joblib.dump(chroms, "./stitch/stitch_temp/group_by_chrom.pkl")
    
    return chroms, len(groups.groups.keys())

def merge(paths=[]):
    '''
    Merge several joblib of dataframes into a single dataframe.
    '''
    result = pd.DataFrame()
    for i in range(len(paths)):
        se = joblib.load(paths[i])
        result = result.append(se)
    result.drop(['locus', 'separation', 'chrom'], axis = 1, inplace = True)
    
    joblib.dump(result, "./stitch/stitched_all.pkl")
    return result    


if __name__ == "__main__":
    data = joblib.load("group_by_chrom.pkl")
    print(len(data))
