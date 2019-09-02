#========================================================================================
#======================= Library of Data Preprocessing Functions ========================
#========================================================================================

import joblib 
import pandas as pd
import numpy as np
import re, os.path

def time_course_prep(Infilepath = "../data/raw/HES3_GFP_ESC.tsv", Outfilepath = "../data/"):
    '''
    Do: group by days, rep
    Infilepath: relative path of time-course data .tsv file
    Outfile: .tsv file of processed pandas DataFrame
    '''
    data = pd.read_csv(Infilepath+".tsv", sep = "\t", index_col = 1)
        
    # p1@, p2@, ..., are not enhancer peaks.
    not_enhancer = [index for index in data.index if index[:2] != "p@"]
    not_promoter = [index for index in data.index if index[:2] == "p@"]

    enhancer = data.drop(not_enhancer)
    enhancer = enhancer.drop('00Annotation', axis=1)
    enhancer.index = [index.split(",")[0][2:]  for index in enhancer.index] 
    
    gene = data.drop(not_promoter)   
    gene_locus = [g[:-2] for g in gene['00Annotation']]
    
    # Retain gene original peak ranks
    peak_rank = pd.DataFrame(gene_locus, index=gene.index)
    peak_rank['short_description'] = peak_rank.index
    #peak_rank.to_csv(os.path.join(Outfilepath, "genes_peak_rank.csv"), index=True)
    
    # Drop low-rank genes
    thres = 2
    print('Before dropping low-rank(rank > {}) genes: {}'.format(thres, gene.shape))
    #peak_rank = pd.read_csv(os.path.join('..', 'data', celltype, 'genes_peak_rank.csv'), index_col=1)
    peak_rank['is_low_rank'] = peak_rank['short_description'].apply(lambda x: is_low_peaks(x, thres=thres))
    low_ranks = [x for x in peak_rank.index if peak_rank.loc[x,'is_low_rank'] == True]   
    gene = gene.drop(low_ranks)
    print('After dropping low-rank(rank > {}) genes: {}'.format(thres, gene.shape))

    gene_locus = [g[:-2] for g in gene['00Annotation']]
    gene = gene.drop('00Annotation', axis=1)
    gene.index = [re.sub('p\w@|p\w\w@', '', index) for index in gene.index]
    
    # shrink column names into only days and rep 
    # construct header of days, rep 
    header = [[], []]
    for col in list(enhancer.columns.values): 
        col = col.split(",")
        days = col[2]
        rep = col[-1].split(".")[0]
        header[0].append(days)
        header[1].append(rep)
        
    enhancer.columns = pd.MultiIndex.from_tuples(list(zip(*header)))
    gene.columns = pd.MultiIndex.from_tuples(list(zip(*header)))
    
    gene2symbol = pd.Series([sym for sym in gene.index], index=gene_locus)
    gene.index = [locus for locus in gene_locus]

    # Save cleaned gene and enhancer profiles
    enhancer.to_csv(os.path.join(Outfilepath, "enhancers_" + Infilepath.split("/")[-1]+'.tsv'), sep="\t")    
    gene.to_csv(os.path.join(Outfilepath, "genes_" + Infilepath.split("/")[-1]+'.tsv'), sep="\t")
    
    # clean symbol before save
    tmp = pd.DataFrame()
    tmp['symbol'] = gene2symbol
    gene2symbol = clean_sym(tmp)
    gene2symbol.to_csv(os.path.join(Outfilepath, "symbol_" + Infilepath.split("/")[-1]+'.tsv'), sep="\t",
                      header=True)    
    
    print("Processed files are stored in " + Outfilepath)
    return enhancer, gene

def clean_sym(data):
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    # Suppress R script warnings
    import warnings
    from rpy2.rinterface import RRuntimeWarning
    warnings.filterwarnings("ignore", category=RRuntimeWarning)

    # only keep the first symbol if multiple exists
    data = data['symbol'].str.split(',').apply(lambda x: x[0])

    # split out into hgnc symbols and ensembl transcript ID
    ensembl_t = data[data.str.startswith('ENST')]
    hgnc_sym = data[~data.str.startswith('ENST')]
    
    # Convert ensembl_t to hgnc_symbol if possible
    with localconverter(robjects.default_converter + pandas2ri.converter):
        ensembl_t_in_r = robjects.conversion.py2rpy(ensembl_t)
    #ensembl_t_in_r = pandas2ri.py2ri(ensembl_t)
    r_sess = robjects.r
    r_sess.source('./misc/utils.R')
    converted_sym = r_sess.convert_ensembl_to_hgnc(ensembl_t_in_r)
    
    with localconverter(robjects.default_converter + pandas2ri.converter):
        converted_sym = robjects.conversion.rpy2py(converted_sym)
    #converted_sym = pandas2ri.ri2py(converted_sym)
    converted_sym.index = converted_sym['ensembl_transcript_id']
    
    for i in range(ensembl_t.shape[0]):
        enst = ensembl_t[i]
        if enst in list(converted_sym.index):
            hgnc = converted_sym.loc[enst]['hgnc_symbol']
            ensembl_t[i] = hgnc

    # Combine back  
    cleared = pd.concat([hgnc_sym, ensembl_t], axis=0)
    return cleared


def is_low_peaks(index, thres=5):
    '''
    Given peaks + locus info, return True when ranking > thres.
    Example: input: p1@LOXL4, thres=5, since p'1' < 5, return False.
    '''
    if index[1:].split('@')[0] != '':
        if int(index[1:].split('@')[0]) < thres:
            return False
    return True


def ToBED(data, filepath, loci_in = "index"):
    '''
    Convert expression data into BED format.
    loci_in: given column name to extract loci data, defaut is from index.
    Return "filename"(any given string).bed 
    '''
    if loci_in == "index":
        # Extract chromosome, start location, end location from dataframe index
        print('Extracting loci data from index.')
        chrom = [index.split(":")[0] for index in data.index]
        if len(data.index[0].split("-")) > 1:
            start = [int(index.split(":")[1].split("-")[0]) for index in data.index]
            end = [int(index.split("-")[1]) for index in data.index]
        else:
            start = [int(index.split(":")[1].split("..")[0]) for index in data.index]
            end = [int(index.split("..")[1]) for index in data.index]

        bed = pd.DataFrame(list(zip(chrom, start, end)), columns=['#CHROM', 'START', 'STOP'])
    
    else:
        try:
            # loci info. is stored in the given column
            print('Extracting loci data from given column name.')
            chrom = [value.split(":")[0] for value in list(data[loci_in])]
            start = [int(value.split(":")[1].split("..")[0]) for value in list(data[loci_in])]
            end = [int(value.split("..")[1]) for value in list(data[loci_in])]
            bed = pd.DataFrame(list(zip(chrom, start, end)), columns=['#CHROM', 'START', 'STOP'])
        except ValueError:
            return "Column name not exists"

    # sort by chrom and start site
    bed = bed.sort_values(['#CHROM', 'START'])

    # save the dataframe into .bed file
    bed.to_csv(filepath, sep="\t", index=False)
    return bed

def convert2ROSE(data, outfilepath='../data/'):
    '''
    Convert data into ROSE-acceptable format.
    data: stitched enhancer expression profile	
    columns: chrom, start, end, region, expression sum
    rows: stitched enhancer loci
    '''
    chrom = [index.split(":")[0] for index in data.index]
    try:
        start = [index.split(":")[1].split("-")[0] for index in data.index]
        end = [index.split("-")[1] for index in data.index]
    except IndexError:
        start = [index.split(":")[1].split("..")[0] for index in data.index]
        end = [index.split("..")[1] for index in data.index]

    exp_sum = list(data.sum(axis=1))
    bed = pd.DataFrame(list(zip(chrom, start, end, data.index, exp_sum)), columns=['CHROM', 'START', 'STOP', 'REGION_ID', 'Expression_Sum'])

    bed.to_csv(os.path.join(outfilepath, "stitched_forcall.tsv"), sep="\t", index=False, index_label=False)
    return bed

def log_transform(data):
    '''
    Apply log2 transformation.
    Plus 1 to all values in case 0 exists in the data.
    '''
    data = data + 1
    data = np.log2(data)  
    return data


if __name__ == "__main__":   
    pass 
