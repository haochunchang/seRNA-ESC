import os
import pandas as pd

def mapGeneToEnhancer(genefile, enhancerfile, outfilepath):
    '''
    Map nearest gene(not-overlapping) to each stitched enhancer locus,
    Select downstream genes first, if ties occur, choose the one that is downstream.
    Genefile: ../data/gene_loci.bed
    enhancerfile: ../data/enhancer_loci.bed
    outfilepath: location where output .bed file store
                 Genes and stitched enhancer mapping table.
    '''
    os.system("bedtools closest -D ref -a "+enhancerfile+" -b "+genefile+" -io -fd -t last > "+outfilepath)


def mapHgncToGene(hgncfile, enhancerfile): 
    '''
    After assigning genes to stitched enhancers, map hgnc_id to each gene locus
    hgnc id reference: ../data/gene_locus2hgnc.bed
    Stitched enhancer loci: ../data/SE2gene.bed
    '''    
    hgnc = pd.read_csv(hgncfile, sep = "\t")
    enhancer = pd.read_csv(enhancerfile, sep = "\t")
    enhancer['hgnc_id'] = ["" for i in range(len(enhancer.index))]
    
    gene_loci = enhancer.loc[:, ['gene_locus', 'gene_start', 'gene_end']]
    hgnc_loci = hgnc.loc[:, ['0', '1', '2']]
    
    for i, series in gene_loci.iterrows():
        result = hgnc_loci[hgnc_loci['0'] == series['gene_locus']]
        result = result[result['1'] == series['gene_start']]
        result = result[result['2'] == series['gene_end']]
        
        if len(result.index) == 1:
            idx = result.index[0]
            enhancer.loc[i, 'hgnc_id'] = str(hgnc.iloc[idx, 4])
        else:
            continue 
    enhancer.to_csv("hgnc_mapped.tsv", sep = "\t", index = False)
        

if __name__ == "__main__":
    enhancers = "../data/SE2gene.bed"
    symb = "../data/gene_locus2hgnc.bed"
    mapHgncToGene(symb, enhancers)
