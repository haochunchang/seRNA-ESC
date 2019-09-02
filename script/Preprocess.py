## @package Preprocess
# Extract enhancer and gene expression profile.
# Normalize gene and stitch & get super-enhancers.
#
# 
# Usage: python Preprocess.py [celltype]

from preprocess import data_preprocess as pre
from preprocess import map_gene
from stitch import stitch_main, stitch_utility
import joblib
import sys, argparse, os.path, math
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import pandas as pd
import os

from misc.utils import ensure_dir

pandas2ri.activate()

## Wrapper of stitching enhancer expression profile.
# @param enhancers: enhancer expression profile: pandas DataFrame
# @param celltype: string of celltype, e.g. 'HES3_GFP_ESC'
# Output: .tsv file of stitched enhancer expression profile.
# 
def stitch_enhancers(enhancers, celltype):
    
    # Group enhancers by chromosome
    enhancers, n_chroms = stitch_utility.group_by_chrom(enhancers)
    
    # Load in gene loci data in hg19 whole genome    
    gene = pd.read_csv(os.path.join("..", "data", celltype, "Super_Enhancers", "gene_loci.bed"), sep="\t")
    gene = gene.drop_duplicates()
    print(gene.shape)
    #gene = None

    stitched = stitch_main.main(enhancers, gene, n_chroms)        
    print("Stitched enhancers have: {}".format(stitched.shape))
    stitched.to_csv(os.path.join("..", "data", celltype, "Super_Enhancers","stitched_enhancer.tsv"), sep="\t")
    
    return stitched
    
## Wrapper of mapping genes to stitched enhancers.
# @param stitched: stitched enhancers expression profile: pandas DataFrame
# @param celltype: string of celltype, e.g. 'HES3_GFP_ESC'
# Output: se2gene.bed, genes to stitched enhancers mapping table. 
#
def mapping(stitched, celltype):

    #stitched = pd.read_csv(os.path.join("../data/"+celltype, "stitched_enhancer.tsv"), sep="\t", index_col=0, header=list(range(3)))
    genefilepath = os.path.join("..", "data", celltype, "Super_Enhancers", "gene_loci.bed")
    enhancerfilepath = os.path.join("..", "data", celltype, "Super_Enhancers", "enhancer_loci.bed")
    outfilepath = os.path.join("..", "data", celltype, "Super_Enhancers", "se2gene.bed")
    
    pre.ToBED(stitched, filepath=enhancerfilepath, loci_in='index')
    map_gene.mapGeneToEnhancer(genefile=genefilepath, enhancerfile=enhancerfilepath, outfilepath=outfilepath)    
    
    print("Genes from "+genefilepath+" mapped to "+enhancerfilepath+". Stored in "+outfilepath)

## Wrapper of super-enhancer selection.
# Output: super-enhancer expression profile .csv file 
#
def select_super(celltype):
    r = robjects.r
    r.source("./preprocess/callSuper_main.R", chdir=True)
    
    super_enhancer = r.Select_Super(outFolder=os.path.join("..", "data", celltype, "Super_Enhancers",
                                                           "super-enhancer"), \
                                    enhancerFile=os.path.join("..", "data", celltype, "Super_Enhancers",
                                                              "stitched_forcall.tsv"),\
                                    enhancerName=celltype)
    super_enhancer = super_enhancer['REGION_ID']
    return super_enhancer


## Main function
def main(celltype, stitched=False):
    
    # Select enhancers & genes
    ensure_dir(os.path.join("../data/", celltype, "Super_Enhancers"))
    enhancers, genes = pre.time_course_prep(Infilepath=os.path.join("..", "data", "raw", celltype), \
                                            Outfilepath=os.path.join("..", "data", celltype, "Super_Enhancers")) 
    pre.ToBED(genes, os.path.join("..", "data", celltype, "Super_Enhancers", "gene_loci.bed"))
    print('Enhancers: {}, Genes: {}'.format(enhancers.shape, genes.shape))

    # Stitch enhancers
    if stitched:
        print("Start Stitching")
        stitched = stitch_enhancers(enhancers, celltype)
    else:
        print("Skip Stitching")
        stitched = enhancers

    # Gene mapping
    mapping(stitched, celltype)
     
    # Select Super-enhancers
    pre.convert2ROSE(stitched, outfilepath=os.path.join("..", "data", celltype, "Super_Enhancers"))
    
    super_enhancer = select_super(celltype)
    super_enhancer.to_csv(os.path.join("..", "data", celltype, "Super_Enhancers",
                                       "super-enhancer_loci.tsv"), 
                          index=False, sep='\t')
    
    super_enhancer = pd.read_csv(os.path.join("..", "data", celltype, "Super_Enhancers",
                                              "super-enhancer_loci.tsv"), 
                                 sep='\t', header=None, names=['REGION_ID'])
    
    se_profile = stitched.loc[super_enhancer['REGION_ID'], :]
    se_profile.to_csv(os.path.join("..", "data", celltype, "Super_Enhancers",
                                   "super-enhancer_exp.tsv"), sep="\t", index=True)
    print("Super-enhancer profile stored in " + os.path.join("..", "data", celltype, "super-enhancer_exp.tsv"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='Preprocess')
    parser.add_argument("--stitched", help="Flag indicating whether stitch enhancers or not.", action="store_true")
    parser.add_argument('--celltype', required=True, type=str) 
    args = parser.parse_args()
    
    ensure_dir(os.path.join('..', 'data', args.celltype))
    
    main(args.celltype, args.stitched)
