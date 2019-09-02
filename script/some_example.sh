#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "Usage: bash $0 [celltype]" >&2
  exit 1
fi

CELLTYPE=$1
echo 'cell type: '$1

#---------------
# Preprocessing
#---------------
python3.5 Preprocess.py --celltype $CELLTYPE #--stitched

#-----------------------------------
# Non-negative matrix factorization
#-----------------------------------
python3.5 NMF.py --average $CELLTYPE

exit

NUM_STATE=$(wc -l ../data/$CELLTYPE/NMF/sample2state.csv | awk '{ print $1 }')
NUM_STATE="$NUM_STATE-1"
echo "Number of hidden state: $NUM_STATE"
#--------------------------------------
# Co-expression Network Analysis 
#--------------------------------------

# select seRNA-associated genes from co-expression analysis
export ALLOW_WGCNA_THREADS=3
python3.5 fromW_coxnet.py $CELLTYPE $NUM_STATE

exit


#--------------------------------------
# TF enrichment on stage-specific seRNAs. 
#--------------------------------------
python3.5 enrich_TF.py --celltype $CELLTYPE
exit



#------------------------------------------
# TF enrichment on seRNAs-associated genes. 
#------------------------------------------
python3.5 enrich_genes.py --celltype $CELLTYPE --num_state $NUM_STATE
exit


#------------------------------------------
# Co-regulatory analysis 
#------------------------------------------
python3.5 coreg_TF.py --celltype $CELLTYPE --num_state $NUM_STATE
exit


#------------------------------------------
# GO enrichment analysis
#------------------------------------------
python3.5 GO_SemSim.py --type FEA --celltype $CELLTYPE --num_state $NUM_STATE
exit