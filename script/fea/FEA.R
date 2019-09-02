#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library('topGO'))

dir.name <- dirname(sys.frame(1)$ofile)
source(file.path(dir.name, 'extract.R'), chdir=TRUE)

## Convert module mRNA loci into gene symbols
convert2symbol <- function(celltype, stage, color){
  
  # Load in gene list
  gene_universe = read.table(file.path('../data', celltype, paste('symbol_', celltype, '.tsv', sep='')), 
                             sep='\t', header=FALSE, row.names=1)
  gene_universe = as.factor(gene_universe[['V2']])
  
  gene_symbol = extract_symbol_from_modules(color, celltype, stage)

  allgenes = factor(as.integer(gene_universe %in% gene_symbol))
  names(allgenes) = gene_universe
  
  return(allgenes)
}

## Wrapper function of Functional Enrichment Analysis
# For mRNA-to-mRNA coxnet modules
TopGO <- function(gene_symbol, celltype, stage, color){
    
  gene_universe = read.table(file.path('../data', celltype, 'Super_Enhancers',
                                       paste('symbol_', celltype, '.tsv', sep='')), 
                             sep='\t', header=TRUE, row.names=1)
  gene_universe = as.factor(gene_universe[['symbol']])
  
  allgenes = factor(as.integer(gene_universe %in% gene_symbol))
  names(allgenes) = gene_universe
  

  # Create topGO object
  GOdata <- new("topGOdata", description=paste('GO of stage', stage, 'of module', color), 
                allGenes=allgenes, ontology="BP",
                annot=annFUN.org, mapping='org.Hs.eg.db', ID='alias')

  # Run enrichment analysis
  resultFisher <- runTest(GOdata, algorithm="classic", statistic="fisher")

  resultKS <- runTest(GOdata, algorithm="classic", statistic="ks")
  
  # Plot method comparing table of top20 GO term
  allRes <- GenTable(GOdata, Fis=resultFisher, KS=resultKS, 
                        orderBy="Fis", ranksOf="KS", topNodes=100, numChar=100000)
  
  # Save tables
  dir.create(file.path('../data', celltype, 'FEA'), showWarnings=FALSE)
  dir.create(file.path('../data', celltype, 'FEA', 'stitched'), showWarnings=FALSE)
  dir.create(file.path('../data', celltype, 'FEA', 'stitched', 
                       paste('state', stage, sep='')), showWarnings=FALSE)
  filepath = file.path('../data', celltype, 'FEA', 'stitched', 
                       paste('state', stage, sep=''))

  write.table(allRes, file.path(filepath, paste(color, 'bp_result.csv', sep='_')), 
              row.names=FALSE, sep=',')

  showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')
    
  graph_dir = file.path('..', 'graph', celltype, 'FEA', 'stitched')
  dir.create(graph_dir, showWarnings=FALSE)
  dir.create(file.path(graph_dir, paste('state', stage, sep='')), showWarnings=FALSE)
    
  printGraph(GOdata, resultFisher, firstSigNodes=5, 
             fn.prefix=file.path(graph_dir, paste('state', stage, sep=''), color),
             useInfo = "all", pdfSW = TRUE)
}

analyze <- function(celltype, stage_id, color_list_path){
  module_color = read.table(color_list_path)
  len = length(module_color[[1]])
  allgenes = list()

  for (i in c(seq(1,len))){
    allgenes[[i]] = convert2symbol(celltype=celltype, stage=stage_id, color=module_color[[1]][i])  
  }

  for (i in c(seq(1,len))){
    TopGO(allgenes[[i]], celltype=celltype, stage=stage_id, color=module_color[[1]][i])  
  }
}

# Execute when invoke this script
main <- function() {
  if (length(args) == 3){
    analyze(args[1], args[2], args[3])
    quit(save='no')
  }
  print('Usage: Rscript FEA.R <cell type> <stage number> <file path of module color list>')
  quit(save='no')
}
