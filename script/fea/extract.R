#celltype <- 'HES3_GFP_ESC'
#module_color <- 'turquoise'
#stage_id <- 0
#============================================================
# Extract gene symbols given module color from WGCNA, celltype and 
# stage number. Used for mRNA-to-mRNA coxnet
if(exists('extract_symbol_from_modules')) rm(extract_symbol_from_modules);
extract_symbol_from_modules <- function(module_color, celltype, stage_id){

  loci2symbol = read.table(file.path('../data', celltype, paste('symbol_', celltype, '.tsv', sep='')), 
                            sep='\t', header=FALSE, row.names=1)
  module = read.table(file.path('../data', celltype, 'module', paste('stage', stage_id, sep=''), 
                                paste(module_color, stage_id, '_module_exp.tsv', sep='')), sep='\t')

  module_genes = rownames(module)
  #module_genes = sapply(module_genes, function(x) sub(".", ":", x, fixed=TRUE))
  module_genes = as.character(module_genes)
  symbol = loci2symbol[module_genes,]
  write.table(symbol, file=file.path('../data', celltype, 'module', paste('stage', stage_id, sep=''),
                                    paste(module_color, stage_id, '_module_symbol.tsv', sep='')), sep='\t')
  return(symbol)  
}

#============================================================
# Extract mRNA from W-selected exp profile given state number.
# from: 'top_seRNA'(default): load in top-j mRNA selected based on top-i seRNA from state2seRNA
#       'w': load in top-j mRNA selected based on mRNA NMF W matrix. 
#celltype <- 'HES3_GFP_ESC_UQ'
#state_id <- 0
#
if(exists('extract_symbol_from_selected')) rm(extract_symbol_from_selected);
extract_symbol_from_selected <- function(celltype, state_id, from='top_seRNA'){

  loci2symbol = read.table(file.path('../data', celltype, paste('symbol_', celltype, '.tsv', sep='')), 
                            sep='\t', header=FALSE, row.names=1)
  if(from=='w'){
    loci = read.table(file.path('../data', celltype, paste('selected_mRNA_state', 
                                                           state_id, '.csv', sep='')), sep=',')
    module_genes = loci[-1, dim(loci)[[2]]]
    module_genes = as.character(module_genes)
    symbol = loci2symbol[module_genes,]
    write.table(symbol, file=file.path('../data', celltype, paste('selected_mRNA_state', state_id, '_symbol.csv', sep='')), sep=',')
    
  } else if(from=='top_seRNA'){
    loci = read.table(file.path('../data', celltype, 
                                paste('state', state_id, '_top_mRNA_from_seRNA.csv', sep='')), sep=',')
    module_genes = loci[-1, 1]
    module_genes = as.character(module_genes)
    symbol = loci2symbol[module_genes,]
    write.table(symbol, file=file.path('../data', celltype, 
                                       paste('state', state_id, '_top_mRNA_from_seRNA_symbol.csv', sep='')), sep=',')
  }

 return(symbol)  
}

