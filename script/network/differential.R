library(limma)
library(data.table)
#=============================================================
#source('./network/binet_utils.R', chdir=TRUE)
#celltype <- 'HES3_GFP_ESC_UQ'
#stage_id <- 0
#exp_profile <- Load_gene_stages(celltype, 3)
#p_value <- 1
#=============================================================
if(exists('calculate_differential')) rm(calculate_differential);
calculate_differential <- function(exp_profile, stage_id, p_value_threshold=0.05, lfc=0.6){

  target = exp_profile[[stage_id+1]]
  other_profiles = exp_profile[-c(stage_id+1)]
  others = data.frame(rbindlist(other_profiles))
  profile = rbind(others, target)

  design = cbind(others=1, target=c(rep(0, nrow(others)), rep(1, nrow(target))))
  fit = lmFit(t(profile), design)
  fit = eBayes(fit)
  
  ntop = ncol(profile)
  toptable = topTable(fit, coef="target", number=ntop, sort.by="p", lfc=lfc, p.value=p_value_threshold)
  return(toptable)
  
}

if(exists('show_differential')) rm(show_differential);
show_differential <- function(exp_profile, stage_id, celltype, p_value=0.05){

  all_top = calculate_differential(exp_profile, stage_id, p_value_threshold=p_value)
  write.table(all_top, file=file.path('../data/', celltype, paste('DE_TABLE_', stage_id, '.tsv', sep='')), sep='\t')
  top_gene_names = rownames(all_top)
  exp_profile = exp_profile[[stage_id+1]][,top_gene_names]
 
  # Convert DE gene loci to gene symbols 
  loci2symbol = read.table(file.path('../data/', celltype, paste('symbol_', celltype, '.tsv', sep='')), 
                           sep='\t', header=FALSE, row.names=1)
  
  top_genes = sapply(top_gene_names, function(x) sub(".", ":", x, fixed=TRUE))
  symbol = loci2symbol[top_genes,]
  write.table(symbol, file=file.path('../data/', celltype, paste('DE_stage', stage_id, '_symbol.tsv', sep='')), sep='\t')
 
  return(exp_profile)

}

