# Set global configuration.
# setwd("~/github/Research/script/network/")
#options(stringsAsFactors=F)

# Read in the R libraries 
library(MASS)
library(class)   
library(cluster) 
library(impute)   
library(WGCNA) 
source('./WGCNA/NetworkFunctions.R')

#====================================================================================
# Read in expression data
if(exists('read_exprdata')) rm(read_exprdata);
read_exprdata <- function(infilepath, sep='\t'){

  data = read.csv(infilepath, sep = sep)
  rownames(data) = data$X
  data = data[,c(-1)]
  Exprdata = apply(data, c(1,2), as.numeric)
  Exprdata = t(Exprdata)

  return(data.frame(Exprdata))
}

#====================================================================================
if(exists('choose_power')) rm(choose_power);
choose_power <- function(Exprdata, celltype, stage_id, powers=c(seq(1,10,by=1),seq(12,20,by=2))){
    
  #Investigate soft thesholding with the power adjacency function 
  RpowerTable = pickSoftThreshold(Exprdata, powerVector=powers)
  suggest_power = RpowerTable[[1]]
  RpowerTable = RpowerTable[[2]]

  #  Plot soft threshold relative to signed R^2 and Mean k.
  dir.create(file.path('../graph', celltype), showWarnings=FALSE)

  jpeg(filename=file.path('../graph', celltype, paste('stage ', stage_id, ' soft_threshold.jpg', sep='')), width=1600, height=1200)
  cex1 = 1.5
  par(mfrow=c(1,2))
  plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],xlab="Soft Threshold (power)",
    ylab="Scale Free Topology Model Fit,signed R^2",type="n")
  text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], labels=powers, cex=cex1, col="red")

  # this line corresponds to using an R^2 cut-off of h
  for (i in c(0.9, 0.8, 0.7)){
    abline(h=i,col="red")
    text(2, i, labels=i)
  }

  plot(RpowerTable[,1], RpowerTable[,5],xlab="Soft Threshold (power)",
       ylab="Mean Connectivity", type="n")
  text(RpowerTable[,1], RpowerTable[,5], labels=powers, cex=cex1, col="red")
  dev.off()

  Rpower = list('suggestion'=suggest_power, 'table'=RpowerTable)
  return(Rpower)
}

#====================================================================================
# Create a scale free topology plot.
# The black curve corresponds to scale free topology
if(exists('check_scalefree')) rm(check_scalefree);
check_scalefree <- function(Exprdata, beta, celltype, stage_id, truncated_flag=F){
    
  Connectivity = softConnectivity(Exprdata, power=beta)-1
 
  dir.create(file.path('../graph', celltype), showWarnings=FALSE) 
  jpeg(filename=file.path('../graph', celltype, paste('stage ', stage_id, ' check_scalefree.jpg', sep='')), width=800, height=600)
  par(mfrow=c(1,1))
  scaleFreePlot(Connectivity, main=paste("soft threshold, power=",beta), truncated=truncated_flag)
  dev.off()
}

#====================================================================================
# Calculate adjacency matrix.
# Choose beta to preserve a certain level of k.
if(exists('calculate_TOM')) rm(calculate_TOM);
calculate_TOM <- function(Exprdata, beta){

  ADJ = adjacency(Exprdata, power=beta)
  # Computes the topological overlap matrix based on the adjacency matrix.
  dissTOM = TOMdist(ADJ)
  rm(ADJ)
  gc()  
  # Carry out hierarchical clustering with the TOM matrix. 
  hierTOM = hclust(as.dist(dissTOM),method="average")

  adj = list('hierTOM'=hierTOM, 'dissTOM'=dissTOM)
  return(adj)
}

#====================================================================================
# Module Detection
# The function cutreeStatistColor colors each gene by the branches that
# result from choosing a particular height cut-off.
# GREY IS RESERVED to color genes that are not part of any module.
if(exists('detect_module')) rm(detect_module);
detect_module <- function(ADJ, celltype, stage_id, min_modulesize=50){
  # Dynamic tree cut:
  unmergedLabels = cutreeDynamicTree(dendro=ADJ$hierTOM, minModuleSize = min_modulesize)
  colorh1 = labels2colors(unmergedLabels)

  # Plot dendrogram
  dir.create(file.path('../graph', celltype), showWarnings=FALSE)
  jpeg(filename=file.path('../graph', celltype, paste('stage ', stage_id, ' Module_dendrogram.jpg', sep='')), width=1600, height=1200)
  par(mfrow=c(2,1),mar=c(2,4,1,1))
  plot(ADJ$hierTOM, main="Cluster Dendrogram", labels=F, xlab="", sub="");
  plotColorUnderTree(ADJ$hierTOM,colors=data.frame(colorDynamic=colorh1))
  title("Module (branch) color")
  dev.off()

  print(table(colorh1))
  modules = list('colorh1'=colorh1, 'unmergedLabels'=unmergedLabels)
  
  return(modules)
}

#====================================================================================
# We also propose to use classical multi-dimensional scaling plots 
# for visualizing the network. Here we chose 2 scaling dimensions
# This also takes about 10 minutes...
other_plot <- function(dissTOM, hierTOM, moduleLabels, colorh1, celltype, stage_id){

  moduleColors = labels2colors(moduleLabels)

  dir.create(file.path('../graph', celltype), showWarnings=FALSE)
  #jpeg(filename=file.path('../graph', celltype, paste('stage ', stage_id, ' TOMplot.jpg', sep='')), width=1600, height=1200)
  #TOMplot(dissTOM, hierTOM, moduleColors)
  #dev.off()

  cmd1 = cmdscale(as.dist(dissTOM),2)
  
  jpeg(filename=file.path('../graph', celltype, paste('stage ', stage_id, ' MDSplot.jpg', sep='')), width=1600, height=1200)
  par(mfrow=c(1,1))
  plot(cmd1, col=as.character(colorh1), main="MDS plot",xlab="Scaling Dimension 1",ylab="Scaling Dimension 2")
  dev.off()
}
#====================================================================================
# Relations between eigengenes of each modules.
# Summarize each module by its first eigengene (referred to as principal components).
# and then correlate these module eigengenes with each other.
if(exists('eigengenes')) rm(eigengenes);
eigengenes <- function(Exprdata, unmergedLabels, colorh1, ADJ, celltype, stage_id, height=0.6){
  
  MElist = moduleEigengenes(Exprdata, colorh1)
  datME = MElist[[1]]
  # We define a dissimilarity measure between the module eigengenes that keeps track of the sign of the correlation between the module eigengenes.
  dissimME = 1-(t(cor(datME, method="p")))/2
  hclustdatME = hclust(as.dist(dissimME), method="average" )
  
  dir.create(file.path('../graph', celltype), showWarnings=FALSE)
  jpeg(filename=file.path('../graph', celltype, paste('stage ', stage_id, ' eigengenes.jpg', sep='')))
  par(mfrow=c(1,1))
  plot(hclustdatME, main="Clustering tree based on the module eigengenes of modules")
  signif(cor(datME, use="p"), 2)
  abline(h=height, col='red')
  dev.off()

  # Merge similar modules based on eigengenes similarity.
  merge = mergeCloseModules(Exprdata, unmergedLabels, cutHeight=height, verbose=2)
  moduleLabels = merge$colors
  moduleColors = labels2colors(moduleLabels)
  #consMEs = merge$newMEs
  # Plot unmerged and merged modules.
  jpeg(filename=file.path('../graph', celltype, paste('stage ', stage_id, ' Merged_modules.jpg', sep=''))
       , width=800, height=600)
  plotDendroAndColors(ADJ$hierTOM, moduleColors, c('merged'), rowText=moduleColors,
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()

  # Now we create scatter plots of the samples (arrays) along the module eigengenes.
  datMEordered = datME[,hclustdatME$order]
  
  jpeg(filename=file.path('../graph', celltype, paste('stage ', stage_id, ' eigengene_relation.jpg', sep=''))
       , width=1600, height=1200)
  pairs(datMEordered, upper.panel=panel.smooth, lower.panel=panel.cor, 
        diag.panel=panel.hist ,main="Relation between module eigengenes")
  dev.off()

  print('Eigengenes plots saved and module merged.')
  module = list('label'=moduleLabels, 'color'=moduleColors, 'data'=datME)
  return(module)
}
#====================================================================================
# The following produces heatmap plots for each module.
# Here the rows are genes and the columns are samples.
# Well defined modules results in characteristic band structures since the corresponding genes are 
# highly correlated.
if(exists('module_heatmap')) rm(module_heatmap);
module_heatmap <- function(Exprdata, colorh1, celltype, stage_id){

  # Create necessary directory to store graph
  dir.create(file.path('../graph', celltype), showWarnings=FALSE) 
  dir.create(file.path('../graph', celltype, 'module'), showWarnings=FALSE)
  dir.create(file.path('../graph', celltype, 'module', paste('stage', stage_id, sep='')), showWarnings=FALSE)
  
  # Plot heatmap for each modules
  ClusterSamples = hclust(dist(Exprdata),method="average") 
  module_lst = names(table(colorh1))
  for (x in module_lst){
    which.module = x
    jpeg(filename=file.path('../graph', celltype, 'module', paste('stage', stage_id, sep=''), 
                            paste(which.module, stage_id, '_heatmap.jpg', sep='')), width = 800, height = 600)
    par(mfcol=c(1,1), cex=2)
    moduledata = t(scale(Exprdata[ClusterSamples$order,][,colorh1==which.module]))
    moduledata = as.matrix(data.frame(moduledata)[order(moduledata[,1]),])
    plotMat(moduledata,nrgcols=30,rlabels=T,clabels=T,rcols=which.module, main=which.module)
    dev.off()
  }
  # Plot heatmap of modules concatenated.
  jpeg(filename=file.path('../graph', celltype, 'module', paste('stage', stage_id, sep=''), 
                          paste('Allmodules_heatmap.jpg', sep='')), width = 800, height = 600)
  par(mfcol=c(1,1), cex=2)
  moduledata = t(scale(Exprdata[ClusterSamples$order,]))
  moduledata = as.matrix(data.frame(moduledata)[order(moduledata[,1]),])
  plotMat(moduledata,rlabels=T,clabels=F,rcols=colorh1,main='WGCNA modules')
  dev.off()
  
  # Save loci in each module in separated lists.
  dir.create(file.path('../data', celltype), showWarnings=FALSE) 
  dir.create(file.path('../data', celltype, 'module'), showWarnings=FALSE)
  dir.create(file.path('../data', celltype, 'module', paste('stage', stage_id, sep='')), showWarnings=FALSE)
  for (x in module_lst){
    which.module = x
    write.table(t(Exprdata[ClusterSamples$order,][,colorh1==which.module]), 
                file.path('../data', celltype, 'module', paste('stage', stage_id, sep=''),
                         paste(which.module, stage_id, "_module_exp.tsv", sep="")), sep="\t")
  }
}

#================================================================
if(exists('select_hub_genes')) rm(select_hub_genes);
select_hub_genes <- function(datExpr, colorh, beta, celltype){
  hubs = chooseTopHubInEachModule(datExpr, colorh, power=beta, omitColors='grey', type='unsigned')
  symbol = read.csv(file.path('..', 'data', celltype, paste('symbol_', celltype, '.tsv', sep='')), 
                    sep='\t', header=FALSE, row.names=1)
  hubsymbol = data.frame(symbol[hubs, ])
  rownames(hubsymbol) = rownames(data.frame(hubs))
  return(hubsymbol)
}
#=======================================================
combine_se2gene <- function(celltype){
  # Combine super-enhancer to gene mapping table
  se2gene = read.table(file.path('..', 'data', celltype, 'se2gene.bed'), sep='\t')[,-7]
  se2gene$SE = do.call(paste, c(se2gene[c("V1", "V2")], sep = ".")) 
  se2gene$SE = do.call(paste, c(se2gene[c("SE", "V3")], sep = "."))
  se2gene$gene = do.call(paste, c(se2gene[c("V4", "V5")], sep = ":")) 
  se2gene$gene = do.call(paste, c(se2gene[c("gene", "V6")], sep = "..")) 
  se2gene = se2gene[,c('SE', 'gene')]
  rownames(se2gene) = se2gene$SE
  return(se2gene)
}
