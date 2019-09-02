source("./callSuper_function.R")
#setwd("~/research/Mining Key Regulators/script/")  

#============================================================================
#===================SUPER-ENHANCER CALLING AND PLOTTING======================
#============================================================================

Select_Super <- function(outFolder="../data/super-enhancer/", enhancerFile="../data/stitched_forcall.tsv", enhancerName="ES_to_cardiomyocyte"){

  stitched_regions <- read.delim2(file= enhancerFile, sep="\t")

  rankBy_factor <- colnames(stitched_regions)[5]
  rankBy_vector <- apply(stitched_regions[5], c(1,2), as.numeric)

  #SETTING NEGATIVE VALUES IN THE rankBy_vector to 0
  rankBy_vector[rankBy_vector < 0] <- 0
  
  #FIGURING OUT THE CUTOFF
  cutoff_options <- calculate_cutoff(rankBy_vector, drawPlot=FALSE,xlab=paste(rankBy_factor,'_enhancers'),ylab=paste(rankByFactor,' Signal','- ',wceName),lwd=2,col=4)

  #These are the super-enhancers
  superEnhancerRows <- which(rankBy_vector> cutoff_options$absolute)
  typicalEnhancers <- setdiff(1:nrow(stitched_regions),superEnhancerRows)
  enhancerDescription <- paste(enhancerName," Enhancers\nCreated from ", enhancerFile,"\nRanked by ",rankBy_factor,"\nUsing cutoff of ",cutoff_options$absolute," for Super-Enhancers",sep="",collapse="")


  #MAKING HOCKEY STICK PLOT
  plotFileName = paste(outFolder,enhancerName,'_Plot_points.png',sep='')
  png(filename=plotFileName, height=18, width=15, unit='in', res=800, pointsize=32)
  signalOrder = order(rankBy_vector,decreasing=TRUE)

  plot(length(rankBy_vector):1,rankBy_vector[signalOrder], col='red',
      xlab=paste(rankBy_factor,'_enhancers'),ylab=paste(rankBy_factor,' Signal'),pch=19,cex=2)	

  abline(h=cutoff_options$absolute,col='grey',lty=2)
  abline(v=length(rankBy_vector)-length(superEnhancerRows),col='grey',lty=2)
  lines(length(rankBy_vector):1,rankBy_vector[signalOrder],lwd=4, col='red')
  text(0,0.8*max(rankBy_vector),paste(' Cutoff used: ',cutoff_options$absolute,'\n','Super-Enhancers identified: ',length(superEnhancerRows))     ,pos=4)

  dev.off()


  #This matrix is just the super_enhancers
  true_super_enhancers <- stitched_regions[superEnhancerRows,]

  additionalTableData <- matrix(data=NA,ncol=2,nrow=nrow(stitched_regions))
  colnames(additionalTableData) <- c("enhancerRank","isSuper")
  additionalTableData[,1] <- nrow(stitched_regions)-rank(rankBy_vector,ties.method="first")+1
  additionalTableData[,2] <- 0
  additionalTableData[superEnhancerRows,2] <- 1


  superTableFile = paste(outFolder,enhancerName,'_SuperEnhancers.table.txt',sep='')
  writeSuperEnhancer_table(true_super_enhancers, enhancerDescription,superTableFile, additionalData= additionalTableData[superEnhancerRows,])

  true_super_enhancers

}
