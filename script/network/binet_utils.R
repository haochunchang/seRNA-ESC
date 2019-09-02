# Read in the R libraries 
options(warn=-1)
library(MASS)
library(class)   
library(cluster) 
library(impute)   
suppressMessages(library(WGCNA))
library(limma)
library(psych)
#=============================================================
# Setup environment variable for developing
#celltype <- 'HES3_GFP_ESC_UQ'
#stage_id <- 0
#datExpr1 <- read.csv(file.path('..', 'data', celltype, 'selected_seRNA_state0.csv'))
#datExpr1$locus <- NULL
#datExpr2 <- read.csv(file.path('../data', celltype, paste('genes_', celltype, '_avg.tsv', sep='')), sep='\t')
#datExpr2 <- datExpr2[,-1]
# ===========================================================
# The function PickSoftThreshold allows one to estimate the power parameter when using
# a soft thresholding approach with the use of the power function AF(s)=s^Power
# The removeFirst option removes the first point (k=1, P(k=1)) from the regression fit.
# Assume datExpr is (genes, seRNA) correlation dataframe 
if (exists("myPickSoftThreshold")) rm(myPickSoftThreshold);
myPickSoftThreshold <- function(datExpr1, datExpr2, method, RsquaredCut=0.85, powervector=c(seq(1,10,by=1),seq(12,20,by=2)),
removeFirst=FALSE,no.breaks=10) {
  
  datExpr1 = t(datExpr1)
  datExpr2 = t(datExpr2)

  colname1=c("Power", "scale law R^2" ,"slope", "truncated R^2","mean(k)","median(k)","max(k)")
  datout=data.frame(matrix(666,nrow=length(powervector),ncol=length(colname1) ))
  names(datout)=colname1
  datout[,1]=powervector
  
  if(exists("fun1")) rm(fun1)
  fun1=function(x) {
    if (method=='bicor'){ corx=abs(bicor(x, datExpr2, use='all.obs')) }
    else if (method=='pearson'){ corx=abs(cor(x, datExpr2, use='all.obs', method='pearson')) }
    else if (method=='spearman'){ corx=abs(cor(x, datExpr2, use='all.obs', method='spearman'))}
    out1=rep(NA, length(powervector) )
    for (j in c(1:length(powervector))) {out1[j]=sum(corx^powervector[j])}
    out1
  } # end of fun1
  datk=t(apply(datExpr1,2,fun1))
  
  for (i in c(1:length(powervector) ) ){
    nolinkshelp <- unique(datk[,i]-1)
    cut2=cut(nolinkshelp,no.breaks)
    binned.k=tapply(nolinkshelp,cut2,mean)
    freq1=as.vector(tapply(nolinkshelp,cut2,length)/length(nolinkshelp))
  # The following code corrects for missing values etc
    breaks1=seq(from=min(nolinkshelp, na.rm=TRUE),to=max(nolinkshelp, na.rm=TRUE),length=no.breaks+1)
    hist1=hist(nolinkshelp,breaks=breaks1,equidist=F,plot=FALSE,right=TRUE)
    binned.k2=hist1$mids
    binned.k=ifelse(is.na(binned.k),binned.k2,binned.k)
    binned.k=ifelse(binned.k==0,binned.k2,binned.k)
    freq1=ifelse(is.na(freq1),0,freq1)

    xx= as.vector(log10(binned.k))
    if(removeFirst) {freq1=freq1[-1]; xx=xx[-1]}
    plot(xx,log10(freq1+.000000001),xlab="log10(k)",ylab="log10(p(k))" )
    lm1= lm(as.numeric(log10(freq1+.000000001))~ xx )
    lm2=lm(as.numeric(log10(freq1+.000000001))~ xx+I(10^xx) )
    datout[i,2]=summary(lm1)$adj.r.squared 
    datout[i,3]=summary(lm1)$coefficients[2,1]  
    datout[i,4]=summary(lm2)$adj.r.squared
    datout[i,5]=mean(nolinkshelp, na.rm=TRUE)
    datout[i,6]= median(nolinkshelp, na.rm=TRUE)
    datout[i,7]= max(nolinkshelp, na.rm=TRUE) 
  } 
  datout=signif(datout,3)
  # the cut-off is chosen as smallest cut with R^2 > RsquaredCut 
  ind1=datout[,2]>RsquaredCut
  indcut=NA
  indcut=ifelse(sum(ind1)>0,min(c(1:length(ind1))[ind1]),indcut)
  # this estimates the power value that should be used. 
  # Don't trust it. You need to consider slope and mean connectivity as well!
  power.estimate=powervector[indcut][[1]]
  list(power.estimate, data.frame(datout));
}


#==========================================================
if(exists('choose_power')) rm(choose_power);
choose_power <- function(Exprdata1, Exprdata2, celltype, stage_id, method='pearson', powers=c(seq(1,10,by=1),seq(12,20,by=2))){

  #Investigate soft thesholding with the power adjacency function 
  RpowerTable = myPickSoftThreshold(Exprdata1, Exprdata2, method, powervector=powers)
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
# Calculate correlation matrix
if(exists('calculate_cor')) rm(calculate_cor);
calculate_cor <- function(datExpr1, datExpr2, method='pearson', abs_flag=TRUE){
 
  datExpr1 = t(datExpr1)
  datExpr2 = t(datExpr2)
  
  if (abs_flag==TRUE) {
    if (method=='bicor'){
      corx=abs(bicor(datExpr1, datExpr2, use='all.obs'))
    } else if (method=='pearson'){
      corx=abs(cor(datExpr1, datExpr2, use='all.obs', method='pearson'))
    } else if (method=='spearman'){
      corx=abs(cor(datExpr1, datExpr2, use='all.obs', method='spearman'))
    }
  } else {
    if (method=='bicor'){
      corx=bicor(datExpr1, datExpr2, use='all.obs')
    } else if (method=='pearson'){
      corx=cor(datExpr1, datExpr2, use='all.obs', method='pearson')
    } else if (method=='spearman'){
      corx=cor(datExpr1, datExpr2, use='all.obs', method='spearman')
    }
  }
   
  return(data.frame(corx))
}
#====================================================================================
# Calculate correlation matrix
if(exists('calculate_cor_test')) rm(calculate_cor_test);
calculate_cor_test <- function(datExpr1, datExpr2, p_value=0.05, method='pearson'){
  
  datExpr1 = t(datExpr1)
  datExpr2 = t(datExpr2)

  corx = corr.test(datExpr1, datExpr2, method=method, alpha=p_value, adjust='none', ci=FALSE)
  corr = corx$r
  corr[corx$p > p_value] = 0
  corr = abs(corr)
  return(data.frame(corr))
}

#====================================================================================
# Create a scale free topology plot.
# The black curve corresponds to scale free topology
# Input: datcorr is pre-computed correlation matrix
if(exists('check_scalefree')) rm(check_scalefree);
check_scalefree <- function(datcorr, beta, celltype, stage_id, truncated_flag=F){
  
  Connectivity=colSums(datcorr^beta)-1
 
  dir.create(file.path('../graph', celltype), showWarnings=FALSE) 
  jpeg(filename=file.path('../graph', celltype, paste('stage ', stage_id, ' check_scalefree.jpg', sep='')), width=800, height=600)
  par(mfrow=c(1,1))
  scaleFreePlot(Connectivity, nBreaks=30, main=paste("soft threshold, power=",beta), truncated=truncated_flag)
  dev.off()
}

#==============================================================================
# Utility function 
# Read in gene exp profiles of several stages in a R list.
# Return an R list containing several data.frame
if (exists('Load_gene_stages')) rm(Load_gene_stages);
Load_gene_stages <- function(celltype, nstages){
  
  genes = list() 
  for (j in seq(0, nstages-1)){
    gene = read.table(file.path('..', 'data', celltype, paste('stage', j, '_genes.tsv', sep=''))
                      , sep='\t', stringsAsFactors=FALSE)
    gene_numeric = apply(gene[-1, -1], c(1,2), as.numeric)
    rownames(gene_numeric) = gene[,1][-1]
    colnames(gene_numeric) = gene[1,][-1]
    genes[[j+1]] = data.frame(t(gene_numeric))
  }
  return(genes)
}
