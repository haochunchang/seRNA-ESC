
#============================================================================
#==============SUPER-ENHANCER CALLING AND PLOTTING FUNCTIONS=================
#============================================================================

 #This function calculates the cutoff by sliding a diagonal line and finding where it is tangential (or as close as possible)
 calculate_cutoff <- function(inputVector, drawPlot=TRUE,...){
 	inputVector <- sort(inputVector)
	inputVector[inputVector<0] <- 0 #set those regions with more control than ranking equal to zero

        #This is the slope of the line we want to slide. This is the diagonal.	
        slope <- (max(inputVector)-min(inputVector))/length(inputVector)
	xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
	y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
	
	if(drawPlot){  #if TRUE, draw the plot
		plot(1:length(inputVector), inputVector,type="l",...)
		b <- y_cutoff-(slope* xPt)
		abline(v= xPt,h= y_cutoff,lty=2,col=8)
		points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
		abline(coef=c(b,slope),col=2)
		title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
		axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
	}
	return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector)))
}

#this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
numPts_below_line <- function(myVector,slope,x){
	yPt <- myVector[x]
	b <- yPt-(slope*x)
	xPts <- 1:length(myVector)
	return(sum(myVector<=(xPts*slope+b)))
}


convert_stitched_to_bed <- function(inputStitched,trackName,trackDescription,outputFile,splitSuper=TRUE,score=c(),superRows=c(),baseColor="0,0,0",superColor="255,0,0"){
	outMatrix <- matrix(data="",ncol=4+ifelse(length(score)==nrow(inputStitched),1,0),nrow=nrow(inputStitched))
	
	outMatrix[,1] <- as.character(inputStitched$CHROM)
	outMatrix[,2] <- as.character(inputStitched$START)
	outMatrix[,3] <- as.character(inputStitched$STOP)
	outMatrix[,4] <- as.character(inputStitched$REGION_ID)
	if(length(score)==nrow(inputStitched)){
		score <- rank(score,ties.method="first")
		score <- length(score)-score+1  #Stupid rank only does smallest to largest. 
		outMatrix[,5] <- as.character(score)
	}
	trackDescription <- paste(trackDescription,"\nCreated on ",format(Sys.time(), "%b %d %Y"),collapse="",sep="")
	trackDescription <- gsub("\n","\t", trackDescription)
	tName <- gsub(" ","_",trackName)
	cat('track name="', tName,'" description="', trackDescription,'" itemRGB=On color=',baseColor,"\n",sep="",file=outputFile)
	write.table(file= outputFile,outMatrix,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
	if(splitSuper==TRUE){
		cat("\ntrack name=\"Super_", tName,'" description="Super ', trackDescription,'" itemRGB=On color=', superColor,"\n",sep="",file=outputFile,append=TRUE)
		write.table(file= outputFile,outMatrix[superRows,],sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
	}
}


convert_stitched_to_gateway_bed <- function(inputStitched,outputFileRoot,splitSuper=TRUE,score=c(),superRows=c()){
	outMatrix <- matrix(data="",ncol=6,nrow=nrow(inputStitched))
	
	outMatrix[,1] <- as.character(inputStitched$CHROM)
	outMatrix[,2] <- as.character(inputStitched$START)
	outMatrix[,3] <- as.character(inputStitched$STOP)
	outMatrix[,4] <- as.character(inputStitched$REGION_ID)
	
	if(length(score)==nrow(inputStitched)){
		score <- rank(score,ties.method="first")
		score <- length(score)-score+1  #Stupid rank only does smallest to largest. 
		outMatrix[,5] <- as.character(score)
	}
	
	outMatrix[,6] <- as.character(rep('.',nrow(outMatrix)))
	
	
	outputFile1 = paste(outputFileRoot,'_Gateway_Enhancers.bed',sep='')
	write.table(file= outputFile1,outMatrix,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
	if(splitSuper==TRUE){
		outputFile2 = paste(outputFileRoot,'_Gateway_SuperEnhancers.bed',sep='')

		write.table(file= outputFile2,outMatrix[superRows,],sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
	}
}


writeSuperEnhancer_table <- function(superEnhancer,description,outputFile,additionalData=NA){
	description <- paste("#",description,"\nCreated on ",format(Sys.time(), "%b %d %Y"),collapse="",sep="")
	description <- gsub("\n","\n#",description)
	cat(description,"\n",file=outputFile)
	if(is.matrix(additionalData)){
		if(nrow(additionalData)!=nrow(superEnhancer)){
			warning("Additional data does not have the same number of rows as the number of super enhancers.\n--->>> ADDITIONAL DATA NOT INCLUDED <<<---\n")
		}else{
			superEnhancer <- cbind(superEnhancer,additionalData)
			superEnhancer = superEnhancer[order(superEnhancer$enhancerRank),]
			
		}
	}
	write.table(file=outputFile,superEnhancer,sep="\t",quote=FALSE,row.names=FALSE,append=TRUE)
}









	
