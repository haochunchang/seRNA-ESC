library(limma)
#setwd("~/research/Mining Key Regulators/script")
#Exprdata <- read.delim2("../data/stitched_HES3_GFP_ESC.tsv")

Normalize <- function(infilepath, outfilepath, method = "quantile"){
    Exprdata <- read.delim2(infilepath)
    row.names(Exprdata) <- Exprdata[, 1] 
    Exprdata <- Exprdata[-c(1), -c(1)]
    Exprdata <- apply(Exprdata, c(1,2), as.numeric)
    
    # Normalize by quantile     
    # Other methods: cyclicloess, scale, ...
    n <- voom(Exprdata, normalize.method = method, plot = TRUE)
    n <- n$E

    # Save plot
    par(mfrow=c(1,1))
    png(filename=paste(outfilepath, "_boxplot", sep=""))
    boxplot(n, main=method)
    title(method)
    dev.off()

    write.csv(n, file=paste(outfilepath, "_normalized.tsv", sep=""), sep="\t") 

} 

