options(stringsAsFactors=FALSE)
args <- commandArgs(trailingOnly=TRUE)

raw_preprocess <- function(infile, outfile){
 
  inpath <- file.path('..', 'data', 'raw', paste(infile, '.tsv', sep=''))
  print(paste('input file path:', inpath, sep=' '))
  data <- read.table(inpath, sep='\t')
  
  df <- data.frame(sapply(data[-1, -c(1,2)], function(x) { if(is.character(x)) {
        as.numeric(as.character(x))} else {x}}))
  # Remove genes or enhancers which has >75% zeros
  df2 <- df[rowMeans(df!=0)>0.25, ]
  
  # Upper Quartile normalization
  a <- apply(df2, 2, function (x){quantile(x)[[4]]})
  data.normal <- sweep(df2, 2, a, `/`)

  # Log2-transformation
  data.normal <- log2(data.normal+1)

  # Check normalized distribution
  jpeg('uq.jpeg', width=1200, height=1200)
  boxplot(data.normal)
  dev.off()

  # Save preprocessed file
  sample.name <- data[1,]
  locus <- data[,c(1,2)]
  locus <- locus[rowMeans(df!=0)>0.25,]
  new <- cbind(locus, data.normal)
  new_data <- new
  names(new_data) <- as.character(sample.name[1,])
  
  outpath <- file.path('..', 'data', 'raw', paste(outfile, '.tsv', sep=''))
  write.table(new_data, outpath, sep='\t', row.names=FALSE)
  print(paste('saved to ', outpath, sep=' '))

}
raw_preprocess(args[1], args[2])
