average_super <- function(data, outfilepath='../../data/HES3_GFP_ESC/super-enhancer_avg_exp.tsv'){  
  # Preprocess
  Exprdata <- data[-c(1), -c(1)]
  row.names(Exprdata) <- data[-c(1),1]
  Exprdata <- apply(Exprdata, c(1,2), as.numeric)
  Exprdata <- t(Exprdata)
  rm(data)

  avgExprdata <- data.frame(t(Exprdata[1,]))
  for (i in seq(1, 21, 3)){
    average <- t((Exprdata[i,] + Exprdata[i+1,] + Exprdata[i+2,]) / 3)
    avgExprdata[i,] <- average
  }
  # Day 7 only have 2 replicates
  i <- 22
  average <- t((Exprdata[i,] + Exprdata[i+1,]) / 2)
  avgExprdata[i,] <- data.frame(average)
  # Day 8 ~ Day 12
  for (i in seq(24, 38, 3)){
    average <- t((Exprdata[i,] + Exprdata[i+1,] + Exprdata[i+2,]) / 3)
    avgExprdata[i,] <- data.frame(average)
  }
  
  avgExprdata <- na.omit(avgExprdata)
  row.names(avgExprdata) <- row.names(Exprdata)[c(seq(1,21,3), 22, seq(24,38,3))]
  
  # Save averaged data
  write.table(avgExprdata, outfilepath, sep = "\t")
  return(avgExprdata)
}

# Load in super-enhancer data
data <- read.csv('../../data/HES3_GFP_ESC/super-enhancer_exp.tsv', sep = "\t")
average_super(data) 
test <- read.csv('../../data/HES3_GFP_ESC/super-enhancer_avg_exp.tsv', sep = "\t")
View(test)
