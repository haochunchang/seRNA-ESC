source("http://bioconductor.org/biocLite.R")

# Check packages from Bioconductor
bio_package.list <- c("AnnotationDbi", "GO.db", "preprocessCore", "topGO", 
                      "GOSemSim", "org.Hs.eg.db", "biomaRt", "impute", 
		      "biomartr", "limma", "Rgraphviz", "pathview", "RGSEA")

for (i in 1:length(bio_package.list)){
  if (require(bio_package.list[i], character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)){
    print(paste('Requirement: ', bio_package.list[i], ' is satisfied.', sep=''))
  }
  else {
    biocLite(bio_package.list[i])  
  }
}

# Check packages from CRAN
cran_package.list <- c("MASS", "class", "cluster", "data.table", "psych", 
                        "doParallel", "WGCNA", "dplyr", "parallel",
                      "pathfindR")

for (i in 1:length(cran_package.list)){
  if (require(cran_package.list[i], character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)){
    print(paste('Requirement: ', cran_package.list[i], ' is satisfied.', sep=''))
  }
  else {
    install.packages(cran_package.list[i], repos="http://cran.csie.ntu.edu.tw/")
  }
}

quit(save='no')
