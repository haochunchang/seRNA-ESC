
convert_ensembl_to_hgnc <- function(enst) {  
    library(biomaRt)
    library(GO.db)
    ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", dataset="hsapiens_gene_ensembl")
    sym = getBM(attributes=c('ensembl_transcript_id', 'hgnc_symbol'), 
                        filters='ensembl_transcript_id', values=enst, mart=ensembl) 
    return(sym)
}
