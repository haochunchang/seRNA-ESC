library(parallel)
library(GOSemSim)

# Calculate pairwise GOSemSim of a single gene set
# Return a single value descibing how functional compact the set is
# Based on GO term semantic similarity
compactness <- function(single_sym, hsGO, m) {

    sim = mgeneSim(single_sym, semData=hsGO, measure=m, verbose=FALSE, drop=NULL)
    # do not consider itself
    if (is.matrix(sim)) {
        sim[upper.tri(sim, diag=TRUE)] = 0
        rowmax = apply(sim, MARGIN=1, FUN=max)
        semsim = mean(rowmax, na.rm=TRUE)
        valid = nrow(sim)
    } else {
        semsim = mean(sim, na.rm=TRUE)
        valid = 1
    }
    return(list('semsim'=semsim, 'valid'=valid))
}

# Calculate each seRNA's functional consistency in terms of their co-expressed mRNA GO semantic similarity
# Return a mapping of seRNA: p-value
gosim <- function(mrna, hsGO, sym_pool) {
   
    mrna = unlist(mrna)
    # extract GO terms of every single mrna
    single_sym = sym_pool[sym_pool %in% mrna]
    # calculate GOsemsim and combine into a single value `Sn`
    m = "Rel"
    
    sn = compactness(single_sym, hsGO, m)
    n_sample = sn$valid
    # Random sample k genes and calculate their GOsemsim
    # Check if random_gosemsim > `Sn`, if yes, count++
    k = 500
    count = 1
    random_score = c(rep.int(0, k))
    while (count <= k) {
        count = count + 1
        tryCatch({
                random_idx = sample(1:length(sym_pool), n_sample, replace=FALSE)    
                random_score[count] = compactness(sym_pool[random_idx], hsGO, m)$semsim
            }, error=function(e) {
                count = count - 1
            }, finally={
                next
            }
        )
    }
    random_score[random_score <= sn$semsim] = 0
    random_score[random_score > sn$semsim] = 1

    # p-value = count of random_gosemsim > `Sn` / k
    pv = mean(random_score, na.rm=TRUE)
    return(list("consist"=sn$semsim, "p_val"=pv))
}

only_hgnc <- function(sym_pool) {
    library('org.Hs.eg.db')
    sym_pool = t(sym_pool)[,1]
    sympool = select(org.Hs.eg.db, sym_pool, c('GENENAME'), "ALIAS")["ALIAS"]
    colnames(sympool) = c('symbol')
    sympool = na.omit(sympool[!(sympool == "NA")])
    return(unique(sympool))
}
