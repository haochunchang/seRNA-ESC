# Brwose many GO results
library(GOSemSim)

main <- function() {
    celltype = "HES3_GFP_ESC"
    state_id = 0
    paths = Sys.glob(file.path('..', 'data', celltype, 'FEA', paste('stage', state_id, sep=''),  
                        '*', 'top20GO_bp.csv'))
    # Rank the seRNA with the most # of mRNA interaction
    #coexpress = auto.get_coexpress_mrna(celltype, state_id)
    hsGO = godata('org.Hs.eg.db', ont="BP")
    
    all_go = list()
    for (p in paths) {
        se = strsplit(p, '/')[[1]][6]
        go = as.vector(read.csv(p)$GO.ID)
        all_go[[se]] = go
    }
   
    m = "Rel"
    b = "BMA"
    goSim = data.frame()
    for (se in names(all_go)) {
        go1 = all_go[[se]]
        for (se2 in names(all_go)) {
            go2 = all_go[[se2]]
            goSim[se, se2] = mgoSim(go1, go2, semData=hsGO, measure=m, combine=b)
        } 
    }
    write.table(goSim, file=paste('goSim_state', state_id, sep=''))
    data = as.matrix(goSim)
    heatmap(data, labRow=rownames(goSim), labCol=colnames(goSim), scale='column')

}
main()
