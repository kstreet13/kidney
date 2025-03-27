
source('setup.R')

require(StabMap)
all(rownames(roi1) == rownames(roi2))
rownames(sce)[rownames(sce)=='MARCH11'] <- 'MARCHF11'

# Genes to keep:
# - 2000 most variable, 
# - ~250 from spatial, 
# - TFs from paper
# - L-R pair genes from CellChat
source('load_TFs_LRs.R')

keep <- which(rowData(sce)$binomial_deviance >= sort(rowData(sce)$binomial_deviance, decreasing = TRUE)[2000] |
                  rownames(sce) %in% rownames(roi1) |
                  rownames(sce) %in% TFs |
                  rownames(sce) %in% LRs)

assay_list <- list(RNA = assay(sce,'logcounts')[keep,],
                   roi1 = assay(roi1,'logcounts'),
                   roi2 = assay(roi2,'logcounts'))

stab <- readRDS('data/stabMap_results.rds')



require(BiocNeighbors)
require(Matrix); require(matrixStats)
knn1 <- queryKNN(X = stab$stab[which(stab$mod=='RNA'), ],
                 query = stab$stab[which(stab$mod=='roi1'), ],
                 k = 3, get.distance = FALSE)
imp1 <- apply(knn1$index, 1, function(idx){
    rowMeans(assay_list[['RNA']][,idx])
})



# impute spatial data (get imputed signal for all missing genes)
imp1 <- imputeEmbedding(
    assay_list[c('RNA','roi1')],
    stab$stab[which(stab$mod %in% c('RNA','roi1')),],
    reference = colnames(assay_list[["RNA"]]),
    query = colnames(assay_list$roi1))

imp2 <- imputeEmbedding(
    assay_list[c('RNA','roi2')],
    stab$stab[which(stab$mod %in% c('RNA','roi2')),],
    reference = colnames(assay_list[["RNA"]]),
    query = colnames(assay_list$roi2))



