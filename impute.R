
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
                 k = 5, get.distance = FALSE)
imp1 <- apply(knn1$index, 1, function(idx){
    rowMeans(assay_list[['RNA']][,idx])
})


knn2 <- queryKNN(X = stab$stab[which(stab$mod=='RNA'), ],
                 query = stab$stab[which(stab$mod=='roi2'), ],
                 k = 5, get.distance = FALSE)
imp2 <- apply(knn2$index, 1, function(idx){
    rowMeans(assay_list[['RNA']][,idx])
})

mypal <- brewer.pal(n = 9, name = "Reds")

g <- 'AQP2'
plot(roi2$center_x, roi2$center_y, asp=1, cex = .1, main = g,
     col = colorby(log1p(imp2[g,]), colors = mypal))

g <- 'NPHS2'
plot(roi2$center_x, roi2$center_y, asp=1, pch=19, cex=.25, col = colorby(assay(roi2,'logcounts')[g,], colors = c('white',2,'darkred')),
     main = paste(g,'- Actual'), axes=FALSE); box()

plot(roi2$center_x, roi2$center_y, asp=1, pch=19, cex=.25, col = colorby(imp2[g,], colors = c('white',4,'darkblue')),
     main = paste(g,'- Imputed'), axes=FALSE); box()


layout(matrix(1:12, nrow=3, byrow = TRUE))

par(mar=c(.5,.5,3,.5))
for(g in c('NPHS2','SLC12A1','PTPRQ','MECOM','MEIS1','PTPRO')){
    plot(roi2$center_x, roi2$center_y, asp=1, pch=19, cex=.25, col = colorby(assay(roi2,'logcounts')[g,], colors = c('white',2,'darkred')),
         main = paste(g,'- Actual'), axes=FALSE); box()
    
    plot(roi2$center_x, roi2$center_y, asp=1, pch=19, cex=.25, col = colorby(imp2[g,], colors = c('white',4,'darkblue')),
         main = paste(g,'- Imputed'), axes=FALSE); box()
}

par(mar=c(5,4,4,2)+.1)

