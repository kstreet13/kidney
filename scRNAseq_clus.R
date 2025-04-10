##
require(Seurat)
seu <- readRDS('data/MergeOfKidneysv2.rds')


# QC
######################
# VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident')
# keep <- which(seu$nFeature_RNA < 7000,
#               seu$nCount_RNA < 25000,
#               seu$percent.mt <5)
# seu <- seu[,keep]

npcs.seq <- c(seq(10,20,2), seq(25,50,5))
clusmat <- NULL
for(npcs in npcs.seq){
    #seu <- FindNeighbors(seu, reduction = "integrated.cca", dims = 1:npcs)
    #seu <- FindClusters(seu, algorithm = 3, resolution = .5)
    #clusmat <- cbind(clusmat, seu$seurat_clusters)
    
    kmsmall <- kmeans(seu@reductions$integrated.cca@cell.embeddings[,1:npcs], centers = 10, iter.max = 100)$cluster
    kmlarge <- kmeans(seu@reductions$integrated.cca@cell.embeddings[,1:npcs], centers = 25, iter.max = 100)$cluster
    clusmat <- cbind(clusmat, kmsmall)
    clusmat <- cbind(clusmat, kmlarge)
}

agreement1 <- sapply(1:ncol(clusmat), function(i){
    sapply(1:ncol(clusmat), function(j){
        adjustedRandIndex(clusmat[,i], clusmat[,j])
    })
})

dune <- Dune(clusmat, verbose = TRUE, parallel = TRUE)
agreement2 <- sapply(1:ncol(clusmat), function(i){
    sapply(1:ncol(clusmat), function(j){
        adjustedRandIndex(dune$currentMat[,i], dune$currentMat[,j])
    })
})

which.max(colMeans(agreement1)) # 8
which.max(colMeans(agreement2)) # 8

clus <- clusmat[,30]

require(SingleCellExperiment)
require(scuttle)
sce <- as.SingleCellExperiment(seu)
rm(seu)
assay(sce,'logcounts') <- log1p(t(t(counts(sce)) / colSums(counts(sce)))*1000)
rownames(sce) <- make.names(rownames(sce)) # spatial gene names will be adjusted like this on import, so it's easier this way
sce$clus <- clus


require(BiocNeighbors)
knn <- findKNN(reducedDim(sce,'INTEGRATED.CCA')[,1:30], k=10)$index
clusQual <- matrix(sce$clus[as.numeric(knn)], nrow = nrow(knn), ncol = ncol(knn))
clusQual <- rowMeans(clusQual == sce$clus)
sce$clus[clusQual <= .5] <- -1

old <- sum(sce$clus == -1)
converged <- FALSE
while(!converged){
    clusQual <- matrix(sce$clus[as.numeric(knn)], nrow = nrow(knn), ncol = ncol(knn))
    clusQual <- rowMeans(clusQual == sce$clus)
    sce$clus[clusQual <= .5] <- -1
    new <- sum(sce$clus == -1)
    converged <- new == old
    old <- new
}



plot(reducedDim(sce,'UMAP'), asp=1, col = colorby(factor(sce$clus)), pch=c(16,1)[1+(sce$clus==-1)])


layout(matrix(1:25,5,5))
par(mar=c(1,1,1,1))
for(i in unique(sce$clus)){
    if(i != -1){
        plot(range(reducedDim(sce,'UMAP')[,1]),range(reducedDim(sce,'UMAP')[,2]), asp=1, col = 'white', main=i)
        points(reducedDim(sce,'UMAP')[which(sce$clus==i), ], col=colorby(factor(sce$clus))[which(sce$clus==i)])
    }
}
layout(1)
par(mar=c(5,4,4,2)+.1)



boxplot(log1p(sce$nCount_RNA)~sce$clus, col = colorby(factor(sort(unique(sce$clus)))))
boxplot(log1p(sce$nFeature_RNA)~sce$clus, col = colorby(factor(sort(unique(sce$clus)))))


plot(log1p(sce$nCount_RNA), log1p(sce$nFeature_RNA), col = colorby(factor(sce$clus)))


# fill everything we can with "100%" confidence
clusQual <- matrix(sce$clus[as.numeric(knn)], nrow = nrow(knn), ncol = ncol(knn))




