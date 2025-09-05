pca <- readRDS('~/Downloads/cheese_pca.rds')
umap <- readRDS("~/Downloads/cheese_umap_annotated_v2.rds")
clus2transfer <- umap$banana_celltype
K <- 15
# for clus2transfer, only the scRNAseq labels matter


# build a vector for identifying the assay each cell came from
assay <- rep('RNA', nrow(pca))
assay[which(umap$og_celltype == 'roi1')] <- 'roi1'
assay[which(umap$og_celltype == 'roi2')] <- 'roi2'

# transfer labels for roi1,roi2
# get KNN for roi1,roi2 in scRNAseq
# transfer most common label (unless there is no plurality)

require(BiocNeighbors)
require(Matrix); require(matrixStats)

# find nearest neighbors
knn <- queryKNN(X = pca[which(assay=='RNA'), ],
                query = pca[which(assay %in% c('roi1','roi2')), ],
                k = K, get.distance = FALSE)
# build vector of transfered labels
tran <- apply(knn$index, 1, function(idx){
    nnLabels <- clus2transfer[idx]
    tab <- sort(table(nnLabels),decreasing = TRUE)
    if(length(tab)==1){
        return(names(tab)[1])
    }
    if(tab[1] == tab[2]){
        return('unclear')
    }else{
        return(names(tab)[1])
    }
})

# construct the final set of cluster labels
umap$transfered <- clus2transfer
umap$transfered[which(assay %in% c('roi1','roi2'))] <- tran






# summarize methods in 6 sentences, with citations (PMID)

