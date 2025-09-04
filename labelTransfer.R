pca <- readRDS('~/Downloads/cheese_pca.rds')
umap <- readRDS("~/Downloads/cheese_umap_annotated_v2.rds")
clus2transfer <- umap$banana_celltype
# for clus2transfer, only the scRNAseq labels matter


# build a vector for identifying the assay each cell came from
assay <- rep('RNA', nrow(pca))
assay[which(umap$og_celltype == 'roi1')] <- 'roi1'
assay[which(umap$og_celltype == 'roi2')] <- 'roi2'

# transfer labels for roi1,roi2
# get KNN (k=3) for roi1,roi2 in scRNAseq
# if 2 NNs in scRNAseq share label, transfer it

require(BiocNeighbors)
require(Matrix); require(matrixStats)

# find nearest neighbors
knn <- queryKNN(X = pca[which(assay=='RNA'), ],
                query = pca[which(assay %in% c('roi1','roi2')), ],
                k = 3, get.distance = FALSE)
# build vector of transfered labels
tran <- apply(knn$index, 1, function(idx){
    nnLabels <- clus2transfer[idx]
    if(any(duplicated(nnLabels))){
        return(nnLabels[duplicated(nnLabels)][1])
    }else{
        return('unclear')
    }
})

# construct the final set of cluster labels
umap$transfered <- clus2transfer
umap$transfered[which(assay %in% c('roi1','roi2'))] <- tran






# summarize methods in 6 sentences, with citations (PMID)

