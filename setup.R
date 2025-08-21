
require(Seurat)
seu <- readRDS('data/MergeOfKidneysv2.rds')
require(SingleCellExperiment)
require(scuttle)
sce <- as.SingleCellExperiment(seu)
rm(seu)
assay(sce,'logcounts') <- log1p(t(t(counts(sce)) / colSums(counts(sce)))*1000)
rownames(sce) <- make.names(rownames(sce)) # spatial gene names will be adjusted like this on import, so it's easier this way
require(scry)
sce <- devianceFeatureSelection(sce, batch = factor(sce$orig.ident))

# plot(reducedDim(sce,'UMAP'), asp=1,
#      col=colorby(assay(sce,'logcounts')['MARCH11',]))


# move spatial setup to its own script
source('setup_spatial.R')
