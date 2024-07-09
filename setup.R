
require(Seurat)
seu <- readRDS('data/MergeOfKidneysv2.rds')
require(SingleCellExperiment)
sce <- as.SingleCellExperiment(seu)
rm(seu)




plot(reducedDim(sce,'UMAP'), asp=1,
     col=colorby(assay(sce,'logcounts')['MARCH11',]))


