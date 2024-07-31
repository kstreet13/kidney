
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

roi1 <- readRDS('data/roi1.rds')
roi1 <- SingleCellExperiment(assays = list(counts = t(roi1$counts)), colData = roi1$coords)
assay(roi1,'logcounts') <- log1p(t(t(counts(roi1)) / colSums(counts(roi1)))*1000)
assay(roi1,'logcounts')[is.na(assay(roi1,'logcounts'))] <- 0
roi2 <- readRDS('data/roi2.rds')
roi2 <- SingleCellExperiment(assays = list(counts = t(roi2$counts)), colData = roi2$coords)
assay(roi2,'logcounts') <- log1p(t(t(counts(roi2)) / colSums(counts(roi2)))*1000)
assay(roi2,'logcounts')[is.na(assay(roi2,'logcounts'))] <- 0

