
require(Seurat)
seu <- readRDS('data/MergeOfKidneysv2.rds')
require(SingleCellExperiment)
sce <- as.SingleCellExperiment(seu)
rm(seu)




plot(reducedDim(sce,'UMAP'), asp=1,
     col=colorby(assay(sce,'logcounts')['MARCH11',]))






# making the spatial data more portable

# coords <- read.csv("~/Projects/kidney/data/spatial/spatialgenomics_lindstrom_roi1_2024-04-20_2103-001/Lindstrom_roi1_13p5wk_BasicResults/roi1_CellCoordinates.csv")
# counts <- read.csv("~/Projects/kidney/data/spatial/spatialgenomics_lindstrom_roi1_2024-04-20_2103-001/Lindstrom_roi1_13p5wk_BasicResults/roi1_CellxGene.csv")
# names(counts)[1] <- 'cell'
# counts <- as.matrix(counts)
# rownames(counts) <- counts[,1]
# counts <- counts[,-1]
# require(Matrix)
# counts <- as(counts, 'dgCMatrix')
# saveRDS(list(coords = coords, counts = counts), file='~/Desktop/roi1.rds')
# 
# coords <- read.csv("~/Projects/kidney/data/spatial/spatialgenomics_lindstrom_roi2_2024-04-20_2116-002/Lindstrom_roi2_16p5wk_BasicResults/roi2_CellCoordinates.csv")
# counts <- read.csv("~/Projects/kidney/data/spatial/spatialgenomics_lindstrom_roi2_2024-04-20_2116-002/Lindstrom_roi2_16p5wk_BasicResults/roi2_CellxGene.csv")
# names(counts)[1] <- 'cell'
# counts <- as.matrix(counts)
# rownames(counts) <- counts[,1]
# counts <- counts[,-1]
# counts <- as(counts, 'dgCMatrix')
# saveRDS(list(coords = coords, counts = counts), file='~/Desktop/roi2.rds')
