require(SingleCellExperiment)


roi1 <- readRDS('data/spatial/roi1.rds')
roi1 <- SingleCellExperiment(assays = list(counts = t(roi1$counts)), colData = roi1$coords)
assay(roi1,'logcounts') <- log1p(t(t(counts(roi1)) / colSums(counts(roi1)))*1000)
assay(roi1,'logcounts')[is.na(assay(roi1,'logcounts'))] <- 0
colnames(roi1) <- paste0('roi1-', colnames(roi1))

roi2 <- readRDS('data/spatial/roi2.rds')
roi2 <- SingleCellExperiment(assays = list(counts = t(roi2$counts)), colData = roi2$coords)
assay(roi2,'logcounts') <- log1p(t(t(counts(roi2)) / colSums(counts(roi2)))*1000)
assay(roi2,'logcounts')[is.na(assay(roi2,'logcounts'))] <- 0
colnames(roi2) <- paste0('roi2-', colnames(roi2))




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
