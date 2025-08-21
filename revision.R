require(Seurat)
require(SingleCellExperiment)
require(scry)
require(StabMap)


#######################
# SET UP RNA-SEQ DATA #
#######################

scDataset <- "cuttlefish"
# scDataset <- "cheese"

if(scDataset == "cuttlefish"){
    sce <- readRDS('data/scRNAseq_rds/Cuttlefish_kidmerge_NE_UE4_allNPC_RV_UE_dp.rds')
    sce <- as.SingleCellExperiment(sce)
    sce <- devianceFeatureSelection(sce, batch = factor(sce$orig.ident))
}



#######################
# SET UP SPATIAL DATA #
#######################

source('setup_spatial.R')

# check
all(rownames(roi1) == rownames(roi2))

# fix some names
rownames(sce)[rownames(sce)=='MARCH11'] <- 'MARCHF11'
rownames(sce)[rownames(sce)=='TMEM246-AS1'] <- 'TMEM246.AS1'
rownames(sce)[rownames(sce)=='MSC-AS1'] <- 'MSC.AS1'
rownames(sce)[rownames(sce)=='HNF4A-AS1'] <- 'HNF4A.AS1'
rownames(sce)[rownames(sce)=='CLDN10-AS1'] <- 'CLDN10.AS1'

##################
# GENE SELECTION #
##################

# Genes to keep:
# - 2000 most variable, 
# - ~250 from spatial, 
# - TFs from paper
# - L-R pair genes from CellChat

# CellChat L-R pairs
LRs <- readRDS('data/LRs.rds')

# TFs from Ng et al.
TFs <- readRDS('data/TFs.rds')


keep <- which(rowData(sce)$binomial_deviance >= sort(rowData(sce)$binomial_deviance, decreasing = TRUE)[2000] |
                  rownames(sce) %in% rownames(roi1) |
                  rownames(sce) %in% TFs |
                  rownames(sce) %in% LRs)

assay_list <- list(RNA = assay(sce,'logcounts')[keep,],
                   roi1 = assay(roi1,'logcounts'),
                   roi2 = assay(roi2,'logcounts'))

mosaicDataUpSet(assay_list, plot = FALSE)

# get joint embedding
stab <- stabMap(assay_list, reference_list = c("RNA"),
                suppressMessages = FALSE, maxFeatures = nrow(assay_list$RNA),
                plot = FALSE)

mod <- factor(rep(c('RNA','roi1','roi2'),
                  times = c(ncol(assay_list$RNA), ncol(assay_list$roi1), ncol(assay_list$roi2))))








rownames(roi1)[!rownames(roi1) %in% rownames(sce)]

rownames(sce)[grep('HNF4A', rownames(sce))]




