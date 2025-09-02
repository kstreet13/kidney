require(Seurat)
require(SingleCellExperiment)
require(scry)
require(StabMap)


#######################
# SET UP RNA-SEQ DATA #
#######################

# scDataset <- "cuttlefish"
scDataset <- "cheese"

if(scDataset == "cuttlefish"){
    sce <- readRDS('data/scRNAseq_rds/Cuttlefish_kidmerge_NE_UE4_allNPC_RV_UE_dp.rds')
    sce <- as.SingleCellExperiment(sce)
    sce <- devianceFeatureSelection(sce, batch = factor(sce$orig.ident))
}
if(scDataset == 'cheese'){
    sce <- readRDS('data/scRNAseq_rds/Kim_2024_Dev_Cell_wk15_17_sub-001.rds')
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
if(scDataset == 'cuttlefish'){
    rownames(sce)[rownames(sce)=='MARCH11'] <- 'MARCHF11'
    rownames(sce)[rownames(sce)=='TMEM246-AS1'] <- 'TMEM246.AS1'
    rownames(sce)[rownames(sce)=='MSC-AS1'] <- 'MSC.AS1'
    rownames(sce)[rownames(sce)=='HNF4A-AS1'] <- 'HNF4A.AS1'
    rownames(sce)[rownames(sce)=='CLDN10-AS1'] <- 'CLDN10.AS1'
    # probably wrong:
    # rownames(sce)[rownames(sce)=='HNF4A'] <- 'HNF4A.AS1'
    # rownames(sce)[rownames(sce)=='MIR3142HG'] <- 'MIR31HG'
    # rownames(sce)[rownames(sce)=='THEMIS2'] <- 'THEMIS'
    # couldn't find match in RNA:  'INSYN2B' 'LINC00210' 'LINC02798' 'LINC02147' 'XCR1' 'FENDRR'
}

if(scDataset == 'cheese'){
    rownames(sce)[rownames(sce)=='MARCH11'] <- 'MARCHF11'
    rownames(sce)[rownames(sce)=='TMEM246-AS1'] <- 'TMEM246.AS1'
    rownames(sce)[rownames(sce)=='MSC-AS1'] <- 'MSC.AS1'
    rownames(sce)[rownames(sce)=='HNF4A-AS1'] <- 'HNF4A.AS1'
    rownames(sce)[rownames(sce)=='CLDN10-AS1'] <- 'CLDN10.AS1'
}


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

###############
# RUN STABMAP #
###############
# get joint embedding (not using RNA as reference, because it contains fewer cell types)
stab <- stabMap(assay_list, #reference_list = c("RNA"),
                suppressMessages = FALSE, maxFeatures = nrow(assay_list$RNA),
                plot = FALSE)

mod <- factor(rep(c('RNA','roi1','roi2'),
                  times = c(ncol(assay_list$RNA), ncol(assay_list$roi1), ncol(assay_list$roi2))))

# follow-up PCA
require(BiocSingular)
pca <- BiocSingular::runPCA(stab, rank=50)

# UMAP
umap <- uwot::umap(pca$x)

# stab <- list(stab = stab, mod = mod, umap = umap)
# saveRDS(stab, file='data/stabMap_results.rds')

ord <- sample(nrow(umap))

layout(matrix(1:4,2,2,byrow = TRUE))
plot(umap[ord,], asp=1, col = alpha(brewer.pal(4,'Set1')[1:3],.4)[mod][ord], main='Combined', pch=16, cex=.5)

plot(range(umap[,1]),range(umap[,2]), asp=1, col='white', main = 'scRNA-seq')
points(umap[which(mod=='RNA'), ], col=alpha(brewer.pal(9,'Set1')[1],.4), cex=.5)

plot(range(umap[,1]),range(umap[,2]), asp=1, col='white', main = 'Spatial 1')
points(umap[which(mod=='roi1'), ], col=alpha(brewer.pal(9,'Set1')[2],.4), cex=.5)

plot(range(umap[,1]),range(umap[,2]), asp=1, col='white', main = 'Spatial 2')
points(umap[which(mod=='roi2'), ], col=alpha(brewer.pal(9,'Set1')[3],.4), cex=.5)




# specific gene
plot(range(umap[,1]),range(umap[,2]), asp=1, col='white')
points(umap[which(mod=='RNA'), ], col=colorby(assay(sce,'logcounts')['MAFB',]))



# saving
saveRDS(pca$x, file='~/Desktop/cheese_pca.rds')
saveRDS(umap, file='~/Desktop/cheese_umap.rds')



# for fixing names

rownames(roi1)[!rownames(roi1) %in% rownames(sce)]

rownames(sce)[grep('MARCH', rownames(sce))]




