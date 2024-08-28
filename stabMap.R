
source('setup.R')

require(StabMap)
# integrate dissociated scRNAseq data with ROI1 spatial data

all(rownames(roi1) == rownames(roi2))

# all genes or just most variable genes?
# check: I think MARCHF11 (spatial) = MARCH11 (RNA). Have same NCBI gene ID
rownames(sce)[rownames(sce)=='MARCH11'] <- 'MARCHF11'

# Genes to keep:
# - 2000 most variable, 
# - ~250 from spatial, 
# - TFs from paper
# - L-R pair genes from CellChat
source('load_TFs_LRs.R')

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
### or ###
mod <- factor(rep(c('RNA','roi1','roi2'),
                  times = c(55297, 192255, 443965)))

umap <- uwot::umap(stab)

ord <- sample(nrow(umap))


layout(matrix(1:4,2,2,byrow = TRUE))
plot(umap[ord,], asp=1, col = alpha(brewer.pal(4,'Set1')[1:3],.4)[mod][ord], main='Combined')

plot(range(umap[,1]),range(umap[,2]), asp=1, col='white', main = 'scRNA-seq')
points(umap[which(mod=='RNA'), ], col=alpha(brewer.pal(9,'Set1')[1],.4))

plot(range(umap[,1]),range(umap[,2]), asp=1, col='white', main = 'Spatial 1')
points(umap[which(mod=='roi1'), ], col=alpha(brewer.pal(9,'Set1')[2],.4))

plot(range(umap[,1]),range(umap[,2]), asp=1, col='white', main = 'Spatial 2')
points(umap[which(mod=='roi2'), ], col=alpha(brewer.pal(9,'Set1')[3],.4))



plot(range(umap[,1]),range(umap[,2]), asp=1, col='white')
points(umap[which(mod=='RNA'), ], col=colorby(assay(sce,'logcounts')['MAFB',]))



# impute spatial data (get imputed signal for all missing genes)
imp <- imputeEmbedding(
    assay_list,
    stab,
    reference = colnames(assay_list[["RNA"]]),
    query = c(colnames(assay_list$roi1),colnames(assay_list$roi2)))




