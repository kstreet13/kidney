require(StabMap)


# integrate dissociated scRNAseq data with ROI1 spatial data

# all genes or just most variable genes?
# check: I think MARCHF11 (spatial) = MARCH11 (RNA). Have same NCBI gene ID
rownames(sce)[rownames(sce)=='MARCH11'] <- 'MARCHF11'

keep <- which(rowData(sce)$binomial_deviance >= sort(rowData(sce)$binomial_deviance, decreasing = TRUE)[2000] |
                  rownames(sce) %in% rownames(roi1))

assay_list <- list(RNA = assay(sce,'logcounts')[keep,],
                   spatial = assay(roi1,'logcounts'))

mosaicDataUpSet(assay_list, plot = FALSE)

sm <- stabMap(assay_list, reference_list = c("RNA"),
              suppressMessages = FALSE, maxFeatures = 2000,
              plot = FALSE)

mod <- factor(rep(c('RNA','spatial'), times=c(ncol(assay_list$RNA),ncol(assay_list$spatial))))

umap <- uwot::umap(sm)

ord <- sample(nrow(umap))

plot(umap[ord,], asp=1, col = alpha(brewer.pal(4,'Set1')[1:2],.4)[mod][ord])

plot(range(umap[,1]),range(umap[,2]), asp=1, col='white')
points(umap[which(mod=='RNA'), ], col=alpha(brewer.pal(9,'Set1')[1],.4))

plot(range(umap[,1]),range(umap[,2]), asp=1, col='white')
points(umap[which(mod=='spatial'), ], col=alpha(brewer.pal(9,'Set1')[2],.4))
