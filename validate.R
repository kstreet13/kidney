require(Seurat)
seu <- readRDS('data/MergeOfKidneysv2.rds')
require(SingleCellExperiment)
require(scuttle)
require(scry)
sce <- as.SingleCellExperiment(seu)
rm(seu)
assay(sce,'logcounts') <- log1p(t(t(counts(sce)) / colSums(counts(sce)))*1000)
rownames(sce) <- make.names(rownames(sce)) # spatial gene names will be adjusted like this on import, so it's easier this way
# require(scry)
sce <- devianceFeatureSelection(sce, batch = factor(sce$orig.ident))
HVGs <- rownames(sce)[which(rowData(sce)$binomial_deviance >= sort(rowData(sce)$binomial_deviance, decreasing = TRUE)[2000])]
# HVGs <- readRDS('data/deviance_top2k.rds')

roi1 <- readRDS('data/roi1.rds')
roi1 <- SingleCellExperiment(assays = list(counts = t(roi1$counts)), colData = roi1$coords)
assay(roi1,'logcounts') <- log1p(t(t(counts(roi1)) / colSums(counts(roi1)))*1000)
assay(roi1,'logcounts')[is.na(assay(roi1,'logcounts'))] <- 0
colnames(roi1) <- paste0('roi1-', colnames(roi1))

###

require(StabMap)
rownames(sce)[rownames(sce)=='MARCH11'] <- 'MARCHF11'
source('load_TFs_LRs.R')
keep <- which(rownames(sce) %in% HVGs |
                  rownames(sce) %in% rownames(roi1) |
                  rownames(sce) %in% TFs |
                  rownames(sce) %in% LRs)

# train = 7:246, test = 1:6
assay_list <- list(RNA = assay(sce,'logcounts')[keep, ],
                   roi1 = assay(roi1,'logcounts')[7:246, ])

rm(sce, TFs, LRs, HVGs, keep)
gc()

# get joint embedding
stab <- stabMap(assay_list, reference_list = c("RNA"),
                suppressMessages = FALSE, maxFeatures = nrow(assay_list$RNA),
                plot = FALSE)


mod <- factor(rep(c('RNA','roi1'),
                  times = c(55297, 192255)))

imp.i <- imputeEmbedding(
    assay_list,
    stab,
    reference = colnames(assay_list[["RNA"]]),
    query = colnames(assay_list$roi1)[1:1000])
imp <- imp.i$RNA[rownames(roi1)[1:6], ]
for(i in 2:192){
    print(i)
    imp.i <- imputeEmbedding(
        assay_list,
        stab,
        reference = colnames(assay_list[["RNA"]]),
        query = colnames(assay_list$roi1)[(1000*i-999):(1000*i)])
    imp <- cbind(imp, imp.i$RNA[rownames(roi1)[1:6], ])
}
imp.i <- imputeEmbedding(
    assay_list,
    stab,
    reference = colnames(assay_list[["RNA"]]),
    query = colnames(assay_list$roi1)[192001:192255])
imp <- cbind(imp, imp.i$RNA[rownames(roi1)[1:6], ])




png(filename = '~/Desktop/impute_validation.png', width = 2000, height = 3500, res = 150)

layout(matrix(1:6, ncol=2))
for(i in 1:6){
    plot(assay(roi1,'logcounts')[rownames(imp)[i], 1:ncol(imp)], imp[i,], asp=1, col = rgb(0,0,0,.5), main = rownames(imp)[i])
    abline(0,1)
}

dev.off()



real <- assay(roi1,'logcounts')[rownames(imp), 1:ncol(imp)]

rowMeans(real)
rowMeans(imp)

cor(real[1,],imp[1,])
cor(real[2,],imp[2,])
cor(real[3,],imp[3,])
cor(real[4,],imp[4,])
cor(real[5,],imp[5,])
cor(real[6,],imp[6,])





