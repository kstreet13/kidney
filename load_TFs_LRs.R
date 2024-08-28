# compile gene lists

# 4 sources:
# - 2000 most variable, 
# - ~250 from spatial, 
# - TFs from paper
# - L-R pair genes from CellChat

# source('setup.R')

# CellChat L-R pairs
load('data/CellChatDB.human.rda')
interxns <- CellChatDB.human$interaction
rm(CellChatDB.human)

# some interactions contain multiple ligrands or multiple receptors
# first column contains the full names for all unique L/R elements
# "receptor" column contains weird mixed names like TGFbR1_R2
LRs <- unique(unlist(strsplit(interxns$interaction_name, split="_")))


# QC
# mean(make.names(unique(interxns$ligand)) %in% rownames(sce))
# mean(make.names(unique(interxns$receptor)) %in% rownames(sce))


interxns$ligand <- make.names(interxns$ligand)
interxns$receptor <- make.names(interxns$receptor)


LRs <- unique(c(interxns$ligand, interxns$receptor))

rm(geneInfo,interxns)


# TFs from Ng et al.
library(readxl)
TFs <- read_excel("data/41587_2020_742_MOESM3_ESM.xlsx")
TFs <- TFs$`Supplementary Table 1: TFs in the Human Tfome`
TFs <- make.names(TFs)
