library(batchelor)
library(harmony)
library(Seurat)
library(scater)
library(irlba)
library(BiocNeighbors)
library(SingleCellExperiment)
library(Matrix)
library(umap)
library(dplyr)
library(Rcpp)
library(dplyr)
library(cowplot)
library(patchwork)
library(rliger)

setwd('./')
rna <- readRDS('./RNA/X.rds')
atac <- readRDS('./ATAC/atac.rds')
atac_co <- readRDS('./ATAC/co.rds')

lig <- createLiger(list(rna = as(as.matrix(rna),'dgCMatrix'), atac = as(as.matrix(atac_co),'dgCMatrix')))
lig <- normalize(lig) %>%
    selectGenes(useDatasets = "rna", thresh = 0) %>%
    scaleNotCenter()

unshareNormed <- normalize(as(as.matrix(atac),'dgCMatrix'))
se <- CreateSeuratObject(unshareNormed)
se <- FindVariableFeatures(se, selection.method = "vst", nfeatures = 2000)
top2000 <- VariableFeatures(se)
unshareScaled <- scaleNotCenter(unshareNormed[top2000,])

varUnsharedFeatures(lig, "atac") <- top2000
scaleUnsharedData(lig, "atac") <- unshareScaled

rm(rna,atac,atac_co,unshareNormed,se,unshareScaled)
gc()

lig <- runUINMF(lig, k = 20, nIteration = 30)
lig <- quantileNorm(lig)
uinmf_res <- as.matrix(lig@H.norm)

out_dir <- './kidney/'
write.csv(uinmf_res ,
          paste0(out_dir,'UINMF.csv'),
          row.names = FALSE)




