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
# codex <- readRDS('./CODEX/codex.rds')
rna_co <- readRDS('./RNA/co_all.rds')
codex_co <- readRDS('./CODEX/co_all.rds')
co_feat <- readRDS('./co_feat.rds')

remian_gene <- setdiff(rownames(rna),co_feat[,1])

rna <- rna[remian_gene,]
rna <- as(as.matrix(rna),'dgCMatrix')
rna_co <- as(as.matrix(rna_co),'dgCMatrix')
codex_co <- as(as.matrix(codex_co),'dgCMatrix')
rna_input <- rbind(rna,rna_co)
lig <- createLiger(list(rna = rna_input, cytof = codex_co))
lig <- normalize(lig)
lig@datasets$cytof@normData <- lig@datasets$cytof@rawData
lig <- selectGenes(lig, useUnsharedDatasets = "rna", unsharedThresh = 0)
lig <- scaleNotCenter(lig)
lig <- runIntegration(lig, k = 20, method = "UINMF")
lig <- quantileNorm(lig, reference = "rna")
uinmf_res <- as.matrix(lig@H.norm)
out_dir <- './'
write.csv(uinmf_res ,
          paste0(out_dir,'UINMF.csv'),
          row.names = FALSE)