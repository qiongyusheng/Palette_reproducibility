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
rna_co <- readRDS('./RNA/co.rds')
cytof <- readRDS('./cytof/co.rds')
rna <- readRDS('./RNA/rna.rds')
meta <- readRDS('./meta.rds')

gene <- c('PTPRC','PTPRCAP','CD3E','CD3D','CD3G','CD19','IL3RA','ITGAM','CD4','CD8A',
          'CD8B','ITGAX','CD14','FCER1A','NONE',
          'HAVCR2','CD274','CD27','TBX21','CTLA4','FOXP3','CD33',
        'IL7R','CCR7','MKI67','IL2RA','NONE','CD38','HLA-DRA','HLA-DRB1','HLA-DRB3',
        'PDCD1','NCAM1','FCGR3A')
remian_gene <- setdiff(rownames(rna),gene)

rm(co_gene)
rna <- rna[remian_gene,]

rna <- as(as.matrix(rna),'dgCMatrix')
rna_co <- as(as.matrix(rna_co),'dgCMatrix')
cytof <- as(as.matrix(cytof),'dgCMatrix')
rna_input <- rbind(rna,rna_co)

lig <- createLiger(list(rna = rna_input, cytof = cytof))
lig <- normalize(lig)
lig <- selectGenes(lig, useUnsharedDatasets = "rna", unsharedThresh = 0)
lig <- scaleNotCenter(lig)
lig <- runIntegration(lig, k = 20, method = "UINMF")
lig <- quantileNorm(lig, reference = "rna")

uinmf_res <- as.matrix(lig@H.norm)

out_dir <- './'
write.csv(uinmf_res ,
          paste0(out_dir,'UINMF.csv'),
          row.names = FALSE)
