library(bindSC)
library(Seurat)
library(data.table)
library(Matrix)
library(matrixStats)
library(Signac)
library(StabMap)
library(scran)
library(patchwork)

setwd('./')
rna <- readRDS('./RNA/X.rds')
atac <- readRDS('./ATAC/atac.rds')
atac_co <- readRDS('./ATAC/co.rds')

rna_n <- NormalizeData(rna)
atac_n <- RunTFIDF(atac)
atac_co_n <- NormalizeData(atac_co)

rm(rna,atac,atac_co)
gc()

assay_list = list(RNA=rna_n,ATAC=rbind(atac_n,atac_co_n))
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("ATAC"),
               maxFeatures = 10000,
              ncomponentsReference = 20,
              ncomponentsSubset = 20)

meta <- readRDS('./meta.rds')

stab <- stab[rownames(meta),]

data.umap <- uwot::umap(as.matrix(stab),min_dist = 0.3)
dim.plot0(data.umap,meta,"modality")
dim.plot0(data.umap,meta,"label")

out_dir <- './kidney/'

write.csv(as.matrix(stab),
          paste0(out_dir,'/stab.csv'),
          row.names = FALSE)
