library(bindSC)
library(Seurat)
library(data.table)
library(Matrix)
library(matrixStats)
library(Signac)
library(StabMap)
library(scran)
library(patchwork)

setwd('./MERFISH_ATAC')

rna <- readRDS('./MERFISH/merfish.rds')
atac <- readRDS('./ATAC/LSI.rds')
atac_co <- readRDS('./ATAC/rna.rds')
top200 <- readRDS('./RF200.rds')

rna_n <- as.matrix(rna)
atac_co_n <- NormalizeData(atac_co[top200[,1],])

rm(rna,atac_co)
gc()

assay_list = list(RNA=as.matrix(rna_n),ATAC=rbind(t(atac),atac_co_n))
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("RNA"),
               maxFeatures = 10000,
              ncomponentsReference = 20,
              ncomponentsSubset = 20)

meta <- readRDS('./meta.rds')
stab <- stab[rownames(meta),]

out_dir <- './Reduced_feat/RF200/'

write.csv(as.matrix(stab),
          paste0(out_dir,'/stab.csv'),
          row.names = FALSE)


