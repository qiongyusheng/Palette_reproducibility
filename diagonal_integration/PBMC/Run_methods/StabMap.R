library(bindSC)
library(Seurat)
library(data.table)
library(Matrix)
library(matrixStats)
library(StabMap)
library(scran)
library(patchwork)

setwd('./')

rna <- as(as.matrix(readRDS('./RNA/rna.rds')),'dgCMatrix')
rna_co <- as(as.matrix(readRDS('./RNA/co.rds')),'dgCMatrix')
cytof_co <- as.matrix(readRDS('./cytof/co.rds'))

gene <- c('PTPRC','PTPRCAP','CD3E','CD3D','CD3G','CD19','IL3RA','ITGAM','CD4','CD8A',
          'CD8B','ITGAX','CD14','FCER1A','NONE',
          'HAVCR2','CD274','CD27','TBX21','CTLA4','FOXP3','CD33',
        'IL7R','CCR7','MKI67','IL2RA','NONE','CD38','HLA-DRA','HLA-DRB1','HLA-DRB3',
        'PDCD1','NCAM1','FCGR3A')
idx <- which(rownames(rna) %in% gene)
rna <- rna[-idx,]

rna_n <- NormalizeData(rbind(rna,rna_co))
cytof_n <- NormalizeData(cytof_co)

assay_list = list(RNA=rna_n,CyTOF=cytof_n)
mosaicDataUpSet(assay_list, plot = F)
stab = stabMap(assay_list,
               plot = FALSE,
               reference_list = c("RNA"),
               maxFeatures = 10000,
              ncomponentsReference = 20,
              ncomponentsSubset = 20)

out_dir <- './'

write.csv(as.matrix(stab),
          paste0(out_dir,'/stab.csv'),
          row.names = FALSE)
