library(Seurat)
library(data.table)
library(Matrix)
library(matrixStats)
setwd('./')

meta <- readRDS('./Human_kidney/meta.rds')

simba <- fread('./simba.csv',data.table = F,header = T)[,1:20]
cell_name <- readRDS('./simba_cell_name.rds')
rownames(simba) <- cell_name[,1]

rownames(meta)[1:19985] <- paste0(meta$barcode[1:19985],'_',meta$modality[1:19985])
simba <- simba[rownames(meta),]
write.csv(as.matrix(simba),
         './SIMBA.csv',
          row.names = FALSE)