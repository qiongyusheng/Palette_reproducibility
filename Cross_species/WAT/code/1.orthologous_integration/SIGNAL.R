library(scInt)
library(RSpectra)
source("/home/server/zy/group_scripts/all_methods/utils.R")
library(SIGNAL)

setwd('/home/server/sqy/data/mosaic/input/human_mouse_WAT')
meta <- readRDS('./meta.rds')
human_homo <- as(as.matrix(readRDS('./human/rna_homo.rds')),'dgCMatrix')
mouse_homo <- as(as.matrix(readRDS('./mouse/rna_homo.rds')),'dgCMatrix')
X <- cbind(human_homo,mouse_homo)

X <- Seurat::NormalizeData(X)

res_SIGNAL = Run.gcPCA(X, meta, g_factor = "CellType", b_factor = "Batch", 
                       npcs = 20, lambda = 50,do.scale = F,do.cosine = T)

out_dir <- './'
write.csv(as.matrix(t(res_SIGNAL)),
          paste0(out_dir,'SIGNAL.csv'),
          row.names = FALSE)