library(RSpectra)
library(SIGNAL)

ct_label = "CellType"

setwd('../')
meta <- readRDS('./meta.rds')
data1 <- as(as.matrix(readRDS('./data1/X.rds')),'dgCMatrix')
data2 <- as(as.matrix(readRDS('./data2/X.rds')),'dgCMatrix')
X <- cbind(data1,data2)

X <- Seurat::NormalizeData(X)
res_SIGNAL = Run.gcPCA(X, meta, g_factor = "CellType", b_factor = "Batch", 
                       npcs = 30, lambda = 50,do.scale = F,do.cosine = T)

out_dir <- '../'
write.csv(as.matrix(t(res_SIGNAL)),
          paste0(out_dir,'SIGNAL.csv'),
          row.names = FALSE)