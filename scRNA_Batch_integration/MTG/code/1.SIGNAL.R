library(scInt)
library(RSpectra)
library(SIGNAL)

setwd('./')
data <- readRDS('./dataset.rds')
varg <- readRDS('./vargs.rds')
data <- lapply(data,function(x) x[varg,])
meta <- readRDS('./meta.rds')
X = do.call(cbind,data)
X <- Seurat::NormalizeData(X)               

res_SIGNAL = Run.gcPCA(X, meta, g_factor = "SubClass", b_factor = "Batch", 
                       npcs = 20, lambda = 50,do.scale = F,do.cosine = T)

out_dir <- './'
write.csv(as.matrix(t(res_SIGNAL)),
          paste0(out_dir,'SIGNAL.csv'),
          row.names = FALSE)

