library(Seurat)
library(dplyr)
library(ggplot2)
library(STACAS)

setwd('./')
data <- readRDS('./dataset.rds')
varg <- readRDS('./vargs.rds')
data <- lapply(data,function(x) x[varg,])
meta <- readRDS('./meta.rds')

object <- CreateSeuratObject(do.call(cbind, data), meta.data = meta, min.cells = 0, min.features = 0)
rm(data)
object_integrated <- object %>% SplitObject(split.by = "Batch") %>%
      Run.STACAS(dims = 1:20, anchor.features = 2000)

res <- Embeddings(object_integrated, reduction = 'pca')
out_dir <- './'
write.csv(as.matrix(res),
          paste0(out_dir,'STACAS.csv'),
          row.names = FALSE)

object_integrated_ss <- object %>% SplitObject(split.by = "Batch") %>%
  Run.STACAS(dims = 1:20, anchor.features = 2000, cell.labels = "CellType")

res <- Embeddings(object_integrated_ss, reduction = 'pca')
write.csv(as.matrix(res),
          paste0(out_dir,'ssSTACAS.csv'),
          row.names = FALSE)
