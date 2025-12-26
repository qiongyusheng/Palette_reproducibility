library(Seurat)
library(dplyr)
library(ggplot2)
library(STACAS)

setwd('./')
data <- readRDS('./dataset.rds')
varg <- readRDS('./vargs.rds')
data <- lapply(data,function(x) x[varg,])
meta <- readRDS("./meta.rds")

task <- c('label_miss5','label_miss10','label_miss20')
mod <- c('_miss5','_miss10','_miss20')

object <- CreateSeuratObject(do.call(cbind, data), meta.data = meta, min.cells = 0, min.features = 0) %>%  NormalizeData()
obj.list <- SplitObject(object,split.by = "Batch")
rm(object)
for(i in 1:3){
  object_integrated_ss <- obj.list %>%
    Run.STACAS(dims = 1:20, anchor.features = 4000, cell.labels = task[i])
  res <- Embeddings(object_integrated_ss, reduction = 'pca')
  out_dir <- './label_miss/'
  write.csv(as.matrix(res),
            paste0(out_dir,'ssSTACAS',mod[i],'.csv'),
            row.names = FALSE)
}

