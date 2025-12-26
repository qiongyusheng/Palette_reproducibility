library(Seurat)
library(cowplot)
library(dplyr)
library(data.table)
library(tibble)
library(Matrix)

setwd('./')
data <- readRDS('./dataset.rds')
varg <- readRDS('./vargs.rds')
meta <- readRDS('./meta.rds')
data <- lapply(data,function(x) x[varg,])

source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
options(future.globals.maxSize = 500 * 1024^3)
library(Matrix)

names(data) <- NULL
musa_obj <- Create.Musa.Object(data.list = data, 
                               samples = paste0('B',1:9), 
                               modals = rep('rna',9),
                               filter=T,
                               min.cells = c(0),  
                               min.features = c(1),
                               modal.order = c('rna'), 
                               sparce = TRUE)
rm(data)
for(i in 1:length(musa_obj@Assays)){
    musa_obj@Assays[[i]]@data <- musa_obj@Assays[[i]]@Raw
}
musa_obj <- Add.HVFs(musa_obj, modal = "rna")
meta <- meta[musa_obj@Meta$cell.name,]
musa_obj@Meta[['celltype']] <- meta$CellType
all.equal(musa_obj@Meta$cell.name,rownames(meta))

musa_obj <- DownSample(musa_obj,
                       modal = c("rna"),
                       supervised = T,
                       group = 'celltype',
                       dims = list(1:50),
                       method = c("PCA"))

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c('rna'),
                           lambda = list(0.8),
                           supervised = T,
                           L2.norm = F,
                           sub.dims = list(1:20),
                           cos.dims = list(1:50),
                           Angle.var = c(15),
                           max.Angle = c(50))

emb <- list()
for(i in 1:9){
    emb[[i]] <- musa_obj@Assays[[i]]@embedding
}
emb <- Reduce(cbind,emb)
emb <- emb[,rownames(meta)]

out_dir <- './'
write.csv(as.matrix(t(emb)),
          paste0(out_dir,'/Palette_supervised.csv'),
          row.names = FALSE)
