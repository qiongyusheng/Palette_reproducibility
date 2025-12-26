library(Seurat)
library(cowplot)
library(dplyr)
library(data.table)
library(tibble)
library(Matrix)

setwd('./')
data <- readRDS('./dataset.rds')
varg <- readRDS('./vargs.rds')
data <- lapply(data,function(x) x[varg,])
meta <- readRDS('./meta.rds')

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
                               samples = c('B1','B2','B3','B4','B5'), 
                               modals = c('rna','rna','rna','rna','rna'),
                               filter=T,
                               min.cells = c(0),  
                               min.features = c(1),
                               modal.order = c('rna'), 
                               sparce = TRUE)
rm(data)
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                          normal.method = "LogNormalize")
musa_obj <- Add.HVFs(musa_obj, modal = "rna")

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c('rna'),
                         sample = NULL,
                         dim.reduce = TRUE,
                         dims = list(1:50),
                         method = c("PCA"),
                         joint = F,
                         resolution = 1,
                         nn.k = 20,
                         nn.method = "annoy",
                         annoy.metric = "euclidean",
                         verbose = TRUE)
musa_obj <- Find.Subspace3(musa_obj,
                           modal = c('rna'),
                           lambda = list(0.8),
                           supervised = F,
                           L2.norm = T,
                           sub.dims = list(1:20),
                           cos.dims = list(1:50),
                           Angle.var = c(15),
                           max.Angle = c(50))

emb <- list()
for(i in 1:5){
    emb[[i]] <- musa_obj@Assays[[i]]@embedding
}
emb <- Reduce(cbind,emb)
emb <- emb[,rownames(meta)]

out_dir <- './'
write.csv(as.matrix(t(emb)),
          paste0(out_dir,'/Palette_unsupervised.csv'),
          row.names = FALSE)
