library(data.table)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Matrix)
rna1 <- readRDS('./10x_Visium_human_tonsil/slice1/rna.rds')
rna2 <- readRDS('./10x_Visium_human_tonsil/slice2/rna.rds')
rna1 <- NormalizeData(rna1)
rna2 <- NormalizeData(rna2)
adt1 <- readRDS('./10x_Visium_human_tonsil/slice1/adt.rds')
adt3 <- readRDS('./10x_Visium_human_tonsil/slice3/adt.rds')
adt1 <- NormalizeData(adt1,normalization.method = 'CLR', margin = 2)
adt3 <- NormalizeData(adt3,normalization.method = 'CLR', margin = 2)
meta1 <- readRDS('./10x_Visium_human_tonsil/slice1/meta.rds')
meta2 <- readRDS('./10x_Visium_human_tonsil/slice2/meta.rds')
meta3 <- readRDS('./10x_Visium_human_tonsil/slice3/meta.rds')

source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
options(future.globals.maxSize = 50 * 1024^3)

musa_obj <- Create.Musa.Object(data.list = list(rna1,rna2,adt1,adt3), 
                               samples = c('b1','b2','b1','b3'), 
                               modals = c("rna","rna",
                                          "adt",'adt'),
                               filter=F,
                               min.cells = c(0,0),  
                               min.features = c(1,1),
                               modal.order = c("rna",'adt'), 
                               sparce = TRUE)
musa_obj@Assays[[1]]@data <- musa_obj@Assays[[1]]@Raw
musa_obj@Assays[[2]]@data <- musa_obj@Assays[[2]]@Raw
musa_obj@Assays[[3]]@data <- musa_obj@Assays[[3]]@Raw
musa_obj@Assays[[4]]@data <- musa_obj@Assays[[4]]@Raw

musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj, modal = "adt")

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'adt'),
                         dims = list(1:50,1:20),
                         method = c("PCA",'PCA'))

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt'),
                           lambda = list(0.8,0.8),
                           sub.dims = list(1:20,1:15),
                           cos.dims = list(1:50,1:20))

musa_obj <- Run.Palette2(musa_obj)

meta1$batch <- 'Batch1'
meta2$batch <- 'Batch2'
meta3$batch <- 'Batch3'
meta <- rbind(meta1,meta2,meta3)

emb <- musa_obj@Int.result[["bind"]]
emb <- emb[rownames(meta),]
write.csv(as.matrix(emb),
          './10xVisium_tonsil/Palette.csv',
          row.names = FALSE)

