source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(Matrix)

setwd('./')

rna <- as(as.matrix(readRDS('./RNA/rna.rds')),'dgCMatrix')
cytof <- as(as.matrix(readRDS('./cytof/cytof.rds')),'dgCMatrix')
rna_co <- as.matrix(readRDS('./RNA/co_q.rds'))
cytof_co <- as.matrix(readRDS('./cytof/co_q.rds'))

musa_obj <- Create.Musa.Object(data.list = list(rna,rna_co,cytof,cytof_co), 
                               samples = c('B1','B1','B2','B2'), 
                               modals = c("rna","co","adt","co"),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna",'co','adt'), 
                               sparce = TRUE)

musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                          normal.method = "LogNormalize")
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "adt",
                          normal.method = "CLR",
                          margin = 2)
musa_obj@Assays[[2]]@data <- rna_co
musa_obj@Assays[[4]]@data <- cytof_co

musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj, modal = "adt")
musa_obj <- Add.HVFs(musa_obj, modal = "co")

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("co"),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:15),# 使用哪些维数
                         method = c("PCA"),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

musa_obj <- Find.Subspace3(musa_obj,
                          modal = c("co"),
                          lambda = list(0.9),
                          alpha = list(0),
                         sub.dims = list(1:15),
                         cos.dims = list(1:15),
                          Angle.var = c(30),
                         max.Angle = c(60))

musa_obj <- IsolatedModalEmbedding(musa_obj,modal = c("adt",'rna'),
                                   dims = list(1:20,1:30))
a <- Run.Palette2(musa_obj,nn = 2,lambda = 0.5,nDims = 20,modal.norm = F)

out_dir <- './'

write.csv(as.matrix(a@Int.result[["bind"]]),
          paste0(out_dir,'Palette.csv'),
          row.names = FALSE)