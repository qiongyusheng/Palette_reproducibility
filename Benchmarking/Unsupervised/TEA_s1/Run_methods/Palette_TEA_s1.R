######################### Random 1 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
setwd("./PBMC_TEA")
rna.list <- list()
atac.list <- list()
adt.list <- list()
meta.list <- list()
for(i in 2:5){
  rna.list[[i-1]] <- readRDS(paste0("./B",i,"/rna.rds"))
  meta.list[[i-1]] <- readRDS(paste0("./B",i,"/meta.rds"))
  atac.list[[i-1]] <- readRDS(paste0("./B",i,"/atac.rds"))
  adt.list[[i-1]] <- readRDS(paste0("./B",i,"/adt.rds"))
}
musa_obj <- Create.Musa.Object(data.list = c(atac.list[c(1,4)],rna.list[c(3,4)],adt.list[c(2,4)]), 
                               samples = c('B1','B4','B3','B4','B2','B4'), 
                               modals = c("atac","atac","rna","rna",'adt','adt'),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna",'atac','adt'), 
                               sparce = TRUE)

musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "adt",
                           normal.method = "CLR",
                           margin = 2)
musa_obj <- IDFLog_try(musa_obj,modal = "atac")
musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj, modal = "adt")
musa_obj <- Add.HVFs(musa_obj, modal = "atac")
musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'adt','atac'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:20,1:50),# 使用哪些维数
                         method = c("PCA",'PCA','LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)
options(future.globals.maxSize = 50 * 1024^3)
musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:20,1:30),
                           cos.dims = list(1:50,1:20,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
mm <- data.frame(batch = c(rep("B1",6194),rep("B2",6381),rep("B3",6412),rep("B4",6530)),
                 celltype = c(as.character(meta.list[[1]][,2]),
                              as.character(meta.list[[2]][,2]),as.character(meta.list[[3]][,2]),
                              as.character(meta.list[[4]][,2]))
)
meta <- Reduce(rbind,meta.list)
rownames(mm) <- rownames(meta)
out_dir <- './Benchmarking/Unsupervised/TEA_s1/output'
d <- a@Int.result[["bind"]]
d <- d[rownames(mm),]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_v1.csv'),
          row.names = FALSE)

saveRDS(mm,paste0(out_dir,'/meta_v1.rds'))

######################### Random 2 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
setwd("./PBMC_TEA")
rna.list <- list()
atac.list <- list()
adt.list <- list()
meta.list <- list()
for(i in 2:5){
  rna.list[[i-1]] <- readRDS(paste0("./B",i,"/rna.rds"))
  meta.list[[i-1]] <- readRDS(paste0("./B",i,"/meta.rds"))
  atac.list[[i-1]] <- readRDS(paste0("./B",i,"/atac.rds"))
  adt.list[[i-1]] <- readRDS(paste0("./B",i,"/adt.rds"))
}
musa_obj <- Create.Musa.Object(data.list = c(atac.list[1:2],rna.list[c(1,3)],adt.list[c(1,4)]), 
                               samples = c('B1','B2','B1','B3','B1','B4'), 
                               modals = c("atac","atac","rna","rna",'adt','adt'),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna",'atac','adt'), 
                               sparce = TRUE)

musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "adt",
                           normal.method = "CLR",
                           margin = 2)
musa_obj <- IDFLog_try(musa_obj,modal = "atac")
musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj, modal = "adt")
musa_obj <- Add.HVFs(musa_obj, modal = "atac")
musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'adt','atac'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:20,1:50),# 使用哪些维数
                         method = c("PCA",'PCA','LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)
options(future.globals.maxSize = 50 * 1024^3)
musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:20,1:30),
                           cos.dims = list(1:50,1:20,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
mm <- data.frame(batch = c(rep("B1",6194),rep("B2",6381),rep("B3",6412),rep("B4",6530)),
                 celltype = c(as.character(meta.list[[1]][,2]),
                              as.character(meta.list[[2]][,2]),as.character(meta.list[[3]][,2]),
                              as.character(meta.list[[4]][,2]))
)
meta <- Reduce(rbind,meta.list)
rownames(mm) <- rownames(meta)
out_dir <- './Benchmarking/Unsupervised/TEA_s1/output'
d <- a@Int.result[["bind"]]
d <- d[rownames(mm),]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette.csv'),
          row.names = FALSE)

saveRDS(mm,paste0(out_dir,'/meta.rds'))

######################### Random 3 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
setwd("./PBMC_TEA")
rna.list <- list()
atac.list <- list()
adt.list <- list()
meta.list <- list()
for(i in 2:5){
  rna.list[[i-1]] <- readRDS(paste0("./B",i,"/rna.rds"))
  meta.list[[i-1]] <- readRDS(paste0("./B",i,"/meta.rds"))
  atac.list[[i-1]] <- readRDS(paste0("./B",i,"/atac.rds"))
  adt.list[[i-1]] <- readRDS(paste0("./B",i,"/adt.rds"))
}

musa_obj <- Create.Musa.Object(data.list = c(atac.list[c(3,4)],rna.list[c(1,3)],adt.list[c(2,3)]), 
                               samples = c('B3','B4','B1','B3','B2','B3'), 
                               modals = c("atac","atac","rna","rna",'adt','adt'),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna",'atac','adt'), 
                               sparce = TRUE)
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "adt",
                           normal.method = "CLR",
                           margin = 2)
musa_obj <- IDFLog_try(musa_obj,modal = "atac")
musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj, modal = "adt")
musa_obj <- Add.HVFs(musa_obj, modal = "atac")
musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'adt','atac'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:20,1:50),# 使用哪些维数
                         method = c("PCA",'PCA','LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)
options(future.globals.maxSize = 50 * 1024^3)
musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:20,1:30),
                           cos.dims = list(1:50,1:20,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
mm <- data.frame(batch = c(rep("B1",6194),rep("B2",6381),rep("B3",6412),rep("B4",6530)),
                 celltype = c(as.character(meta.list[[1]][,2]),
                              as.character(meta.list[[2]][,2]),as.character(meta.list[[3]][,2]),
                              as.character(meta.list[[4]][,2]))
)
meta <- Reduce(rbind,meta.list)
rownames(mm) <- rownames(meta)
out_dir <- './Benchmarking/Unsupervised/TEA_s1/output'
d <- a@Int.result[["bind"]]
d <- d[rownames(mm),]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_v2.csv'),
          row.names = FALSE)

saveRDS(mm,paste0(out_dir,'/meta_v2.rds'))

######################### Random 4 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
setwd("./PBMC_TEA")
rna.list <- list()
atac.list <- list()
adt.list <- list()
meta.list <- list()
for(i in 2:5){
  rna.list[[i-1]] <- readRDS(paste0("./B",i,"/rna.rds"))
  meta.list[[i-1]] <- readRDS(paste0("./B",i,"/meta.rds"))
  atac.list[[i-1]] <- readRDS(paste0("./B",i,"/atac.rds"))
  adt.list[[i-1]] <- readRDS(paste0("./B",i,"/adt.rds"))
}
musa_obj <- Create.Musa.Object(data.list = c(atac.list[c(1,3)],rna.list[c(1,4)],adt.list[c(1,2)]), 
                               samples = c('B1','B3','B1','B4','B1','B2'), 
                               modals = c("atac","atac","rna","rna",'adt','adt'),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna",'atac','adt'), 
                               sparce = TRUE)

musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "adt",
                           normal.method = "CLR",
                           margin = 2)
musa_obj <- IDFLog_try(musa_obj,modal = "atac")
musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj, modal = "adt")
musa_obj <- Add.HVFs(musa_obj, modal = "atac")
musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'adt','atac'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:20,1:50),# 使用哪些维数
                         method = c("PCA",'PCA','LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)
options(future.globals.maxSize = 50 * 1024^3)
musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:20,1:30),
                           cos.dims = list(1:50,1:20,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
mm <- data.frame(batch = c(rep("B1",6194),rep("B2",6381),rep("B3",6412),rep("B4",6530)),
                 celltype = c(as.character(meta.list[[1]][,2]),
                              as.character(meta.list[[2]][,2]),as.character(meta.list[[3]][,2]),
                              as.character(meta.list[[4]][,2]))
)
meta <- Reduce(rbind,meta.list)
rownames(mm) <- rownames(meta)
out_dir <- './Benchmarking/Unsupervised/TEA_s1/output'
d <- a@Int.result[["bind"]]
d <- d[rownames(mm),]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_v3.csv'),
          row.names = FALSE)

saveRDS(mm,paste0(out_dir,'/meta_v3.rds'))

######################### Random 5 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
setwd("./PBMC_TEA")
rna.list <- list()
atac.list <- list()
adt.list <- list()
meta.list <- list()
for(i in 2:5){
  rna.list[[i-1]] <- readRDS(paste0("./B",i,"/rna.rds"))
  meta.list[[i-1]] <- readRDS(paste0("./B",i,"/meta.rds"))
  atac.list[[i-1]] <- readRDS(paste0("./B",i,"/atac.rds"))
  adt.list[[i-1]] <- readRDS(paste0("./B",i,"/adt.rds"))
}

musa_obj <- Create.Musa.Object(data.list = c(atac.list[c(1,2)],rna.list[c(1,4)],adt.list[c(1,3)]), 
                               samples = c('B1','B2','B1','B4','B1','B3'), 
                               modals = c("atac","atac","rna","rna",'adt','adt'),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna",'atac','adt'), 
                               sparce = TRUE)
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "adt",
                           normal.method = "CLR",
                           margin = 2)
musa_obj <- IDFLog_try(musa_obj,modal = "atac")
musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj, modal = "adt")
musa_obj <- Add.HVFs(musa_obj, modal = "atac")

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'adt','atac'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:20,1:50),# 使用哪些维数
                         method = c("PCA",'PCA','LSI'),
                         joint = F,
                         resolution = 1,
                         nn.k = 20,
                         nn.method = "annoy",
                         annoy.metric = "euclidean",
                         verbose = TRUE)
options(future.globals.maxSize = 50 * 1024^3)
musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:20,1:30),
                           cos.dims = list(1:50,1:20,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
mm <- data.frame(batch = c(rep("B1",6194),rep("B2",6381),rep("B3",6412),rep("B4",6530)),
                 celltype = c(as.character(meta.list[[1]][,2]),
                              as.character(meta.list[[2]][,2]),as.character(meta.list[[3]][,2]),
                              as.character(meta.list[[4]][,2]))
)
meta <- Reduce(rbind,meta.list)
rownames(mm) <- rownames(meta)
out_dir <- './Benchmarking/Unsupervised/TEA_s1/output'
d <- a@Int.result[["bind"]]
d <- d[rownames(mm),]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_v4.csv'),
          row.names = FALSE)

saveRDS(mm,paste0(out_dir,'/meta_v4.rds'))

