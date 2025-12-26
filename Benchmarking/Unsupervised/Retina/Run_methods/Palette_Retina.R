######################### Random 1 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
options(future.globals.maxSize = 50 * 1024^3)
library(Matrix)
setwd('./Human_Retina')
path <- c('./LGS1OD','./LGS1OS','./LGS2OD','./LGS2OS','./LGS3OD','./LGS3OS','./LVG1OD','./LVG1OS')
library(data.table)
rna.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:8){
  if(i %in% c(3,4,5,6,8)){
    rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds')) 
  }else{
    rna.list[[i]] <- NULL 
  }
  
  if(i %in% c(1,2,4,6,7)){
    peak <- fread(paste0(path[i],'/peak.csv'),data.table = F)
    cells <- fread(paste0(path[i],'/bcd.csv'),data.table = F)
    atac.list[[i]] <- readMM(paste0(path[i],'/atac.mtx'))
    colnames(atac.list[[i]]) <- cells[,1]
    rownames(atac.list[[i]]) <- peak[,1]
    # atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
  }else{
    atac.list[[i]] <- NULL 
  }
  # rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds'))
  # atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
}
options(future.globals.maxSize = 500 * 1024^3)
musa_obj <- Create.Musa.Object(data.list = c(rna.list[c(3,4,5,6,8)],atac.list[c(1,2,4,6,7)]), 
                               samples = c('b3','b4','b5','b6','b8',
                                           'b1','b2','b4','b6','b7'), 
                               modals = c("rna","rna",'rna','rna','rna',
                                          'atac','atac','atac','atac','atac'),
                               filter=T,
                               min.cells = c(0,0),  
                               min.features = c(1,1),
                               modal.order = c("rna",'atac'), 
                               sparce = TRUE)
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- IDFLog_try(musa_obj,modal = "atac")
musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj,modal = "atac")

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'atac'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:50),# 使用哪些维数
                         method = c("PCA",'LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)
options(future.globals.maxSize = 50 * 1024^3)
musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'atac'),
                           lambda = list(0.8,0.5),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30),
                           cos.dims = list(1:50,1:50),
                           Angle.var = c(15,20),
                           max.Angle = c(50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
mm <- data.frame(batch = c(rep("B1",8455),rep("B2",9182),rep("B3",6234),rep("B4",5955),
                           rep("B5",5707),rep("B6",5013),rep("B7",5037),rep("B8",4729)),
                 modality = c(rep('atac',8455+9182),
                              rep('rna',6234),
                              rep('multiome',5955),
                              rep('rna',5707),
                              rep('multiome',5013),
                              rep('atac',5037),
                              rep('rna',4729)),
                 celltype = c(as.character(meta.list[[1]][,4]),
                              as.character(meta.list[[2]][,4]),
                              as.character(meta.list[[3]][,4]),
                              as.character(meta.list[[4]][,4]),
                              as.character(meta.list[[5]][,4]),
                              as.character(meta.list[[6]][,4]),
                              as.character(meta.list[[7]][,4]),
                              as.character(meta.list[[8]][,4]))
)

meta <- Reduce(rbind,meta.list)
rownames(mm) <- rownames(meta)
out_dir <- './Benchmarking/Unsupervised/Retina/output'
d <- a@Int.result[["bind"]]
d <- d[rownames(mm),]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_s1.csv'),
          row.names = FALSE)

saveRDS(mm,paste0(out_dir,'/meta_s1.rds'))

######################### Random 2 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
options(future.globals.maxSize = 50 * 1024^3)
library(Matrix)
setwd('./Human_Retina')
path <- c('./LGS1OD','./LGS1OS','./LGS2OD','./LGS2OS','./LGS3OD','./LGS3OS','./LVG1OD','./LVG1OS')
library(data.table)

rna.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:8){
  if(i %in% c(2,3,4,5,7)){
    rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds')) 
  }else{
    rna.list[[i]] <- NULL 
  }
  
  if(i %in% c(1,2,3,6,8)){
    peak <- fread(paste0(path[i],'/peak.csv'),data.table = F)
    cells <- fread(paste0(path[i],'/bcd.csv'),data.table = F)
    atac.list[[i]] <- readMM(paste0(path[i],'/atac.mtx'))
    colnames(atac.list[[i]]) <- cells[,1]
    rownames(atac.list[[i]]) <- peak[,1]
    # atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
  }else{
    atac.list[[i]] <- NULL 
  }
  # rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds'))
  # atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
}
options(future.globals.maxSize = 500 * 1024^3)
musa_obj <- Create.Musa.Object(data.list = c(rna.list[c(2,3,4,5,7)],atac.list[c(1,2,3,6,8)]), 
                               samples = c('b2','b3','b4','b5','b7',
                                           'b1','b2','b3','b6','b8'), 
                               modals = c("rna","rna",'rna','rna','rna',
                                          'atac','atac','atac','atac','atac'),
                               filter=T,
                               min.cells = c(0,0),  
                               min.features = c(1,1),
                               modal.order = c("rna",'atac'), 
                               sparce = TRUE)

musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- IDFLog_try(musa_obj,modal = "atac")
musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj,modal = "atac")
musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'atac'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:50),# 使用哪些维数
                         method = c("PCA",'LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)
options(future.globals.maxSize = 50 * 1024^3)
musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'atac'),
                           lambda = list(0.8,0.5),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30),
                           cos.dims = list(1:50,1:50),
                           Angle.var = c(15,20),
                           max.Angle = c(50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
mm <- data.frame(batch = c(rep("B1",8455),rep("B2",9182),rep("B3",6234),rep("B4",5955),
                           rep("B5",5707),rep("B6",5013),rep("B7",5037),rep("B8",4729)),
                 modality = c(rep('rna',8455),
                              rep('multiome',9182+6234),
                              rep('rna',5955+5707),
                              rep('atac',5013),
                              rep('rna',5037),
                              rep('atac',4729)),
                 celltype = c(as.character(meta.list[[1]][,4]),
                              as.character(meta.list[[2]][,4]),
                              as.character(meta.list[[3]][,4]),
                              as.character(meta.list[[4]][,4]),
                              as.character(meta.list[[5]][,4]),
                              as.character(meta.list[[6]][,4]),
                              as.character(meta.list[[7]][,4]),
                              as.character(meta.list[[8]][,4]))
)

meta <- Reduce(rbind,meta.list)
rownames(mm) <- rownames(meta)
out_dir <- './Benchmarking/Unsupervised/Retina/output'
d <- a@Int.result[["bind"]]
d <- d[rownames(mm),]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_s2.csv'),
          row.names = FALSE)

saveRDS(mm,paste0(out_dir,'/meta_s2.rds'))


######################### Random 3 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
options(future.globals.maxSize = 50 * 1024^3)
library(Matrix)
setwd('./Human_Retina')
path <- c('./LGS1OD','./LGS1OS','./LGS2OD','./LGS2OS','./LGS3OD','./LGS3OS','./LVG1OD','./LVG1OS')
library(data.table)
rna.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:8){
  if(i %in% c(1,2,3,4,5)){
    rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds')) 
  }else{
    rna.list[[i]] <- NULL 
  }
  
  if(i %in% c(1,2,6,7,8)){
    peak <- fread(paste0(path[i],'/peak.csv'),data.table = F)
    cells <- fread(paste0(path[i],'/bcd.csv'),data.table = F)
    atac.list[[i]] <- readMM(paste0(path[i],'/atac.mtx'))
    colnames(atac.list[[i]]) <- cells[,1]
    rownames(atac.list[[i]]) <- peak[,1]
    # atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
  }else{
    atac.list[[i]] <- NULL 
  }
  # rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds'))
  # atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
}
options(future.globals.maxSize = 500 * 1024^3)
musa_obj <- Create.Musa.Object(data.list = c(rna.list[c(1,2,3,4,5)],atac.list[c(1,2,6,7,8)]), 
                               samples = c('b1','b2','b3','b4','b5',
                                           'b1','b2','b6','b7','b8'), 
                               modals = c("rna","rna",'rna','rna','rna',
                                          'atac','atac','atac','atac','atac'),
                               filter=T,
                               min.cells = c(0,0),  
                               min.features = c(1,1),
                               modal.order = c("rna",'atac'), 
                               sparce = TRUE)
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- IDFLog_try(musa_obj,modal = "atac")
musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj,modal = "atac")
musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'atac'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:50),# 使用哪些维数
                         method = c("PCA",'LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)
options(future.globals.maxSize = 50 * 1024^3)
musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'atac'),
                           lambda = list(0.8,0.5),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30),
                           cos.dims = list(1:50,1:50),
                           Angle.var = c(15,20),
                           max.Angle = c(50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
mm <- data.frame(batch = c(rep("B1",8455),rep("B2",9182),rep("B3",6234),rep("B4",5955),
                           rep("B5",5707),rep("B6",5013),rep("B7",5037),rep("B8",4729)),
                 modality = c(rep('multiome',8455+9182),rep('rna',6234+5955+5707),rep('atac',5013+5037+4729)),
                 celltype = c(as.character(meta.list[[1]][,4]),
                              as.character(meta.list[[2]][,4]),
                              as.character(meta.list[[3]][,4]),
                              as.character(meta.list[[4]][,4]),
                              as.character(meta.list[[5]][,4]),
                              as.character(meta.list[[6]][,4]),
                              as.character(meta.list[[7]][,4]),
                              as.character(meta.list[[8]][,4]))
)
meta <- Reduce(rbind,meta.list)
rownames(mm) <- rownames(meta)
out_dir <- './Benchmarking/Unsupervised/Retina/output'
d <- a@Int.result[["bind"]]
d <- d[rownames(mm),]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette.csv'),
          row.names = FALSE)

saveRDS(mm,paste0(out_dir,'/meta.rds'))

######################### Random 4 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
options(future.globals.maxSize = 50 * 1024^3)
library(Matrix)
setwd('./Human_Retina')
path <- c('./LGS1OD','./LGS1OS','./LGS2OD','./LGS2OS','./LGS3OD','./LGS3OS','./LVG1OD','./LVG1OS')
library(data.table)
rna.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:8){
  if(i %in% c(2,4,5,6,7)){
    rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds')) 
  }else{
    rna.list[[i]] <- NULL 
  }
  
  if(i %in% c(1,3,4,5,8)){
    peak <- fread(paste0(path[i],'/peak.csv'),data.table = F)
    cells <- fread(paste0(path[i],'/bcd.csv'),data.table = F)
    atac.list[[i]] <- readMM(paste0(path[i],'/atac.mtx'))
    colnames(atac.list[[i]]) <- cells[,1]
    rownames(atac.list[[i]]) <- peak[,1]
    # atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
  }else{
    atac.list[[i]] <- NULL 
  }
  # rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds'))
  # atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
}
options(future.globals.maxSize = 500 * 1024^3)
musa_obj <- Create.Musa.Object(data.list = c(rna.list[c(2,4,5,6,7)],atac.list[c(1,3,4,5,8)]), 
                               samples = c('b2','b4','b5','b6','b7',
                                           'b1','b3','b4','b5','b8'), 
                               modals = c("rna","rna",'rna','rna','rna',
                                          'atac','atac','atac','atac','atac'),
                               filter=T,
                               min.cells = c(0,0),  
                               min.features = c(1,1),
                               modal.order = c("rna",'atac'), 
                               sparce = TRUE)

musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- IDFLog_try(musa_obj,modal = "atac")
musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj,modal = "atac")
musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'atac'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:50),# 使用哪些维数
                         method = c("PCA",'LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)
options(future.globals.maxSize = 50 * 1024^3)
musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'atac'),
                           lambda = list(0.8,0.5),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30),
                           cos.dims = list(1:50,1:50),
                           Angle.var = c(15,20),
                           max.Angle = c(50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)

mm <- data.frame(batch = c(rep("B1",8455),rep("B2",9182),rep("B3",6234),rep("B4",5955),
                           rep("B5",5707),rep("B6",5013),rep("B7",5037),rep("B8",4729)),
                 modality = c(rep('atac',8455),
                              rep('rna',9182),
                              rep('atac',6234),
                              rep('multiome',5955+5707),
                              rep('rna',5013+5037),
                              rep('atac',4729)),
                 celltype = c(as.character(meta.list[[1]][,4]),
                              as.character(meta.list[[2]][,4]),
                              as.character(meta.list[[3]][,4]),
                              as.character(meta.list[[4]][,4]),
                              as.character(meta.list[[5]][,4]),
                              as.character(meta.list[[6]][,4]),
                              as.character(meta.list[[7]][,4]),
                              as.character(meta.list[[8]][,4]))
)


meta <- Reduce(rbind,meta.list)
rownames(mm) <- rownames(meta)
out_dir <- './Benchmarking/Unsupervised/Retina/output'
d <- a@Int.result[["bind"]]
d <- d[rownames(mm),]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_s3.csv'),
          row.names = FALSE)

saveRDS(mm,paste0(out_dir,'/meta_s3.rds'))


######################### Random 5 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
options(future.globals.maxSize = 50 * 1024^3)
library(Matrix)
setwd('./Human_Retina')
path <- c('./LGS1OD','./LGS1OS','./LGS2OD','./LGS2OS','./LGS3OD','./LGS3OS','./LVG1OD','./LVG1OS')
library(data.table)

rna.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:8){
  if(i %in% c(1,3,5,6,7)){
    rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds')) 
  }else{
    rna.list[[i]] <- NULL 
  }
  
  if(i %in% c(2,4,6,7,8)){
    peak <- fread(paste0(path[i],'/peak.csv'),data.table = F)
    cells <- fread(paste0(path[i],'/bcd.csv'),data.table = F)
    atac.list[[i]] <- readMM(paste0(path[i],'/atac.mtx'))
    colnames(atac.list[[i]]) <- cells[,1]
    rownames(atac.list[[i]]) <- peak[,1]
    # atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
  }else{
    atac.list[[i]] <- NULL 
  }
  # rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds'))
  # atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
}
options(future.globals.maxSize = 500 * 1024^3)
musa_obj <- Create.Musa.Object(data.list = c(rna.list[c(1,3,5,6,7)],atac.list[c(2,4,6,7,8)]), 
                               samples = c('b1','b3','b5','b6','b7',
                                           'b2','b4','b6','b7','b8'), 
                               modals = c("rna","rna",'rna','rna','rna',
                                          'atac','atac','atac','atac','atac'),
                               filter=T,
                               min.cells = c(0,0),  
                               min.features = c(1,1),
                               modal.order = c("rna",'atac'), 
                               sparce = TRUE)

musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- IDFLog_try(musa_obj,modal = "atac")
musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj,modal = "atac")
musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'atac'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:50),# 使用哪些维数
                         method = c("PCA",'LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)
options(future.globals.maxSize = 50 * 1024^3)
musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'atac'),
                           lambda = list(0.8,0.5),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30),
                           cos.dims = list(1:50,1:50),
                           Angle.var = c(15,20),
                           max.Angle = c(50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
mm <- data.frame(batch = c(rep("B1",8455),rep("B2",9182),rep("B3",6234),rep("B4",5955),
                           rep("B5",5707),rep("B6",5013),rep("B7",5037),rep("B8",4729)),
                 modality = c(rep('rna',8455),
                              rep('atac',9182),
                              rep('rna',6234),
                              rep('atac',5955),
                              rep('rna',5707),
                              rep('multiome',5013+5037),
                              rep('atac',4729)),
                 celltype = c(as.character(meta.list[[1]][,4]),
                              as.character(meta.list[[2]][,4]),
                              as.character(meta.list[[3]][,4]),
                              as.character(meta.list[[4]][,4]),
                              as.character(meta.list[[5]][,4]),
                              as.character(meta.list[[6]][,4]),
                              as.character(meta.list[[7]][,4]),
                              as.character(meta.list[[8]][,4]))
)
meta <- Reduce(rbind,meta.list)
rownames(mm) <- rownames(meta)
out_dir <- '/home/server/sqy/data/mosaic/result/bench/new_bench/Retina'
d <- a@Int.result[["bind"]]
d <- d[rownames(mm),]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_s4.csv'),
          row.names = FALSE)

saveRDS(mm,paste0(out_dir,'/meta_s4.rds'))

