######################### Random 1 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

setwd('./Young_Age')
out_dir <- './Benchmarking/Unsupervised/Ab-seq/output'

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in paths){
  adt.list[[i]] <- readRDS(paste0('./',i,'/adt.rds'))
  rna.list[[i]] <- readRDS(paste0('./',i,'/rna.rds'))
  meta.list[[i]] <- readRDS(paste0('./',i,'/meta.rds'))
}
names(adt.list) <- paste0(names(adt.list),'-adt')
names(rna.list) <- paste0(names(rna.list),'-rna')
musa_obj <- Create.Musa.Object(data.list = c(adt.list[c(1,3,4,5)],rna.list[c(1,2,5,6)]), 
                               samples = c('Age1','Age3','BM1','BM2',
                                           'Age1','Age2','BM2','BM3'), 
                               modals = c("adt","adt","adt","adt",
                                          "rna","rna","rna","rna"),
                               filter=T,
                               min.cells = c(0,0),  
                               min.features = c(1,1),
                               modal.order = c("rna",'adt'), 
                               sparce = TRUE)
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "adt",
                           normal.method = "CLR",
                           margin = 2)

musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj,modal = "adt")

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'adt'),
                         sample = NULL,
                         dim.reduce = TRUE,
                         dims = list(1:50,1:50),
                         method = c("PCA","PCA"),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt'),
                           lambda = list(0.8,0.8),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30),
                           cos.dims = list(1:50,1:50),
                           Angle.var = c(15,15),
                           max.Angle = c(50,50))

a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5, emb_L2 = T)

cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))

d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_s1.csv'),
          row.names = FALSE)
mm <- data.frame(batch = c(rep("Age1",6316),rep("Age2",4071),rep("Age3",10639),
                           rep("BM1",9751),rep("BM2",6963),rep("BM3",11317)),
                 ID = c(as.character(meta.list[[1]][,1]),
                        as.character(meta.list[[2]][,1]),
                        as.character(meta.list[[3]][,1]),
                        as.character(meta.list[[4]][,1]),
                        as.character(meta.list[[5]][,1]),
                        as.character(meta.list[[6]][,1])),
                 celltype.l2 = c(as.character(meta.list[[1]][,4]),
                                 as.character(meta.list[[2]][,4]),
                                 as.character(meta.list[[3]][,4]),
                                 as.character(meta.list[[4]][,4]),
                                 as.character(meta.list[[5]][,4]),
                                 as.character(meta.list[[6]][,4])),
                 modality = c(as.character(meta.list[[1]][,6]),
                              as.character(meta.list[[2]][,6]),
                              as.character(meta.list[[3]][,6]),
                              as.character(meta.list[[4]][,6]),
                              as.character(meta.list[[5]][,6]),
                              as.character(meta.list[[6]][,6]))
)
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

setwd('./Young_Age')
out_dir <- './Benchmarking/Unsupervised/Ab-seq/output'

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in paths){
  adt.list[[i]] <- readRDS(paste0('./',i,'/adt.rds'))
  rna.list[[i]] <- readRDS(paste0('./',i,'/rna.rds'))
  meta.list[[i]] <- readRDS(paste0('./',i,'/meta.rds'))
}
names(adt.list) <- paste0(names(adt.list),'-adt')
names(rna.list) <- paste0(names(rna.list),'-rna')
musa_obj <- Create.Musa.Object(data.list = c(adt.list[c(1,3,5,6)],rna.list[c(2,3,4,6)]), 
                               samples = c('Age1','Age3','BM2','BM3',
                                           'Age2','Age3','BM1','BM3'), 
                               modals = c("adt","adt","adt","adt",
                                          "rna","rna","rna","rna"),
                               filter=T,
                               min.cells = c(0,0),  
                               min.features = c(1,1),
                               modal.order = c("rna",'adt'), 
                               sparce = TRUE)
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "adt",
                           normal.method = "CLR",
                           margin = 2)

musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj,modal = "adt")

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'adt'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:50),# 使用哪些维数
                         method = c("PCA","PCA"),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt'),
                           lambda = list(0.8,0.8),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30),
                           cos.dims = list(1:50,1:50),
                           Angle.var = c(15,15),
                           max.Angle = c(50,50))

a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5, emb_L2 = T)

cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))

d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_s2.csv'),
          row.names = FALSE)
mm <- data.frame(batch = c(rep("Age1",6316),rep("Age2",4071),rep("Age3",10639),
                           rep("BM1",9751),rep("BM2",6963),rep("BM3",11317)),
                 ID = c(as.character(meta.list[[1]][,1]),
                        as.character(meta.list[[2]][,1]),
                        as.character(meta.list[[3]][,1]),
                        as.character(meta.list[[4]][,1]),
                        as.character(meta.list[[5]][,1]),
                        as.character(meta.list[[6]][,1])),
                 celltype.l2 = c(as.character(meta.list[[1]][,4]),
                                 as.character(meta.list[[2]][,4]),
                                 as.character(meta.list[[3]][,4]),
                                 as.character(meta.list[[4]][,4]),
                                 as.character(meta.list[[5]][,4]),
                                 as.character(meta.list[[6]][,4])),
                 modality = c(as.character(meta.list[[1]][,6]),
                              as.character(meta.list[[2]][,6]),
                              as.character(meta.list[[3]][,6]),
                              as.character(meta.list[[4]][,6]),
                              as.character(meta.list[[5]][,6]),
                              as.character(meta.list[[6]][,6]))
)
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

setwd('./Young_Age')
out_dir <- './Benchmarking/Unsupervised/Ab-seq/output'

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in paths){
  adt.list[[i]] <- readRDS(paste0('./',i,'/adt.rds'))
  rna.list[[i]] <- readRDS(paste0('./',i,'/rna.rds'))
  meta.list[[i]] <- readRDS(paste0('./',i,'/meta.rds'))
}
names(adt.list) <- paste0(names(adt.list),'-adt')
names(rna.list) <- paste0(names(rna.list),'-rna')
musa_obj <- Create.Musa.Object(data.list = c(adt.list[c(2,3,4,5)],rna.list[c(1,2,4,6)]), 
                               samples = c('Age2','Age3','BM1','BM2',
                                           'Age1','Age2','BM1','BM3'), 
                               modals = c("adt","adt","adt","adt",
                                          "rna","rna","rna","rna"),
                               filter=T,
                               min.cells = c(0,0),  
                               min.features = c(1,1),
                               modal.order = c("rna",'adt'), 
                               sparce = TRUE)
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "adt",
                           normal.method = "CLR",
                           margin = 2)

musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj,modal = "adt")

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'adt'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:50),# 使用哪些维数
                         method = c("PCA","PCA"),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt'),
                           lambda = list(0.8,0.8),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30),
                           cos.dims = list(1:50,1:50),
                           Angle.var = c(15,15),
                           max.Angle = c(50,50))

a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5, emb_L2 = T)

cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))

d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_s3.csv'),
          row.names = FALSE)
mm <- data.frame(batch = c(rep("Age1",6316),rep("Age2",4071),rep("Age3",10639),
                           rep("BM1",9751),rep("BM2",6963),rep("BM3",11317)),
                 ID = c(as.character(meta.list[[1]][,1]),
                        as.character(meta.list[[2]][,1]),
                        as.character(meta.list[[3]][,1]),
                        as.character(meta.list[[4]][,1]),
                        as.character(meta.list[[5]][,1]),
                        as.character(meta.list[[6]][,1])),
                 celltype.l2 = c(as.character(meta.list[[1]][,4]),
                                 as.character(meta.list[[2]][,4]),
                                 as.character(meta.list[[3]][,4]),
                                 as.character(meta.list[[4]][,4]),
                                 as.character(meta.list[[5]][,4]),
                                 as.character(meta.list[[6]][,4])),
                 modality = c(as.character(meta.list[[1]][,6]),
                              as.character(meta.list[[2]][,6]),
                              as.character(meta.list[[3]][,6]),
                              as.character(meta.list[[4]][,6]),
                              as.character(meta.list[[5]][,6]),
                              as.character(meta.list[[6]][,6]))
)
saveRDS(mm,paste0(out_dir,'/meta_s3.rds'))

######################### Random 4 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

setwd('./Young_Age')
out_dir <- './Benchmarking/Unsupervised/Ab-seq/output'

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in paths){
  adt.list[[i]] <- readRDS(paste0('./',i,'/adt.rds'))
  rna.list[[i]] <- readRDS(paste0('./',i,'/rna.rds'))
  meta.list[[i]] <- readRDS(paste0('./',i,'/meta.rds'))
}
names(adt.list) <- paste0(names(adt.list),'-adt')
names(rna.list) <- paste0(names(rna.list),'-rna')

musa_obj <- Create.Musa.Object(data.list = c(adt.list[c(1,2,4,6)],rna.list[c(1,3,4,5)]), 
                               samples = c('Age1','Age2','BM1','BM3',
                                           'Age1','Age3','BM1','BM2'), 
                               modals = c("adt","adt","adt","adt",
                                          "rna","rna","rna","rna"),
                               filter=T,
                               min.cells = c(0,0),  
                               min.features = c(1,1),
                               modal.order = c("rna",'adt'), 
                               sparce = TRUE)

musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "adt",
                           normal.method = "CLR",
                           margin = 2)

musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj,modal = "adt")

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'adt'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:50),# 使用哪些维数
                         method = c("PCA","PCA"),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt'),
                           lambda = list(0.8,0.8),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30),
                           cos.dims = list(1:50,1:50),
                           Angle.var = c(15,15),
                           max.Angle = c(50,50))

a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5, emb_L2 = T)

cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))

d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette.csv'),
          row.names = FALSE)
mm <- data.frame(batch = c(rep("Age1",6316),rep("Age2",4071),rep("Age3",10639),
                           rep("BM1",9751),rep("BM2",6963),rep("BM3",11317)),
                 ID = c(as.character(meta.list[[1]][,1]),
                        as.character(meta.list[[2]][,1]),
                        as.character(meta.list[[3]][,1]),
                        as.character(meta.list[[4]][,1]),
                        as.character(meta.list[[5]][,1]),
                        as.character(meta.list[[6]][,1])),
                 celltype.l2 = c(as.character(meta.list[[1]][,4]),
                                 as.character(meta.list[[2]][,4]),
                                 as.character(meta.list[[3]][,4]),
                                 as.character(meta.list[[4]][,4]),
                                 as.character(meta.list[[5]][,4]),
                                 as.character(meta.list[[6]][,4])),
                 modality = c(as.character(meta.list[[1]][,6]),
                              as.character(meta.list[[2]][,6]),
                              as.character(meta.list[[3]][,6]),
                              as.character(meta.list[[4]][,6]),
                              as.character(meta.list[[5]][,6]),
                              as.character(meta.list[[6]][,6]))
)
saveRDS(mm,paste0(out_dir,'/meta.rds'))


######################### Random 5 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

setwd('./Young_Age')
out_dir <- './Benchmarking/Unsupervised/Ab-seq/output'

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in paths){
  adt.list[[i]] <- readRDS(paste0('./',i,'/adt.rds'))
  rna.list[[i]] <- readRDS(paste0('./',i,'/rna.rds'))
  meta.list[[i]] <- readRDS(paste0('./',i,'/meta.rds'))
}
names(adt.list) <- paste0(names(adt.list),'-adt')
names(rna.list) <- paste0(names(rna.list),'-rna')

musa_obj <- Create.Musa.Object(data.list = c(adt.list[c(1,2,4,5)],rna.list[c(1,3,5,6)]), 
                               samples = c('Age1','Age2','BM1','BM2',
                                           'Age1','Age3','BM2','BM3'), 
                               modals = c("adt","adt","adt","adt",
                                          "rna","rna","rna","rna"),
                               filter=T,
                               min.cells = c(0,0),  
                               min.features = c(1,1),
                               modal.order = c("rna",'adt'), 
                               sparce = TRUE)
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "adt",
                           normal.method = "CLR",
                           margin = 2)

musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj,modal = "adt")

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'adt'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:50),# 使用哪些维数
                         method = c("PCA","PCA"),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt'),
                           lambda = list(0.8,0.8),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30),
                           cos.dims = list(1:50,1:50),
                           Angle.var = c(15,15),
                           max.Angle = c(50,50))

a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5, emb_L2 = T)

cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))

d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_s4.csv'),
          row.names = FALSE)
mm <- data.frame(batch = c(rep("Age1",6316),rep("Age2",4071),rep("Age3",10639),
                           rep("BM1",9751),rep("BM2",6963),rep("BM3",11317)),
                 ID = c(as.character(meta.list[[1]][,1]),
                        as.character(meta.list[[2]][,1]),
                        as.character(meta.list[[3]][,1]),
                        as.character(meta.list[[4]][,1]),
                        as.character(meta.list[[5]][,1]),
                        as.character(meta.list[[6]][,1])),
                 celltype.l2 = c(as.character(meta.list[[1]][,4]),
                                 as.character(meta.list[[2]][,4]),
                                 as.character(meta.list[[3]][,4]),
                                 as.character(meta.list[[4]][,4]),
                                 as.character(meta.list[[5]][,4]),
                                 as.character(meta.list[[6]][,4])),
                 modality = c(as.character(meta.list[[1]][,6]),
                              as.character(meta.list[[2]][,6]),
                              as.character(meta.list[[3]][,6]),
                              as.character(meta.list[[4]][,6]),
                              as.character(meta.list[[5]][,6]),
                              as.character(meta.list[[6]][,6]))
)
saveRDS(mm,paste0(out_dir,'/meta_s4.rds'))