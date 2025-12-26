######################### Random 1 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
options(future.globals.maxSize = 5000 * 1024^3)
setwd('./Young_Age')
out_dir <- './Benchmarking/CellType_missing/output/Ab-seq'
meta <- readRDS('./Young_Age/meta.rds')

set.seed(123)
rm_ct <- sample(unique(meta$celltype.l2),size = 6,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 6,replace = T)
rm_ct

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in 1:length(paths)){
  tmp <- readRDS(paste0('./',paths[i],'/meta.rds'))
  meta.list[[paths[i]]] <- dplyr::filter(tmp,celltype.l2 != rm_ct[i])
  adt.list[[paths[i]]] <- readRDS(paste0('./',paths[i],'/adt.rds'))[,rownames(meta.list[[i]])]
  rna.list[[paths[i]]] <- readRDS(paste0('./',paths[i],'/rna.rds'))[,rownames(meta.list[[i]])]
  
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

mm <- data.frame(batch = c(rep("Age1",6192),rep("Age2",4035),rep("Age3",9380),
                           rep("BM1",8707),rep("BM2",5829),rep("BM3",11170)),
                 celltype.l2 = c(as.character(meta.list[[1]][,4]),
                                 as.character(meta.list[[2]][,4]),
                                 as.character(meta.list[[3]][,4]),
                                 as.character(meta.list[[4]][,4]),
                                 as.character(meta.list[[5]][,4]),
                                 as.character(meta.list[[6]][,4])),
                 modality = c(rep('ABseq',6192),rep('RNA',4035),rep('ADT',9380),
                              rep('ADT',8707),rep('ABseq',5829),rep('RNA',11170))
)

d <- a@Int.result[["bind"]]
d <- d[cellname,]

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
options(future.globals.maxSize = 5000 * 1024^3)
setwd('./Young_Age')
out_dir <- './Benchmarking/CellType_missing/output/Ab-seq'
meta <- readRDS('./Young_Age/meta.rds')

set.seed(123)
rm_ct <- sample(unique(meta$celltype.l2),size = 6,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 6,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 6,replace = T)
rm_ct

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in 1:length(paths)){
  tmp <- readRDS(paste0('./',paths[i],'/meta.rds'))
  meta.list[[paths[i]]] <- dplyr::filter(tmp,celltype.l2 != rm_ct[i])
  adt.list[[paths[i]]] <- readRDS(paste0('./',paths[i],'/adt.rds'))[,rownames(meta.list[[i]])]
  rna.list[[paths[i]]] <- readRDS(paste0('./',paths[i],'/rna.rds'))[,rownames(meta.list[[i]])]
  
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
mm <- data.frame(batch = c(rep("Age1",4811),rep("Age2",3916),rep("Age3",9865),
                           rep("BM1",9595),rep("BM2",6963),rep("BM3",10154)),
                 celltype.l2 = c(as.character(meta.list[[1]][,4]),
                                 as.character(meta.list[[2]][,4]),
                                 as.character(meta.list[[3]][,4]),
                                 as.character(meta.list[[4]][,4]),
                                 as.character(meta.list[[5]][,4]),
                                 as.character(meta.list[[6]][,4])),
                 modality = c(rep('ADT',4811),rep('RNA',3916),rep('ABseq',9865),
                              rep('RNA',9595),rep('ADT',6963),rep('ABseq',10154))
)
d <- a@Int.result[["bind"]]
d <- d[cellname,]

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
options(future.globals.maxSize = 5000 * 1024^3)
setwd('./Young_Age')
out_dir <- './Benchmarking/CellType_missing/output/Ab-seq'
meta <- readRDS('./Young_Age/meta.rds')

set.seed(123)
rm_ct <- sample(unique(meta$celltype.l2),size = 6,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 6,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 6,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 6,replace = T)
rm_ct

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in 1:length(paths)){
  tmp <- readRDS(paste0('./',paths[i],'/meta.rds'))
  meta.list[[paths[i]]] <- dplyr::filter(tmp,celltype.l2 != rm_ct[i])
  adt.list[[paths[i]]] <- readRDS(paste0('./',paths[i],'/adt.rds'))[,rownames(meta.list[[i]])]
  rna.list[[paths[i]]] <- readRDS(paste0('./',paths[i],'/rna.rds'))[,rownames(meta.list[[i]])]
  
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

mm <- data.frame(batch = c(rep("Age1",6259 ),rep("Age2",3312),rep("Age3",10616),
                           rep("BM1",8821),rep("BM2",6924),rep("BM3",11170)),
                 celltype.l2 = c(as.character(meta.list[[1]][,4]),
                                 as.character(meta.list[[2]][,4]),
                                 as.character(meta.list[[3]][,4]),
                                 as.character(meta.list[[4]][,4]),
                                 as.character(meta.list[[5]][,4]),
                                 as.character(meta.list[[6]][,4])),
                 modality = c(rep('RNA',6259 ),rep('ABseq',3312),rep('ADT',10616),
                              rep('ABseq',8821),rep('ADT',6924),rep('RNA',11170))
)

d <- a@Int.result[["bind"]]
d <- d[cellname,]

write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_s3.csv'),
          row.names = FALSE)
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
options(future.globals.maxSize = 5000 * 1024^3)
setwd('./Young_Age')
out_dir <- './Benchmarking/CellType_missing/output/Ab-seq'
meta <- readRDS('./Young_Age/meta.rds')
set.seed(123)
rm_ct <- sample(unique(meta$celltype.l2),size = 6,replace = T)
rm_ct

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in 1:length(paths)){
  tmp <- readRDS(paste0('./',paths[i],'/meta.rds'))
  meta.list[[paths[i]]] <- dplyr::filter(tmp,celltype.l2 != rm_ct[i])
  adt.list[[paths[i]]] <- readRDS(paste0('./',paths[i],'/adt.rds'))[,rownames(meta.list[[i]])]
  rna.list[[paths[i]]] <- readRDS(paste0('./',paths[i],'/rna.rds'))[,rownames(meta.list[[i]])]
  
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

mm <- data.frame(batch = c(rep("Age1",6256),rep("Age2",3949),rep("Age3",10350),
                           rep("BM1",8821),rep("BM2",6382),rep("BM3",10337)),
                 celltype.l2 = c(as.character(meta.list[[1]][,4]),
                                 as.character(meta.list[[2]][,4]),
                                 as.character(meta.list[[3]][,4]),
                                 as.character(meta.list[[4]][,4]),
                                 as.character(meta.list[[5]][,4]),
                                 as.character(meta.list[[6]][,4])),
                 modality = c(rep('ABseq',6256),rep('ADT',3949),rep('RNA',10350),
                              rep('ABseq',8821),rep('RNA',6382),rep('ADT',10337))
)
d <- a@Int.result[["bind"]]
d <- d[cellname,]

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
options(future.globals.maxSize = 5000 * 1024^3)
setwd('./Young_Age')
out_dir <- './Benchmarking/CellType_missing/output/Ab-seq'
meta <- readRDS('./Young_Age/meta.rds')

set.seed(123)
rm_ct <- sample(unique(meta$celltype.l2),size = 6,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 6,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 6,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 6,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 6,replace = T)
rm_ct

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in 1:length(paths)){
  tmp <- readRDS(paste0('./',paths[i],'/meta.rds'))
  meta.list[[paths[i]]] <- dplyr::filter(tmp,celltype.l2 != rm_ct[i])
  adt.list[[paths[i]]] <- readRDS(paste0('./',paths[i],'/adt.rds'))[,rownames(meta.list[[i]])]
  rna.list[[paths[i]]] <- readRDS(paste0('./',paths[i],'/rna.rds'))[,rownames(meta.list[[i]])]
  
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
mm <- data.frame(batch = c(rep("Age1",4793 ),rep("Age2",3916),rep("Age3",8488),
                           rep("BM1",9082),rep("BM2",6295),rep("BM3",10165 )),
                 celltype.l2 = c(as.character(meta.list[[1]][,4]),
                                 as.character(meta.list[[2]][,4]),
                                 as.character(meta.list[[3]][,4]),
                                 as.character(meta.list[[4]][,4]),
                                 as.character(meta.list[[5]][,4]),
                                 as.character(meta.list[[6]][,4])),
                 modality = c(rep('ABseq',4793),rep('ADT',3916),rep('RNA',8488),
                              rep('ADT',9082),rep('ABseq',6295),rep('RNA',10165 ))
)

d <- a@Int.result[["bind"]]
d <- d[cellname,]

write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_s4.csv'),
          row.names = FALSE)
saveRDS(mm,paste0(out_dir,'/meta_s4.rds'))



