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
meta <- readRDS('./Benchmarking/Unsupervised/BMMC_s1/output/random1/meta_drop5.rds')
setwd('./BMMC_CITE_Multiome')
out_dir <- './Benchmarking/CellType_missing/output/BMMC_s1'

set.seed(123)
rm_ct <- sample(unique(meta$celltype.l2),size = 5,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 5,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 5,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 5,replace = T)
rm_ct

path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:3){
  tmp <- readRDS(paste0(path1[i],'/meta.rds'))
  meta1.list[[i]] <- dplyr::filter(tmp,celltype.l2 != rm_ct[i])
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))[,rownames(meta1.list[[i]])]
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))[,rownames(meta1.list[[i]])]
}
path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:2){
  tmp <- readRDS(paste0(path2[i],'/meta.rds'))
  meta2.list[[i]] <- dplyr::filter(tmp,celltype.l2 != rm_ct[i+3])
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))[,rownames(meta2.list[[i]])]
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))[,rownames(meta2.list[[i]])]
}
meta.list <- c(meta1.list,meta2.list)

musa_obj <- Create.Musa.Object(data.list = c(rna1.list,
                                             rna2.list,
                                             adt.list,
                                             atac.list), 
                               samples = c('CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome3',
                                           'CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome3'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna",
                                          "adt",'adt','adt',
                                          "atac","atac"),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna","atac",'adt'), 
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
                         dims = list(1:50,1:50,1:50),# 使用哪些维数
                         method = c("PCA",'PCA','LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))

mm <- data.frame(batch = c(rep("CITE1",nrow(meta1.list[[1]])),
                           rep("CITE2",nrow(meta1.list[[2]])),
                           rep("CITE3",nrow(meta1.list[[3]])),
                           rep("Multiome1",nrow(meta2.list[[1]])),
                           rep("Multiome3",nrow(meta2.list[[2]]))),
                 celltype.l2 = c(as.character(meta1.list[[1]][,5]),
                                 as.character(meta1.list[[2]][,5]),
                                 as.character(meta1.list[[3]][,5]),
                                 as.character(meta2.list[[1]][,5]),
                                 as.character(meta2.list[[2]][,5])),
                 celltype.l1 = c(as.character(meta1.list[[1]][,4]),
                                 as.character(meta1.list[[2]][,4]),
                                 as.character(meta1.list[[3]][,4]),
                                 as.character(meta2.list[[1]][,4]),
                                 as.character(meta2.list[[2]][,4])),
                 modality = c(rep('CITE-seq',(nrow(meta1.list[[1]])+nrow(meta1.list[[2]])+nrow(meta1.list[[3]]))),
                              rep('Multiome',(nrow(meta2.list[[1]])+nrow(meta2.list[[2]]))))
)
d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_drop5.csv'),
          row.names = FALSE)
saveRDS(mm,paste0(out_dir,'/meta_drop5.rds'))


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
meta <- readRDS('./Benchmarking/Unsupervised/BMMC_s1/output/random2/meta_drop1.rds')
setwd('./BMMC_CITE_Multiome')
out_dir <- './Benchmarking/CellType_missing/output/BMMC_s1'

set.seed(123)
rm_ct <- sample(unique(meta$celltype.l2),size = 5,replace = T)
rm_ct
path1 <- c('./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:2){
  tmp <- readRDS(paste0(path1[i],'/meta.rds'))
  meta1.list[[i]] <- dplyr::filter(tmp,celltype.l2 != rm_ct[i])
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))[,rownames(meta1.list[[i]])]
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))[,rownames(meta1.list[[i]])]
}
path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:3){
  tmp <- readRDS(paste0(path2[i],'/meta.rds'))
  meta2.list[[i]] <- dplyr::filter(tmp,celltype.l2 != rm_ct[i+2])
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))[,rownames(meta2.list[[i]])]
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))[,rownames(meta2.list[[i]])]
}
meta.list <- c(meta1.list,meta2.list)

musa_obj <- Create.Musa.Object(data.list = c(rna1.list,
                                             rna2.list,
                                             adt.list,
                                             atac.list), 
                               samples = c('CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3'), 
                               modals = c("rna","rna",
                                          "rna","rna","rna",
                                          'adt','adt',
                                          "atac","atac","atac"),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna","atac",'adt'), 
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
                         dims = list(1:50,1:50,1:50),# 使用哪些维数
                         method = c("PCA",'PCA','LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))
mm <- data.frame(batch = c(rep("CITE1",nrow(meta1.list[[1]])),
                           rep("CITE2",nrow(meta1.list[[2]])),
                           rep("Multiome1",nrow(meta2.list[[1]])),
                           rep("Multiome2",nrow(meta2.list[[2]])),
                           rep("Multiome3",nrow(meta2.list[[3]]))),
                 celltype.l2 = c(as.character(meta1.list[[1]][,5]),
                                 as.character(meta1.list[[2]][,5]),
                                 as.character(meta2.list[[1]][,5]),
                                 as.character(meta2.list[[2]][,5]),
                                 as.character(meta2.list[[3]][,5])),
                 celltype.l1 = c(as.character(meta1.list[[1]][,4]),
                                 as.character(meta1.list[[2]][,4]),
                                 as.character(meta2.list[[1]][,4]),
                                 as.character(meta2.list[[2]][,4]),
                                 as.character(meta2.list[[3]][,4])),
                 modality = c(rep('CITE-seq',(nrow(meta1.list[[1]])+nrow(meta1.list[[2]]))),
                              rep('Multiome',(nrow(meta2.list[[1]])+nrow(meta2.list[[2]])+nrow(meta2.list[[3]]))))
)
d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_drop1.csv'),
          row.names = FALSE)
saveRDS(mm,paste0(out_dir,'/meta_drop1.rds'))


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
meta <- readRDS('./Benchmarking/Unsupervised/BMMC_s1/output/random3/meta_drop6.rds')
setwd('./BMMC_CITE_Multiome')
out_dir <- './Benchmarking/CellType_missing/output/BMMC_s1'

set.seed(123)
rm_ct <- sample(unique(meta$celltype.l2),size = 5,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 5,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 5,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 5,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 5,replace = T)
rm_ct

path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:3){
  tmp <- readRDS(paste0(path1[i],'/meta.rds'))
  meta1.list[[i]] <- dplyr::filter(tmp,celltype.l2 != rm_ct[i])
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))[,rownames(meta1.list[[i]])]
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))[,rownames(meta1.list[[i]])]
}
path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:2){
  tmp <- readRDS(paste0(path2[i],'/meta.rds'))
  meta2.list[[i]] <- dplyr::filter(tmp,celltype.l2 != rm_ct[i+3])
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))[,rownames(meta2.list[[i]])]
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))[,rownames(meta2.list[[i]])]
}
meta.list <- c(meta1.list,meta2.list)

musa_obj <- Create.Musa.Object(data.list = c(rna1.list,
                                             rna2.list,
                                             adt.list,
                                             atac.list), 
                               samples = c('CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2',
                                           'CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna",
                                          "adt",'adt','adt',
                                          "atac","atac"),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna","atac",'adt'), 
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
                         dims = list(1:50,1:50,1:50),# 使用哪些维数
                         method = c("PCA",'PCA','LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)

cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))
mm <- data.frame(batch = c(rep("CITE1",nrow(meta1.list[[1]])),
                           rep("CITE2",nrow(meta1.list[[2]])),
                           rep("CITE3",nrow(meta1.list[[3]])),
                           rep("Multiome1",nrow(meta2.list[[1]])),
                           rep("Multiome3",nrow(meta2.list[[2]]))),
                 celltype.l2 = c(as.character(meta1.list[[1]][,5]),
                                 as.character(meta1.list[[2]][,5]),
                                 as.character(meta1.list[[3]][,5]),
                                 as.character(meta2.list[[1]][,5]),
                                 as.character(meta2.list[[2]][,5])),
                 celltype.l1 = c(as.character(meta1.list[[1]][,4]),
                                 as.character(meta1.list[[2]][,4]),
                                 as.character(meta1.list[[3]][,4]),
                                 as.character(meta2.list[[1]][,4]),
                                 as.character(meta2.list[[2]][,4])),
                 modality = c(rep('CITE-seq',(nrow(meta1.list[[1]])+nrow(meta1.list[[2]])+nrow(meta1.list[[3]]))),
                              rep('Multiome',(nrow(meta2.list[[1]])+nrow(meta2.list[[2]]))))
)
d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_drop6.csv'),
          row.names = FALSE)
saveRDS(mm,paste0(out_dir,'/meta_drop6.rds'))


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
meta <- readRDS('./Benchmarking/Unsupervised/BMMC_s1/output/random4/meta_drop4.rds')
setwd('./BMMC_CITE_Multiome')
out_dir <- './Benchmarking/CellType_missing/output/BMMC_s1'

set.seed(123)
rm_ct <- sample(unique(meta$celltype.l2),size = 5,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 5,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 5,replace = T)
rm_ct

path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:3){
  tmp <- readRDS(paste0(path1[i],'/meta.rds'))
  meta1.list[[i]] <- dplyr::filter(tmp,celltype.l2 != rm_ct[i])
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))[,rownames(meta1.list[[i]])]
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))[,rownames(meta1.list[[i]])]
}
path2 <- c('./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:2){
  tmp <- readRDS(paste0(path2[i],'/meta.rds'))
  meta2.list[[i]] <- dplyr::filter(tmp,celltype.l2 != rm_ct[i+3])
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))[,rownames(meta2.list[[i]])]
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))[,rownames(meta2.list[[i]])]
}
meta.list <- c(meta1.list,meta2.list)

musa_obj <- Create.Musa.Object(data.list = c(rna1.list,
                                             rna2.list,
                                             adt.list,
                                             atac.list), 
                               samples = c('CITE1','CITE2','CITE3',
                                           'Multiome2','Multiome3',
                                           'CITE1','CITE2','CITE3',
                                           'Multiome2','Multiome3'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna",
                                          "adt",'adt','adt',
                                          "atac","atac"),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna","atac",'adt'), 
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
                         dims = list(1:50,1:50,1:50),# 使用哪些维数
                         method = c("PCA",'PCA','LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))
mm <- data.frame(batch = c(rep("CITE1",nrow(meta1.list[[1]])),
                           rep("CITE2",nrow(meta1.list[[2]])),
                           rep("CITE3",nrow(meta1.list[[3]])),
                           rep("Multiome1",nrow(meta2.list[[1]])),
                           rep("Multiome3",nrow(meta2.list[[2]]))),
                 celltype.l2 = c(as.character(meta1.list[[1]][,5]),
                                 as.character(meta1.list[[2]][,5]),
                                 as.character(meta1.list[[3]][,5]),
                                 as.character(meta2.list[[1]][,5]),
                                 as.character(meta2.list[[2]][,5])),
                 celltype.l1 = c(as.character(meta1.list[[1]][,4]),
                                 as.character(meta1.list[[2]][,4]),
                                 as.character(meta1.list[[3]][,4]),
                                 as.character(meta2.list[[1]][,4]),
                                 as.character(meta2.list[[2]][,4])),
                 modality = c(rep('CITE-seq',(nrow(meta1.list[[1]])+nrow(meta1.list[[2]])+nrow(meta1.list[[3]]))),
                              rep('Multiome',(nrow(meta2.list[[1]])+nrow(meta2.list[[2]]))))
)
d <- a@Int.result[["bind"]]
d <- d[cellname,]

write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_drop4.csv'),
          row.names = FALSE)
saveRDS(mm,paste0(out_dir,'/meta_drop4.rds'))

######################### Random 5 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
options(future.globals.maxSize = 5000 * 1024^3)
meta <- readRDS('./Benchmarking/Unsupervised/BMMC_s1/output/random5/meta_drop2.rds')
setwd('./BMMC_CITE_Multiome')
out_dir <- './Benchmarking/CellType_missing/output/BMMC_s1'
set.seed(123)
rm_ct <- sample(unique(meta$celltype.l2),size = 5,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 5,replace = T)
rm_ct

path1 <- c('./CITE/s1d1_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:2){
  tmp <- readRDS(paste0(path1[i],'/meta.rds'))
  meta1.list[[i]] <- dplyr::filter(tmp,celltype.l2 != rm_ct[i])
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))[,rownames(meta1.list[[i]])]
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))[,rownames(meta1.list[[i]])]
}

path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:3){
  tmp <- readRDS(paste0(path2[i],'/meta.rds'))
  meta2.list[[i]] <- dplyr::filter(tmp,celltype.l2 != rm_ct[i+2])
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))[,rownames(meta2.list[[i]])]
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))[,rownames(meta2.list[[i]])]
}
meta.list <- c(meta1.list,meta2.list)

musa_obj <- Create.Musa.Object(data.list = c(rna1.list,
                                             rna2.list,
                                             adt.list,
                                             atac.list), 
                               samples = c('CITE1','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'CITE1','CITE3',
                                           'Multiome1','Multiome2','Multiome3'), 
                               modals = c("rna","rna",
                                          "rna","rna","rna",
                                          "adt",'adt',
                                          "atac","atac","atac"),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna","atac",'adt'), 
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
                         dims = list(1:50,1:50,1:50),# 使用哪些维数
                         method = c("PCA",'PCA','LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))

mm <- data.frame(batch = c(rep("CITE1",nrow(meta1.list[[1]])),
                           rep("CITE2",nrow(meta1.list[[2]])),
                           rep("Multiome1",nrow(meta2.list[[1]])),
                           rep("Multiome2",nrow(meta2.list[[2]])),
                           rep("Multiome3",nrow(meta2.list[[3]]))),
                 celltype.l2 = c(as.character(meta1.list[[1]][,5]),
                                 as.character(meta1.list[[2]][,5]),
                                 as.character(meta2.list[[1]][,5]),
                                 as.character(meta2.list[[2]][,5]),
                                 as.character(meta2.list[[3]][,5])),
                 celltype.l1 = c(as.character(meta1.list[[1]][,4]),
                                 as.character(meta1.list[[2]][,4]),
                                 as.character(meta2.list[[1]][,4]),
                                 as.character(meta2.list[[2]][,4]),
                                 as.character(meta2.list[[3]][,4])),
                 modality = c(rep('CITE-seq',(nrow(meta1.list[[1]])+nrow(meta1.list[[2]]))),
                              rep('Multiome',(nrow(meta2.list[[1]])+nrow(meta2.list[[2]])+nrow(meta2.list[[3]]))))
)
d <- a@Int.result[["bind"]]
d <- d[cellname,]

write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_drop2.csv'),
          row.names = FALSE)
saveRDS(mm,paste0(out_dir,'/meta_drop2.rds'))

