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
meta <- readRDS('./Benchmarking/Unsupervised/BMMC_s2/output/random1/meta_34.rds')
setwd('./BMMC_CITE_Multiome')
out_dir <- './Benchmarking/CellType_missing/output/BMMC_s2'

set.seed(123)
rm_ct <- sample(unique(meta$celltype.l2),size = 10,replace = T)
rm_ct

path <- c('CITE/s1d1_cite','CITE/s1d2_cite','CITE/s1d3_cite',
          'Multiome/s1d1_multi','Multiome/s1d2_multi','Multiome/s1d3_multi')
rna.list <- list()
adt.list <- list()
atac.list <- list()
meta.list <- list()
rna.list <- list()
adt.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:6){
  
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
  rna.list[[path[i]]] <- readRDS(paste0(path[i],'/rna.rds'))
  
  if(i %in% c(1,2,3)){
    adt.list[[path[i]]] <- readRDS(paste0(path[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[path[i]]] <- readRDS(paste0(path[i],'/atac.rds'))
  }
}
names(rna.list) <- NULL
names(adt.list) <- NULL
names(atac.list) <- NULL
colnames(adt.list[[1]]) <- paste0(colnames(adt.list[[1]]),'-adt')
colnames(adt.list[[2]]) <- paste0(colnames(adt.list[[2]]),'-adt')

colnames(atac.list[[2]]) <- paste0(colnames(atac.list[[2]]),'-atac')
colnames(atac.list[[3]]) <- paste0(colnames(atac.list[[3]]),'-atac')

colnames(rna.list[[1]]) <- paste0(colnames(rna.list[[1]]),'-rna')
colnames(rna.list[[2]]) <- paste0(colnames(rna.list[[2]]),'-rna')
colnames(rna.list[[5]]) <- paste0(colnames(rna.list[[5]]),'-rna')
colnames(rna.list[[6]]) <- paste0(colnames(rna.list[[6]]),'-rna')


meta_c1rna <- meta.list[[1]]
meta_c1rna$modality <- 'RNA'
meta_c1rna$batch <- 'CITE1_RNA'

meta_c1adt <- meta.list[[1]]
meta_c1adt$modality <- 'ADT'
meta_c1adt$batch <- 'CITE1_ADT'

meta_c2rna <- meta.list[[2]]
meta_c2rna$modality <- 'RNA'
meta_c2rna$batch <- 'CITE2_RNA'

meta_c2adt <- meta.list[[2]]
meta_c2adt$modality <- 'ADT'
meta_c2adt$batch <- 'CITE2_ADT'

meta_c3 <- meta.list[[3]]
meta_c3$modality <- 'CITE'
meta_c3$batch <- 'CITE3'


meta_m1 <- meta.list[[4]]
meta_m1$modality <- 'Multiome'
meta_m1$batch <- 'Multiome1'


meta_m2rna <- meta.list[[5]]
meta_m2rna$modality <- 'RNA'
meta_m2rna$batch <- 'Multiome2_RNA'

meta_m2atac <- meta.list[[5]]
meta_m2atac$modality <- 'ATAC'
meta_m2atac$batch <- 'Multiome2_ATAC'

meta_m3rna <- meta.list[[6]]
meta_m3rna$modality <- 'RNA'
meta_m3rna$batch <- 'Multiome3_RNA'

meta_m3atac <- meta.list[[6]]
meta_m3atac$modality <- 'ATAC'
meta_m3atac$batch <- 'Multiome3_ATAC'
rownames(meta_c1rna) <- colnames(rna.list[[1]])
rownames(meta_c1adt) <- colnames(adt.list[[1]])

rownames(meta_c2rna) <- colnames(rna.list[[2]])
rownames(meta_c2adt) <- colnames(rna.list[[2]])

rownames(meta_c3) <- colnames(rna.list[[3]])

rownames(meta_m1) <- colnames(rna.list[[4]])
rownames(meta_m2atac) <- colnames(atac.list[[2]])
rownames(meta_m2rna) <- colnames(rna.list[[5]])
rownames(meta_m3atac) <- colnames(atac.list[[3]])
rownames(meta_m3rna) <- colnames(rna.list[[6]])
mm <- rbind(meta_c1rna,meta_c1adt,
            meta_c2rna,meta_c2adt,
            meta_c3,
            meta_m1,
            meta_m2rna,meta_m2atac,
            meta_m3rna,meta_m3atac)
meta.list <- list(meta_c1rna,meta_c1adt,
                  meta_c2rna,meta_c2adt,
                  meta_c3,
                  meta_m1,
                  meta_m2rna,meta_m2atac,
                  meta_m3rna,meta_m3atac)
for(i in 1:10){
  meta.list[[i]] <- dplyr::filter(meta.list[[i]],celltype.l2 != rm_ct[i]) 
}
rna.list[[1]] <- rna.list[[1]][,rownames(meta.list[[1]])]

adt.list[[1]] <- adt.list[[1]][,rownames(meta.list[[2]])]

rna.list[[2]] <- rna.list[[2]][,rownames(meta.list[[3]])]

adt.list[[2]] <- adt.list[[2]][,rownames(meta.list[[4]])]

rna.list[[3]] <- rna.list[[3]][,rownames(meta.list[[5]])]
adt.list[[3]] <- adt.list[[3]][,rownames(meta.list[[5]])]

rna.list[[4]] <- rna.list[[4]][,rownames(meta.list[[6]])]
atac.list[[1]] <- atac.list[[1]][,rownames(meta.list[[6]])]

rna.list[[5]] <- rna.list[[5]][,rownames(meta.list[[7]])]

atac.list[[2]] <- atac.list[[2]][,rownames(meta.list[[8]])]

rna.list[[6]] <- rna.list[[6]][,rownames(meta.list[[9]])]

atac.list[[3]] <- atac.list[[3]][,rownames(meta.list[[10]])]

musa_obj <- Create.Musa.Object(data.list = c(rna.list[1:3],
                                             adt.list,
                                             rna.list[4:6],
                                             atac.list), 
                               samples = c('CITE1rna','CITE2rna','CITE3',
                                           'CITE1adt','CITE2adt','CITE3',
                                           'Multiome1','Multiome2rna','Multiome3rna',
                                           'Multiome1','Multiome2atac','Multiome3atac'), 
                               modals = c("rna","rna","rna",
                                          "adt",'adt','adt',
                                          "rna","rna",'rna',
                                          "atac","atac",'atac'),
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
  x <- as.character(rownames(x))
  x
}))

meta <- Reduce(rbind,meta.list)
d <- a@Int.result[["bind"]]
d <- d[cellname,]

write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_34.csv'),
          row.names = FALSE)
saveRDS(meta,paste0(out_dir,'/meta_34.rds'))


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
meta <- readRDS('./Benchmarking/Unsupervised/BMMC_s2/output/random1/meta_34.rds')
setwd('./BMMC_CITE_Multiome')
out_dir <- './Benchmarking/CellType_missing/output/BMMC_s2'
set.seed(123)
rm_ct <- sample(unique(meta$celltype.l2),size = 10,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 10,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 10,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 10,replace = T)
rm_ct

path <- c('CITE/s1d1_cite','CITE/s1d2_cite','CITE/s1d3_cite',
          'Multiome/s1d1_multi','Multiome/s1d2_multi','Multiome/s1d3_multi')
rna.list <- list()
adt.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:6){
  
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
  rna.list[[path[i]]] <- readRDS(paste0(path[i],'/rna.rds'))
  
  if(i %in% c(1,2,3)){
    adt.list[[path[i]]] <- readRDS(paste0(path[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[path[i]]] <- readRDS(paste0(path[i],'/atac.rds'))
  }
}
names(rna.list) <- NULL
names(adt.list) <- NULL
names(atac.list) <- NULL

colnames(adt.list[[1]]) <- paste0(colnames(adt.list[[1]]),'-adt')
colnames(adt.list[[2]]) <- paste0(colnames(adt.list[[2]]),'-adt')

colnames(atac.list[[1]]) <- paste0(colnames(atac.list[[1]]),'-atac')
colnames(atac.list[[3]]) <- paste0(colnames(atac.list[[3]]),'-atac')

colnames(rna.list[[1]]) <- paste0(colnames(rna.list[[1]]),'-rna')
colnames(rna.list[[2]]) <- paste0(colnames(rna.list[[2]]),'-rna')
colnames(rna.list[[4]]) <- paste0(colnames(rna.list[[4]]),'-rna')
colnames(rna.list[[6]]) <- paste0(colnames(rna.list[[6]]),'-rna')

meta_c1rna <- meta.list[[1]]
meta_c1rna$modality <- 'RNA'
meta_c1rna$batch <- 'CITE1_RNA'

meta_c1adt <- meta.list[[1]]
meta_c1adt$modality <- 'ADT'
meta_c1adt$batch <- 'CITE1_ADT'

meta_c2rna <- meta.list[[2]]
meta_c2rna$modality <- 'RNA'
meta_c2rna$batch <- 'CITE2_RNA'

meta_c2adt <- meta.list[[2]]
meta_c2adt$modality <- 'ADT'
meta_c2adt$batch <- 'CITE2_ADT'

meta_c3 <- meta.list[[3]]
meta_c3$modality <- 'CITE'
meta_c3$batch <- 'CITE3'

meta_m1rna <- meta.list[[4]]
meta_m1rna$modality <- 'RNA'
meta_m1rna$batch <- 'Multiome1_RNA'

meta_m1atac <- meta.list[[4]]
meta_m1atac$modality <- 'ATAC'
meta_m1atac$batch <- 'Multiome1_ATAC'

meta_m2 <- meta.list[[5]]
meta_m2$modality <- 'Multiome'
meta_m2$batch <- 'Multiome2'

meta_m3rna <- meta.list[[6]]
meta_m3rna$modality <- 'RNA'
meta_m3rna$batch <- 'Multiome3_RNA'

meta_m3atac <- meta.list[[6]]
meta_m3atac$modality <- 'ATAC'
meta_m3atac$batch <- 'Multiome3_ATAC'

rownames(meta_c1rna) <- colnames(rna.list[[1]])
rownames(meta_c1adt) <- colnames(adt.list[[1]])
rownames(meta_c2rna) <- colnames(rna.list[[2]])
rownames(meta_c2adt) <- colnames(adt.list[[2]])
rownames(meta_c3) <- colnames(rna.list[[3]])

rownames(meta_m1rna) <- colnames(rna.list[[4]])
rownames(meta_m1atac) <- colnames(atac.list[[1]])
rownames(meta_m2) <- colnames(rna.list[[5]])
rownames(meta_m3rna) <- colnames(rna.list[[6]])
rownames(meta_m3atac) <- colnames(atac.list[[3]])

meta.list <- list(meta_c1rna,meta_c1adt,
                  meta_c2rna,meta_c2adt,
                  meta_c3,
                  meta_m1rna,meta_m1atac,
                  meta_m2,
                  meta_m3rna,meta_m3atac)
for(i in 1:10){
  meta.list[[i]] <- dplyr::filter(meta.list[[i]],celltype.l2 != rm_ct[i]) 
}

rna.list[[1]] <- rna.list[[1]][,rownames(meta.list[[1]])]

adt.list[[1]] <- adt.list[[1]][,rownames(meta.list[[2]])]

rna.list[[2]] <- rna.list[[2]][,rownames(meta.list[[3]])]

adt.list[[2]] <- adt.list[[2]][,rownames(meta.list[[4]])]

rna.list[[3]] <- rna.list[[3]][,rownames(meta.list[[5]])]
adt.list[[3]] <- adt.list[[3]][,rownames(meta.list[[5]])]

rna.list[[4]] <- rna.list[[4]][,rownames(meta.list[[6]])]

atac.list[[1]] <- atac.list[[1]][,rownames(meta.list[[7]])]

rna.list[[5]] <- rna.list[[5]][,rownames(meta.list[[8]])]
atac.list[[2]] <- atac.list[[2]][,rownames(meta.list[[8]])]

rna.list[[6]] <- rna.list[[6]][,rownames(meta.list[[9]])]

atac.list[[3]] <- atac.list[[3]][,rownames(meta.list[[10]])]

musa_obj <- Create.Musa.Object(data.list = c(rna.list,
                                             adt.list,
                                             atac.list), 
                               samples = c('cite1rna','cite2rna','cite3','mul1rna','mul2','mul3rna',
                                           'cite1adt','cite2adt','cite3','mul1atac','mul2','mul3atac'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna",'rna',
                                          "adt",'adt','adt',
                                          "atac","atac",'atac'),
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
  x <- as.character(rownames(x))
  x
}))
meta <- Reduce(rbind,meta.list)
d <- a@Int.result[["bind"]]
d <- d[cellname,]

write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_35.csv'),
          row.names = FALSE)
saveRDS(meta,paste0(out_dir,'/meta_35.rds'))


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
meta <- readRDS('./Benchmarking/Unsupervised/BMMC_s2/output/random1/meta_34.rds')
setwd('./BMMC_CITE_Multiome')
out_dir <- './Benchmarking/CellType_missing/output/BMMC_s2'

set.seed(123)
rm_ct <- sample(unique(meta$celltype.l2),size = 10,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 10,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 10,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 10,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 10,replace = T)
rm_ct

path <- c('CITE/s1d1_cite','CITE/s1d2_cite','CITE/s1d3_cite',
          'Multiome/s1d1_multi','Multiome/s1d2_multi','Multiome/s1d3_multi')
rna.list <- list()
adt.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:6){
  
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
  rna.list[[path[i]]] <- readRDS(paste0(path[i],'/rna.rds'))
  
  if(i %in% c(1,2,3)){
    adt.list[[path[i]]] <- readRDS(paste0(path[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[path[i]]] <- readRDS(paste0(path[i],'/atac.rds'))
  }
}
names(rna.list) <- NULL
names(adt.list) <- NULL
names(atac.list) <- NULL

colnames(adt.list[[1]]) <- paste0(colnames(adt.list[[1]]),'-adt')
colnames(adt.list[[2]]) <- paste0(colnames(adt.list[[2]]),'-adt')

colnames(atac.list[[1]]) <- paste0(colnames(atac.list[[1]]),'-atac')
colnames(atac.list[[2]]) <- paste0(colnames(atac.list[[2]]),'-atac')

colnames(rna.list[[1]]) <- paste0(colnames(rna.list[[1]]),'-rna')
colnames(rna.list[[2]]) <- paste0(colnames(rna.list[[2]]),'-rna')
colnames(rna.list[[4]]) <- paste0(colnames(rna.list[[4]]),'-rna')
colnames(rna.list[[5]]) <- paste0(colnames(rna.list[[5]]),'-rna')

meta_c1rna <- meta.list[[1]]
meta_c1rna$modality <- 'RNA'
meta_c1rna$batch <- 'CITE1_RNA'

meta_c1adt <- meta.list[[1]]
meta_c1adt$modality <- 'ADT'
meta_c1adt$batch <- 'CITE1_ADT'

meta_c2rna <- meta.list[[2]]
meta_c2rna$modality <- 'RNA'
meta_c2rna$batch <- 'CITE2_RNA'

meta_c2adt <- meta.list[[2]]
meta_c2adt$modality <- 'ADT'
meta_c2adt$batch <- 'CITE2_ADT'

meta_c3 <- meta.list[[3]]
meta_c3$modality <- 'CITE'
meta_c3$batch <- 'CITE3'

meta_m1rna <- meta.list[[4]]
meta_m1rna$modality <- 'RNA'
meta_m1rna$batch <- 'Multiome1_RNA'

meta_m1atac <- meta.list[[4]]
meta_m1atac$modality <- 'ATAC'
meta_m1atac$batch <- 'Multiome1_ATAC'

meta_m2rna <- meta.list[[5]]
meta_m2rna$modality <- 'RNA'
meta_m2rna$batch <- 'Multiome2_RNA'

meta_m2atac <- meta.list[[5]]
meta_m2atac$modality <- 'ATAC'
meta_m2atac$batch <- 'Multiome2_ATAC'

meta_m3 <- meta.list[[6]]
meta_m3$modality <- 'Multiome'
meta_m3$batch <- 'Multiome3'

rownames(meta_c1rna) <- colnames(rna.list[[1]])
rownames(meta_c1adt) <- colnames(adt.list[[1]])
rownames(meta_c2rna) <- colnames(rna.list[[2]])
rownames(meta_c2adt) <- colnames(adt.list[[2]])
rownames(meta_c3) <- colnames(rna.list[[3]])

rownames(meta_m1rna) <- colnames(rna.list[[4]])
rownames(meta_m1atac) <- colnames(atac.list[[1]])
rownames(meta_m2rna) <- colnames(rna.list[[5]])
rownames(meta_m2atac) <- colnames(atac.list[[2]])
rownames(meta_m3) <- colnames(rna.list[[6]])

meta.list <- list(meta_c1rna,meta_c1adt,
                  meta_c2rna,meta_c2adt,
                  meta_c3,
                  meta_m1rna,meta_m1atac,
                  meta_m2rna,meta_m2atac,
                  meta_m3)

for(i in 1:10){
  meta.list[[i]] <- dplyr::filter(meta.list[[i]],celltype.l2 != rm_ct[i]) 
}

rna.list[[1]] <- rna.list[[1]][,rownames(meta.list[[1]])]

adt.list[[1]] <- adt.list[[1]][,rownames(meta.list[[2]])]

rna.list[[2]] <- rna.list[[2]][,rownames(meta.list[[3]])]

adt.list[[2]] <- adt.list[[2]][,rownames(meta.list[[4]])]

rna.list[[3]] <- rna.list[[3]][,rownames(meta.list[[5]])]
adt.list[[3]] <- adt.list[[3]][,rownames(meta.list[[5]])]

rna.list[[4]] <- rna.list[[4]][,rownames(meta.list[[6]])]

atac.list[[1]] <- atac.list[[1]][,rownames(meta.list[[7]])]

rna.list[[5]] <- rna.list[[5]][,rownames(meta.list[[8]])]

atac.list[[2]] <- atac.list[[2]][,rownames(meta.list[[9]])]

rna.list[[6]] <- rna.list[[6]][,rownames(meta.list[[10]])]
atac.list[[3]] <- atac.list[[3]][,rownames(meta.list[[10]])]

musa_obj <- Create.Musa.Object(data.list = c(rna.list,
                                             adt.list,
                                             atac.list), 
                               samples = c('cite1rna','cite2rna','cite3','mul1rna','mul2rna','mul3',
                                           'cite1adt','cite2adt','cite3','mul1atac','mul2atac','mul3'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna",'rna',
                                          "adt",'adt','adt',
                                          "atac","atac",'atac'),
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
  x <- as.character(rownames(x))
  x
}))
meta <- Reduce(rbind,meta.list)
d <- a@Int.result[["bind"]]
d <- d[cellname,]

write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_36.csv'),
          row.names = FALSE)
saveRDS(meta,paste0(out_dir,'/meta_36.rds'))


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
meta <- readRDS('./Benchmarking/Unsupervised/BMMC_s2/output/random1/meta_34.rds')
setwd('./BMMC_CITE_Multiome')
out_dir <- './Benchmarking/CellType_missing/output/BMMC_s2'

set.seed(123)
rm_ct <- sample(unique(meta$celltype.l2),size = 10,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 10,replace = T)
rm_ct

path <- c('CITE/s1d1_cite','CITE/s1d2_cite','CITE/s1d3_cite',
          'Multiome/s1d1_multi','Multiome/s1d2_multi','Multiome/s1d3_multi')
rna.list <- list()
adt.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:6){
  
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
  rna.list[[path[i]]] <- readRDS(paste0(path[i],'/rna.rds'))
  
  if(i %in% c(1,2,3)){
    adt.list[[path[i]]] <- readRDS(paste0(path[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[path[i]]] <- readRDS(paste0(path[i],'/atac.rds'))
  }
}
names(rna.list) <- NULL
names(adt.list) <- NULL
names(atac.list) <- NULL

colnames(adt.list[[1]]) <- paste0(colnames(adt.list[[1]]),'-adt')
colnames(adt.list[[3]]) <- paste0(colnames(adt.list[[3]]),'-adt')

colnames(atac.list[[2]]) <- paste0(colnames(atac.list[[2]]),'-atac')
colnames(atac.list[[3]]) <- paste0(colnames(atac.list[[3]]),'-atac')

colnames(rna.list[[1]]) <- paste0(colnames(rna.list[[1]]),'-rna')
colnames(rna.list[[3]]) <- paste0(colnames(rna.list[[3]]),'-rna')
colnames(rna.list[[5]]) <- paste0(colnames(rna.list[[5]]),'-rna')
colnames(rna.list[[6]]) <- paste0(colnames(rna.list[[6]]),'-rna')
meta_c1rna <- meta.list[[1]]
meta_c1rna$modality <- 'RNA'
meta_c1rna$batch <- 'CITE1_RNA'

meta_c1adt <- meta.list[[1]]
meta_c1adt$modality <- 'ADT'
meta_c1adt$batch <- 'CITE1_ADT'

meta_c2 <- meta.list[[2]]
meta_c2$modality <- 'CITE'
meta_c2$batch <- 'CITE2'

meta_c3rna <- meta.list[[3]]
meta_c3rna$modality <- 'RNA'
meta_c3rna$batch <- 'CITE3_RNA'

meta_c3adt <- meta.list[[3]]
meta_c3adt$modality <- 'ADT'
meta_c3adt$batch <- 'CITE3_ADT'


meta_m1 <- meta.list[[4]]
meta_m1$modality <- 'Multiome'
meta_m1$batch <- 'Multiome1'


meta_m2rna <- meta.list[[5]]
meta_m2rna$modality <- 'RNA'
meta_m2rna$batch <- 'Multiome2_RNA'

meta_m2atac <- meta.list[[5]]
meta_m2atac$modality <- 'ATAC'
meta_m2atac$batch <- 'Multiome2_ATAC'

meta_m3rna <- meta.list[[6]]
meta_m3rna$modality <- 'RNA'
meta_m3rna$batch <- 'Multiome3_RNA'

meta_m3atac <- meta.list[[6]]
meta_m3atac$modality <- 'ATAC'
meta_m3atac$batch <- 'Multiome3_ATAC'

rownames(meta_c1rna) <- colnames(rna.list[[1]])
rownames(meta_c1adt) <- colnames(adt.list[[1]])
rownames(meta_c2) <- colnames(rna.list[[2]])
rownames(meta_c3rna) <- colnames(rna.list[[3]])
rownames(meta_c3adt) <- colnames(adt.list[[3]])

rownames(meta_m1) <- colnames(rna.list[[4]])
rownames(meta_m2atac) <- colnames(atac.list[[2]])
rownames(meta_m2rna) <- colnames(rna.list[[5]])
rownames(meta_m3atac) <- colnames(atac.list[[3]])
rownames(meta_m3rna) <- colnames(rna.list[[6]])

meta.list <- list(meta_c1rna,meta_c1adt,
                  meta_c2,
                  meta_c3rna,meta_c3adt,
                  meta_m1,
                  meta_m2rna,meta_m2atac,
                  meta_m3rna,meta_m3atac)
for(i in 1:10){
  meta.list[[i]] <- dplyr::filter(meta.list[[i]],celltype.l2 != rm_ct[i]) 
}
rna.list[[1]] <- rna.list[[1]][,rownames(meta.list[[1]])]

adt.list[[1]] <- adt.list[[1]][,rownames(meta.list[[2]])]

rna.list[[2]] <- rna.list[[2]][,rownames(meta.list[[3]])]
adt.list[[2]] <- adt.list[[2]][,rownames(meta.list[[3]])]

rna.list[[3]] <- rna.list[[3]][,rownames(meta.list[[4]])]

adt.list[[3]] <- adt.list[[3]][,rownames(meta.list[[5]])]

rna.list[[4]] <- rna.list[[4]][,rownames(meta.list[[6]])]
atac.list[[1]] <- atac.list[[1]][,rownames(meta.list[[6]])]

rna.list[[5]] <- rna.list[[5]][,rownames(meta.list[[7]])]

atac.list[[2]] <- atac.list[[2]][,rownames(meta.list[[8]])]

rna.list[[6]] <- rna.list[[6]][,rownames(meta.list[[9]])]

atac.list[[3]] <- atac.list[[3]][,rownames(meta.list[[10]])]


musa_obj <- Create.Musa.Object(data.list = c(rna.list,
                                             adt.list,
                                             atac.list), 
                               samples = c('cite1rna','cite2','cite3rna','mul1','mul2rna','mul3rna',
                                           'cite1adt','cite2','cite3','mul1','mul2atac','mul3atac'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna",'rna',
                                          "adt",'adt','adt',
                                          "atac","atac",'atac'),
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
  x <- as.character(rownames(x))
  x
}))

meta <- Reduce(rbind,meta.list)
d <- a@Int.result[["bind"]]
d <- d[cellname,]

write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_24.csv'),
          row.names = FALSE)
saveRDS(meta,paste0(out_dir,'/meta_24.rds'))

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
meta <- readRDS('./Benchmarking/Unsupervised/BMMC_s2/output/random1/meta_34.rds')
setwd('./BMMC_CITE_Multiome')
out_dir <- './Benchmarking/CellType_missing/output/BMMC_s2'
set.seed(123)
rm_ct <- sample(unique(meta$celltype.l2),size = 10,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 10,replace = T)
rm_ct <- sample(unique(meta$celltype.l2),size = 10,replace = T)
rm_ct

path <- c('CITE/s1d1_cite','CITE/s1d2_cite','CITE/s1d3_cite',
          'Multiome/s1d1_multi','Multiome/s1d2_multi','Multiome/s1d3_multi')
rna.list <- list()
adt.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:6){
  
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
  rna.list[[path[i]]] <- readRDS(paste0(path[i],'/rna.rds'))
  
  if(i %in% c(1,2,3)){
    adt.list[[path[i]]] <- readRDS(paste0(path[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[path[i]]] <- readRDS(paste0(path[i],'/atac.rds'))
  }
}
names(rna.list) <- NULL
names(adt.list) <- NULL
names(atac.list) <- NULL
colnames(adt.list[[1]]) <- paste0(colnames(adt.list[[1]]),'-adt')
colnames(adt.list[[3]]) <- paste0(colnames(adt.list[[3]]),'-adt')

colnames(atac.list[[1]]) <- paste0(colnames(atac.list[[1]]),'-atac')
colnames(atac.list[[2]]) <- paste0(colnames(atac.list[[2]]),'-atac')

colnames(rna.list[[1]]) <- paste0(colnames(rna.list[[1]]),'-rna')
colnames(rna.list[[3]]) <- paste0(colnames(rna.list[[3]]),'-rna')
colnames(rna.list[[4]]) <- paste0(colnames(rna.list[[4]]),'-rna')
colnames(rna.list[[5]]) <- paste0(colnames(rna.list[[5]]),'-rna')

meta_c1rna <- meta.list[[1]]
meta_c1rna$modality <- 'RNA'
meta_c1rna$batch <- 'CITE1_RNA'

meta_c1adt <- meta.list[[1]]
meta_c1adt$modality <- 'ADT'
meta_c1adt$batch <- 'CITE1_ADT'

meta_c2 <- meta.list[[2]]
meta_c2$modality <- 'CITE'
meta_c2$batch <- 'CITE2'

meta_c3rna <- meta.list[[3]]
meta_c3rna$modality <- 'RNA'
meta_c3rna$batch <- 'CITE3_RNA'

meta_c3adt <- meta.list[[3]]
meta_c3adt$modality <- 'ADT'
meta_c3adt$batch <- 'CITE3_ADT'


meta_m1rna <- meta.list[[4]]
meta_m1rna$modality <- 'RNA'
meta_m1rna$batch <- 'Multiome1_RNA'

meta_m1atac <- meta.list[[4]]
meta_m1atac$modality <- 'ATAC'
meta_m1atac$batch <- 'Multiome1_ATAC'

meta_m2rna <- meta.list[[5]]
meta_m2rna$modality <- 'RNA'
meta_m2rna$batch <- 'Multiome2_RNA'

meta_m2atac <- meta.list[[5]]
meta_m2atac$modality <- 'ATAC'
meta_m2atac$batch <- 'Multiome2_ATAC'

meta_m3 <- meta.list[[6]]
meta_m3$modality <- 'Multiome'
meta_m3$batch <- 'Multiome3'

rownames(meta_c1rna) <- colnames(rna.list[[1]])
rownames(meta_c1adt) <- colnames(adt.list[[1]])
rownames(meta_c2) <- colnames(rna.list[[2]])
rownames(meta_c3rna) <- colnames(rna.list[[3]])
rownames(meta_c3adt) <- colnames(adt.list[[3]])

rownames(meta_m1rna) <- colnames(rna.list[[4]])
rownames(meta_m1atac) <- colnames(atac.list[[1]])
rownames(meta_m2rna) <- colnames(rna.list[[5]])
rownames(meta_m2atac) <- colnames(atac.list[[2]])
rownames(meta_m3) <- colnames(rna.list[[6]])

meta.list <- list(meta_c1rna,meta_c1adt,
                  meta_c2,
                  meta_c3rna,meta_c3adt,
                  meta_m1rna,meta_m1atac,
                  meta_m2rna,meta_m2atac,
                  meta_m3)

for(i in 1:10){
  meta.list[[i]] <- dplyr::filter(meta.list[[i]],celltype.l2 != rm_ct[i]) 
}

rna.list[[1]] <- rna.list[[1]][,rownames(meta.list[[1]])]

adt.list[[1]] <- adt.list[[1]][,rownames(meta.list[[2]])]

rna.list[[2]] <- rna.list[[2]][,rownames(meta.list[[3]])]
adt.list[[2]] <- adt.list[[2]][,rownames(meta.list[[3]])]

rna.list[[3]] <- rna.list[[3]][,rownames(meta.list[[4]])]

adt.list[[3]] <- adt.list[[3]][,rownames(meta.list[[5]])]

rna.list[[4]] <- rna.list[[4]][,rownames(meta.list[[6]])]

atac.list[[1]] <- atac.list[[1]][,rownames(meta.list[[7]])]

rna.list[[5]] <- rna.list[[5]][,rownames(meta.list[[8]])]

atac.list[[2]] <- atac.list[[2]][,rownames(meta.list[[9]])]

rna.list[[6]] <- rna.list[[6]][,rownames(meta.list[[10]])]
atac.list[[3]] <- atac.list[[3]][,rownames(meta.list[[10]])]

musa_obj <- Create.Musa.Object(data.list = c(rna.list,
                                             adt.list,
                                             atac.list), 
                               samples = c('cite1rna','cite2','cite3rna','mul1rna','mul2rna','mul3',
                                           'cite1adt','cite2','cite3','mul1atac','mul2atac','mul3'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna",'rna',
                                          "adt",'adt','adt',
                                          "atac","atac",'atac'),
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
  x <- as.character(rownames(x))
  x
}))
meta <- Reduce(rbind,meta.list)
d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_26.csv'),
          row.names = FALSE)
saveRDS(meta,paste0(out_dir,'/meta_26.rds'))



