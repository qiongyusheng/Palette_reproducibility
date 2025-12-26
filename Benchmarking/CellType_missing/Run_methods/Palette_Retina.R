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
setwd('./Human_Retina')
out_dir <- './Benchmarking/CellType_missing/output/Retina'
meta <- readRDS('./Human_Retina/meta.rds')

set.seed(123)
rm_ct <- sample(unique(meta$celltype),size = 8,replace = T)
rm_ct <- sample(unique(meta$celltype),size = 8,replace = T)
rm_ct

path <- c('./LGS1OD','./LGS1OS','./LGS2OD','./LGS2OS','./LGS3OD','./LGS3OS','./LVG1OD','./LVG1OS')

rna.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:8){
  tmp <- readRDS(paste0(path[i],'/meta.rds'))
  meta.list[[i]] <- dplyr::filter(tmp,cell_type__custom != rm_ct[i])
  
  if(i %in% c(3,4,5,6,8)){
    rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds'))[,rownames(meta.list[[i]])]
  }else{
    rna.list[[i]] <- NULL 
  }
  
  if(i %in% c(1,2,4,6,7)){
    peak <- fread(paste0(path[i],'/peak.csv'),data.table = F)
    cells <- fread(paste0(path[i],'/bcd.csv'),data.table = F)
    atac.list[[i]] <- readMM(paste0(path[i],'/atac.mtx'))
    colnames(atac.list[[i]]) <- cells[,1]
    rownames(atac.list[[i]]) <- peak[,1]
    atac.list[[i]] <- atac.list[[i]][,rownames(meta.list[[i]])]
    # atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
  }else{
    atac.list[[i]] <- NULL 
  }
}
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
cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))
mm <- data.frame(batch = c(rep("B1",8422),rep("B2",8734),rep("B3",5681),rep("B4",5654),
                           rep("B5",5649),rep("B6",4990),rep("B7",4804),rep("B8",4359)),
                 modality = c(rep('atac',8422),rep('atac',8734),rep('rna',5681),rep('Multiome',5654),
                              rep('rna',5649),rep('Multiome',4990),rep('atac',4804),rep('rna',4359)),
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
d <- a@Int.result[["bind"]][cellname,]

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
setwd('./Human_Retina')
out_dir <- './Benchmarking/CellType_missing/output/Retina'
meta <- readRDS('./Human_Retina/meta.rds')
set.seed(123)
rm_ct <- sample(unique(meta$celltype),size = 8,replace = T)
rm_ct <- sample(unique(meta$celltype),size = 8,replace = T)
rm_ct <- sample(unique(meta$celltype),size = 8,replace = T)
rm_ct

path <- c('./LGS1OD','./LGS1OS','./LGS2OD','./LGS2OS','./LGS3OD','./LGS3OS','./LVG1OD','./LVG1OS')

rna.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:8){
  tmp <- readRDS(paste0(path[i],'/meta.rds'))
  meta.list[[i]] <- dplyr::filter(tmp,cell_type__custom != rm_ct[i])
  
  if(i %in% c(2,3,4,5,7)){
    rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds'))[,rownames(meta.list[[i]])]
  }else{
    rna.list[[i]] <- NULL 
  }
  
  if(i %in% c(1,2,3,6,8)){
    peak <- fread(paste0(path[i],'/peak.csv'),data.table = F)
    cells <- fread(paste0(path[i],'/bcd.csv'),data.table = F)
    atac.list[[i]] <- readMM(paste0(path[i],'/atac.mtx'))
    colnames(atac.list[[i]]) <- cells[,1]
    rownames(atac.list[[i]]) <- peak[,1]
    atac.list[[i]] <- atac.list[[i]][,rownames(meta.list[[i]])]
    # atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
  }else{
    atac.list[[i]] <- NULL 
  }
}

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
cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))
mm <- data.frame(batch = c(rep("B1",8157),rep("B2",8734),rep("B3",5627),rep("B4",5849),
                           rep("B5",5664),rep("B6",4836),rep("B7",4643),rep("B8",1800)),
                 modality = c(rep('atac',8157),rep('Multiome',8734),rep('Multiome',5627),rep('rna',5849),
                              rep('rna',5664),rep('atac',4836),rep('rna',4643),rep('atac',1800)),
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
d <- a@Int.result[["bind"]][cellname,]

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
setwd('./Human_Retina')
out_dir <- './Benchmarking/CellType_missing/output/Retina'
meta <- readRDS('./Human_Retina/meta.rds')
set.seed(123)
rm_ct <- sample(unique(meta$celltype),size = 8,replace = T)
rm_ct

path <- c('./LGS1OD','./LGS1OS','./LGS2OD','./LGS2OS','./LGS3OD','./LGS3OS','./LVG1OD','./LVG1OS')

rna.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:8){
  tmp <- readRDS(paste0(path[i],'/meta.rds'))
  meta.list[[i]] <- dplyr::filter(tmp,cell_type__custom != rm_ct[i])
  
  if(i %in% c(1,2,3,4,5)){
    rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds'))[,rownames(meta.list[[i]])]
  }else{
    rna.list[[i]] <- NULL 
  }
  
  if(i %in% c(1,2,6,7,8)){
    peak <- fread(paste0(path[i],'/peak.csv'),data.table = F)
    cells <- fread(paste0(path[i],'/bcd.csv'),data.table = F)
    atac.list[[i]] <- readMM(paste0(path[i],'/atac.mtx'))
    colnames(atac.list[[i]]) <- cells[,1]
    rownames(atac.list[[i]]) <- peak[,1]
    atac.list[[i]] <- atac.list[[i]][,rownames(meta.list[[i]])]
    # atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
  }else{
    atac.list[[i]] <- NULL 
  }
  # rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds'))
  # atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
}

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
cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))
mm <- data.frame(batch = c(rep("B1",8361),rep("B2",9110),rep("B3",5681),rep("B4",5593),
                           rep("B5",5678),rep("B6",4701),rep("B7",4997),rep("B8",4594)),
                 modality = c(rep('multiome',8361+9110),rep('rna',5681+5593+5678),rep('atac',4701+4997+4594)),
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
d <- a@Int.result[["bind"]][cellname,]

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
setwd('./Human_Retina')
out_dir <- './Benchmarking/CellType_missing/output/Retina'
meta <- readRDS('./Human_Retina/meta.rds')
set.seed(123)
rm_ct <- sample(unique(meta$celltype),size = 8,replace = T)
rm_ct <- sample(unique(meta$celltype),size = 8,replace = T)
rm_ct <- sample(unique(meta$celltype),size = 8,replace = T)
rm_ct <- sample(unique(meta$celltype),size = 8,replace = T)
rm_ct

path <- c('./LGS1OD','./LGS1OS','./LGS2OD','./LGS2OS','./LGS3OD','./LGS3OS','./LVG1OD','./LVG1OS')

rna.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:8){
  tmp <- readRDS(paste0(path[i],'/meta.rds'))
  meta.list[[i]] <- dplyr::filter(tmp,cell_type__custom != rm_ct[i])
  
  if(i %in% c(2,4,5,6,7)){
    rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds'))[,rownames(meta.list[[i]])]
  }else{
    rna.list[[i]] <- NULL 
  }
  
  if(i %in% c(1,3,4,5,8)){
    peak <- fread(paste0(path[i],'/peak.csv'),data.table = F)
    cells <- fread(paste0(path[i],'/bcd.csv'),data.table = F)
    atac.list[[i]] <- readMM(paste0(path[i],'/atac.mtx'))
    colnames(atac.list[[i]]) <- cells[,1]
    rownames(atac.list[[i]]) <- peak[,1]
    atac.list[[i]] <- atac.list[[i]][,rownames(meta.list[[i]])]
    # atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
  }else{
    atac.list[[i]] <- NULL 
  }
}

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

cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))

mm <- data.frame(batch = c(rep("B1",8009),rep("B2",8734),rep("B3",6193),rep("B4",5715),
                           rep("B5",5607),rep("B6",4701),rep("B7",2014),rep("B8",4695)),
                 modality = c(rep('atac',8009),rep('rna',8734),rep('atac',6193),rep('Multiome',5715),
                              rep('Multiome',5607),rep('rna',4701),rep('rna',2014),rep('atac',4695)),
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
d <- a@Int.result[["bind"]][cellname,]


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
options(future.globals.maxSize = 5000 * 1024^3)
setwd('./Human_Retina')
out_dir <- './Benchmarking/CellType_missing/output/Retina'
meta <- readRDS('./Human_Retina/meta.rds')

set.seed(123)
rm_ct <- sample(unique(meta$celltype),size = 8,replace = T)
rm_ct <- sample(unique(meta$celltype),size = 8,replace = T)
rm_ct <- sample(unique(meta$celltype),size = 8,replace = T)
rm_ct <- sample(unique(meta$celltype),size = 8,replace = T)
rm_ct <- sample(unique(meta$celltype),size = 8,replace = T)
rm_ct

path <- c('./LGS1OD','./LGS1OS','./LGS2OD','./LGS2OS','./LGS3OD','./LGS3OS','./LVG1OD','./LVG1OS')

rna.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:8){
  tmp <- readRDS(paste0(path[i],'/meta.rds'))
  meta.list[[i]] <- dplyr::filter(tmp,cell_type__custom != rm_ct[i])
  
  if(i %in% c(1,3,5,6,7)){
    rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds'))[,rownames(meta.list[[i]])]
  }else{
    rna.list[[i]] <- NULL 
  }
  
  if(i %in% c(2,4,6,7,8)){
    peak <- fread(paste0(path[i],'/peak.csv'),data.table = F)
    cells <- fread(paste0(path[i],'/bcd.csv'),data.table = F)
    atac.list[[i]] <- readMM(paste0(path[i],'/atac.mtx'))
    colnames(atac.list[[i]]) <- cells[,1]
    rownames(atac.list[[i]]) <- peak[,1]
    atac.list[[i]] <- atac.list[[i]][,rownames(meta.list[[i]])]
    # atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
  }else{
    atac.list[[i]] <- NULL 
  }
}

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
cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))
mm <- data.frame(batch = c(rep("B1",8157),rep("B2",8625),rep("B3",6157),rep("B4",2910),
                           rep("B5",5274),rep("B6",4647),rep("B7",4643),rep("B8",1800)),
                 modality = c(rep('rna',8157),rep('atac',8625),rep('rna',6157),rep('atac',2910),
                              rep('rna',5274),rep('Multiome',4647),rep('Multiome',4643),rep('atac',1800)),
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
d <- a@Int.result[["bind"]][cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_s4.csv'),
          row.names = FALSE)
saveRDS(mm,paste0(out_dir,'/meta_s4.rds'))
