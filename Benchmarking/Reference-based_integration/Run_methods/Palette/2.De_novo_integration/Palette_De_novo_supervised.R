######################### ADT ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:3){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta1.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}
path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:3){
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))
  meta2.list[[i]] <- readRDS(paste0(path2[i],'/meta.rds'))
}
options(future.globals.maxSize = 500 * 1024^3)
path1 <- c('./CITE/s2d1_cite','./CITE/s2d4_cite','./CITE/s2d5_cite')
q <- list()
meta_q.list <- list()
for(i in 1:3){
  q[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
  # meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$modality,'-',meta_q.list[[i]]$batch)
  meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$batch,'-',meta_q.list[[i]]$modality)
  meta_q.list[[i]]$modality <- 'ADT'
}
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/adt/meta_supervised.rds')
musa_obj <- Create.Musa.Object(data.list = c(rna1.list,
                                             rna2.list,
                                             adt.list,
                                             atac.list,q), 
                               samples = c('CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'ADT1','ADT2','ADT3'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna","rna",
                                          "adt",'adt','adt',
                                          "atac","atac","atac",
                                          'adt','adt','adt'),
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

meta <- meta[musa_obj@Meta$cell.name,]
musa_obj@Meta[['celltype']] <- meta$pred_l2_fastMNN
all.equal(musa_obj@Meta$cell.name,rownames(meta))

musa_obj <- DownSample(musa_obj,
                       modal = c('rna','adt','atac'),
                       supervised = T,
                       group = 'celltype',
                       dims = list(1:50,1:50,1:50),# 使用哪些维数
                       method = c("PCA","PCA","LSI"))

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           alpha = list(0,0,0),
                           supervised = T,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/adt/meta.rds')
cellname <- rownames(meta)
out_dir <- './Benchmarking/Reference-based_integration/output/Result/adt/'
d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'Palette_denovo_supervised.csv'),
          row.names = FALSE)


######################### ADT + ATAC ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:3){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta1.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}
path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:3){
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))
  meta2.list[[i]] <- readRDS(paste0(path2[i],'/meta.rds'))
}
ref <- c(rna1.list,rna2.list,adt.list,atac.list)
rm(rna1.list,rna2.list,adt.list,atac.list)


options(future.globals.maxSize = 500 * 1024^3)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s2d1_cite','./CITE/s2d4_cite','./CITE/s2d5_cite',
           './Multiome/s2d1_multi','./Multiome/s2d4_multi','./Multiome/s2d5_multi')

# rna.list <- list()
adt.list <- list()
atac.list <- list()
meta_q.list <- list()
for(i in 1:6){
  if(i %in% c(1,2,3)){
    adt.list[[i]] <-  readRDS(paste0(path1[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[i-3]] <-  readRDS(paste0(path1[i],'/atac.rds'))
  }
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
  # meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$modality,'-',meta_q.list[[i]]$batch)
  meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$batch,'-',meta_q.list[[i]]$modality)
  
  if(i %in% c(1,2,3)){
    meta_q.list[[i]]$modality <- 'ADT'
  }else{
    meta_q.list[[i]]$modality <- 'ATAC'
  }
  
}
q <- c(adt.list,atac.list)
rm(adt.list,atac.list)

musa_obj <- Create.Musa.Object(data.list = c(ref,q), 
                               samples = c('CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'ADT1','ADT2','ADT3',
                                           'ATAC1','ATAC2','ATAC3'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna","rna",
                                          "adt",'adt','adt',
                                          "atac","atac","atac",
                                          "adt",'adt','adt',
                                          "atac","atac","atac"),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna","atac",'adt'), 
                               sparce = TRUE)

rm(ref,q)

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

meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/adt_atac/meta_supervised.rds')
meta <- meta[musa_obj@Meta$cell.name,]
musa_obj@Meta[['celltype']] <- meta$pred_l2_fastMNN
all.equal(musa_obj@Meta$cell.name,rownames(meta))

musa_obj <- DownSample(musa_obj,
                       modal = c('rna','adt','atac'),
                       supervised = T,
                       group = 'celltype',
                       dims = list(1:50,1:50,1:50),# 使用哪些维数
                       method = c("PCA","PCA","LSI"))

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           alpha = list(0,0,0),
                           supervised = T,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/adt_atac/meta.rds')
cellname <- rownames(meta)
out_dir <- './Benchmarking/Reference-based_integration/output/Result/adt_atac/'
d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'Palette_denovo_supervised.csv'),
          row.names = FALSE)

######################### ATAC ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:3){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta1.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}
path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:3){
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))
  meta2.list[[i]] <- readRDS(paste0(path2[i],'/meta.rds'))
}

setwd('./BMMC_CITE_Multiome')
path1 <- c('./Multiome/s2d1_multi','./Multiome/s2d4_multi','./Multiome/s2d5_multi')
q <- list()
meta_q.list <- list()
for(i in 1:3){
  q[[i]] <- readRDS(paste0(path1[i],'/atac.rds'))
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
  # meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$modality,'-',meta_q.list[[i]]$batch)
  meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$batch,'-',meta_q.list[[i]]$modality)
  meta_q.list[[i]]$modality <- 'ATAC'
}
options(future.globals.maxSize = 500 * 1024^3)
musa_obj <- Create.Musa.Object(data.list = c(rna1.list,
                                             rna2.list,
                                             adt.list,
                                             atac.list,q), 
                               samples = c('CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'ATAC1','ATAC2','ATAC3'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna","rna",
                                          "adt",'adt','adt',
                                          "atac","atac","atac",
                                          'atac','atac','atac'),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna","atac",'adt'), 
                               sparce = TRUE)

rm(rna1.list,
   rna2.list,
   adt.list,
   atac.list,q)
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

meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/atac/meta_supervised.rds')
meta <- meta[musa_obj@Meta$cell.name,]
musa_obj@Meta[['celltype']] <- meta$pred_l2_fastMNN
all.equal(musa_obj@Meta$cell.name,rownames(meta))

musa_obj <- DownSample(musa_obj,
                       modal = c('rna','adt','atac'),
                       supervised = T,
                       group = 'celltype',
                       dims = list(1:50,1:50,1:50),# 使用哪些维数
                       method = c("PCA","PCA","LSI"))

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           alpha = list(0,0,0),
                           supervised = T,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')

meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/atac/meta.rds')
cellname <- rownames(meta)
out_dir <- './Benchmarking/Reference-based_integration/output/Result/atac/'
d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'Palette_denovo_supervised.csv'),
          row.names = FALSE)

######################### CITE ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:3){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta1.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}
path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:3){
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))
  meta2.list[[i]] <- readRDS(paste0(path2[i],'/meta.rds'))
}
ref <- c(rna1.list,rna2.list,adt.list,atac.list)
rm(rna1.list,rna2.list,adt.list,atac.list)
options(future.globals.maxSize = 500 * 1024^3)

setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s2d1_cite','./CITE/s2d4_cite','./CITE/s2d5_cite')
rna1.list <- list()
adt.list <- list()
meta_q.list <- list()
for(i in 1:3){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}
q <- c(rna1.list,adt.list)
rm(rna1.list,adt.list)
gc()

musa_obj <- Create.Musa.Object(data.list = c(ref,q), 
                               samples = c('CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'CITE4','CITE5','CITE6',
                                           'CITE4','CITE5','CITE6'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna","rna",
                                          "adt",'adt','adt',
                                          "atac","atac","atac",
                                          "rna","rna","rna",
                                          "adt",'adt','adt'),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna","atac",'adt'), 
                               sparce = TRUE)

rm(ref,q)

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

meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/cite/meta_supervised.rds')
meta <- meta[musa_obj@Meta$cell.name,]
musa_obj@Meta[['celltype']] <- meta$pred_l2_fastMNN
all.equal(musa_obj@Meta$cell.name,rownames(meta))

musa_obj <- DownSample(musa_obj,
                       modal = c('rna','adt','atac'),
                       supervised = T,
                       group = 'celltype',
                       dims = list(1:50,1:50,1:50),# 使用哪些维数
                       method = c("PCA","PCA","LSI"))

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           alpha = list(0,0,0),
                           supervised = T,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/cite/meta.rds')
cellname <- rownames(meta)
out_dir <- './Benchmarking/Reference-based_integration/output/Result/cite/'
d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'Palette_denovo_supervised.csv'),
          row.names = FALSE)


######################### CITE + Multiome ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:3){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta1.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}
path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:3){
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))
  meta2.list[[i]] <- readRDS(paste0(path2[i],'/meta.rds'))
}
ref <- c(rna1.list,rna2.list,adt.list,atac.list)
rm(rna1.list,rna2.list,adt.list,atac.list)

options(future.globals.maxSize = 500 * 1024^3)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s2d1_cite','./CITE/s2d4_cite','./CITE/s2d5_cite',
           './Multiome/s2d1_multi','./Multiome/s2d4_multi','./Multiome/s2d5_multi')

rna.list <- list()
adt.list <- list()
atac.list <- list()
meta_q.list <- list()
for(i in 1:6){
  rna.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  if(i %in% c(1,2,3)){
    adt.list[[i]] <-  readRDS(paste0(path1[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[i-3]] <-  readRDS(paste0(path1[i],'/atac.rds'))
  }
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
  # meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$modality,'-',meta_q.list[[i]]$batch)
  meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$batch,'-',meta_q.list[[i]]$modality)
}
q <- c(rna.list,adt.list,atac.list)
rm(rna.list,adt.list,atac.list)
gc()

musa_obj <- Create.Musa.Object(data.list = c(ref,q), 
                               samples = c('CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'CITE4','CITE5','CITE6',
                                           'Multiome4','Multiome5','Multiome6',
                                           'CITE4','CITE5','CITE6',
                                           'Multiome4','Multiome5','Multiome6'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna","rna",
                                          "adt",'adt','adt',
                                          "atac","atac","atac",
                                          "rna","rna","rna",
                                          "rna","rna","rna",
                                          "adt",'adt','adt',
                                          "atac","atac","atac"),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna","atac",'adt'), 
                               sparce = TRUE)

rm(ref,q)

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

meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/cite_multi/meta_supervised.rds')
meta <- meta[musa_obj@Meta$cell.name,]
musa_obj@Meta[['celltype']] <- meta$pred_l2_fastMNN
all.equal(musa_obj@Meta$cell.name,rownames(meta))

musa_obj <- DownSample(musa_obj,
                       modal = c('rna','adt','atac'),
                       supervised = T,
                       group = 'celltype',
                       dims = list(1:50,1:50,1:50),# 使用哪些维数
                       method = c("PCA","PCA","LSI"))

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           alpha = list(0,0,0),
                           supervised = T,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/cite_multi/meta.rds')
cellname <- rownames(meta)
out_dir <- './Benchmarking/Reference-based_integration/output/Result/cite_multi/'
d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'Palette_denovo_supervised.csv'),
          row.names = FALSE)

######################### Multiome ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:3){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta1.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}
path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:3){
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))
  meta2.list[[i]] <- readRDS(paste0(path2[i],'/meta.rds'))
}
ref <- c(rna1.list,rna2.list,adt.list,atac.list)
rm(rna1.list,rna2.list,adt.list,atac.list)



setwd('./BMMC_CITE_Multiome')
path1 <- c('./Multiome/s2d1_multi','./Multiome/s2d4_multi','./Multiome/s2d5_multi')
rna1.list <- list()
atac.list <- list()
meta_q.list <- list()
for(i in 1:3){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path1[i],'/atac.rds'))
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}
q <- c(rna1.list,atac.list)
rm(rna1.list,atac.list)
gc()

musa_obj <- Create.Musa.Object(data.list = c(ref,q), 
                               samples = c('CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'Multiome4','Multiome5','Multiome6',
                                           'Multiome4','Multiome5','Multiome6'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna","rna",
                                          "adt",'adt','adt',
                                          "atac","atac","atac",
                                          "rna","rna","rna",
                                          "atac","atac","atac"),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna","atac",'adt'), 
                               sparce = TRUE)

rm(ref,q)

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

meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/multiome/meta_supervised.rds')
meta <- meta[musa_obj@Meta$cell.name,]
musa_obj@Meta[['celltype']] <- meta$pred_l2_fastMNN
all.equal(musa_obj@Meta$cell.name,rownames(meta))
options(future.globals.maxSize = 500 * 1024^3)
musa_obj <- DownSample(musa_obj,
                       modal = c('rna','adt','atac'),
                       supervised = T,
                       group = 'celltype',
                       dims = list(1:50,1:50,1:50),# 使用哪些维数
                       method = c("PCA","PCA","LSI"))

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           alpha = list(0,0,0),
                           supervised = T,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/multiome/meta.rds')
cellname <- rownames(meta)
out_dir <- './Benchmarking/Reference-based_integration/output/Result/multiome/'
d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'Palette_denovo_supervised.csv'),
          row.names = FALSE)


######################### RNA ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:3){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta1.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}
path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:3){
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))
  meta2.list[[i]] <- readRDS(paste0(path2[i],'/meta.rds'))
}
options(future.globals.maxSize = 500 * 1024^3)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s2d1_cite','./CITE/s2d4_cite','./CITE/s2d5_cite',
           './Multiome/s2d1_multi','./Multiome/s2d4_multi','./Multiome/s2d5_multi')
q <- list()
meta_q.list <- list()
for(i in 1:6){
  q[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
  # meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$modality,'-',meta_q.list[[i]]$batch)
  meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$batch,'-',meta_q.list[[i]]$modality)
  meta_q.list[[i]]$modality <- 'RNA'
}

musa_obj <- Create.Musa.Object(data.list = c(rna1.list,
                                             rna2.list,
                                             adt.list,
                                             atac.list,q), 
                               samples = c('CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'RNA1','RNA2','RNA3',
                                           'RNA4','RNA5','RNA6'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna","rna",
                                          "adt",'adt','adt',
                                          "atac","atac","atac",
                                          "rna","rna","rna",
                                          "rna","rna","rna"),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna","atac",'adt'), 
                               sparce = TRUE)

rm(rna1.list,
   rna2.list,
   adt.list,
   atac.list,q)
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

meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/rna/meta_supervised.rds')
meta <- meta[musa_obj@Meta$cell.name,]
musa_obj@Meta[['celltype']] <- meta$pred_l2_fastMNN
all.equal(musa_obj@Meta$cell.name,rownames(meta))

musa_obj <- DownSample(musa_obj,
                       modal = c('rna','adt','atac'),
                       supervised = T,
                       group = 'celltype',
                       dims = list(1:50,1:50,1:50),# 使用哪些维数
                       method = c("PCA","PCA","LSI"))

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           alpha = list(0,0,0),
                           supervised = T,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')

meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/rna/meta.rds')
cellname <- rownames(meta)
out_dir <- './Benchmarking/Reference-based_integration/output/Result/rna/'
d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'Palette_denovo_supervised.csv'),
          row.names = FALSE)


######################### RNA + ADT ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:3){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta1.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}
path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:3){
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))
  meta2.list[[i]] <- readRDS(paste0(path2[i],'/meta.rds'))
}
ref <- c(rna1.list,rna2.list,adt.list,atac.list)
rm(rna1.list,rna2.list,adt.list,atac.list)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s2d1_cite','./CITE/s2d4_cite','./CITE/s2d5_cite',
           './Multiome/s2d1_multi','./Multiome/s2d4_multi','./Multiome/s2d5_multi')

rna.list <- list()
adt.list <- list()
# atac.list <- list()
meta_q.list <- list()
for(i in 1:6){
  if(i %in% c(1,2,3)){
    adt.list[[i]] <-  readRDS(paste0(path1[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    rna.list[[i-3]] <-  readRDS(paste0(path1[i],'/rna.rds'))
  }
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
  # meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$modality,'-',meta_q.list[[i]]$batch)
  meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$batch,'-',meta_q.list[[i]]$modality)
  if(i %in% c(1,2,3)){
    meta_q.list[[i]]$modality <- 'ADT'
  }else{
    meta_q.list[[i]]$modality <- 'RNA'
  }
}
q <- c(adt.list,rna.list)
rm(rna.list,adt.list)
options(future.globals.maxSize = 500 * 1024^3)
musa_obj <- Create.Musa.Object(data.list = c(ref,q), 
                               samples = c('CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'ADT1','ADT2','ADT3',
                                           'RNA1','RNA2','RNA3'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna","rna",
                                          "adt",'adt','adt',
                                          "atac","atac","atac",
                                          "adt",'adt','adt',
                                          "rna","rna","rna"),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna","atac",'adt'), 
                               sparce = TRUE)
rm(ref,q)
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
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/rna_adt/meta_supervised.rds')

meta <- meta[musa_obj@Meta$cell.name,]
musa_obj@Meta[['celltype']] <- meta$pred_l2_fastMNN
all.equal(musa_obj@Meta$cell.name,rownames(meta))

musa_obj <- DownSample(musa_obj,
                       modal = c('rna','adt','atac'),
                       supervised = T,
                       group = 'celltype',
                       dims = list(1:50,1:50,1:50),# 使用哪些维数
                       method = c("PCA","PCA","LSI"))

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           alpha = list(0,0,0),
                           supervised = T,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')

meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/rna_adt/meta.rds')
cellname <- rownames(meta)
out_dir <- './Benchmarking/Reference-based_integration/output/Result/rna_adt/'
d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'Palette_denovo_supervised.csv'),
          row.names = FALSE)

######################### RNA + ADT + ATAC ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:3){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta1.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}
path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:3){
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))
  meta2.list[[i]] <- readRDS(paste0(path2[i],'/meta.rds'))
}
ref <- c(rna1.list,rna2.list,adt.list,atac.list)
rm(rna1.list,rna2.list,adt.list,atac.list)


options(future.globals.maxSize = 500 * 1024^3)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s2d1_cite','./CITE/s2d4_cite','./CITE/s2d5_cite',
           './Multiome/s2d1_multi','./Multiome/s2d4_multi','./Multiome/s2d5_multi')

rna.list <- list()
adt.list <- list()
atac.list <- list()
meta_q.list <- list()
for(i in 1:6){
  if(i %in% c(1,2)){
    adt.list[[i]] <-  readRDS(paste0(path1[i],'/adt.rds'))
  }
  if(i %in% c(3,4)){
    rna.list[[i-2]] <-  readRDS(paste0(path1[i],'/rna.rds'))
  }
  if(i %in% c(5,6)){
    atac.list[[i-4]] <-  readRDS(paste0(path1[i],'/atac.rds'))
  }
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
  # meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$modality,'-',meta_q.list[[i]]$batch)
  meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$batch,'-',meta_q.list[[i]]$modality)
  if(i %in% c(1,2)){
    meta_q.list[[i]]$modality <- 'ADT'
  }
  if(i %in% c(3,4)){
    meta_q.list[[i]]$modality <- 'RNA'
  }
  if(i %in% c(5,6)){
    meta_q.list[[i]]$modality <- 'ATAC'
  }
  
}
q <- c(adt.list,rna.list,atac.list)
rm(adt.list,atac.list,rna.list)

musa_obj <- Create.Musa.Object(data.list = c(ref,q), 
                               samples = c('CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'ADT1','ADT2',
                                           'RNA1','RNA2',
                                           'ATAC1','ATAC2'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna","rna",
                                          "adt",'adt','adt',
                                          "atac","atac","atac",
                                          "adt",'adt',
                                          "rna","rna",
                                          "atac","atac"),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna","atac",'adt'), 
                               sparce = TRUE)

rm(ref,q)

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

meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/rna_adt_atac/meta_supervised.rds')
meta <- meta[musa_obj@Meta$cell.name,]
musa_obj@Meta[['celltype']] <- meta$pred_l2_fastMNN
all.equal(musa_obj@Meta$cell.name,rownames(meta))

musa_obj <- DownSample(musa_obj,
                       modal = c('rna','adt','atac'),
                       supervised = T,
                       group = 'celltype',
                       dims = list(1:50,1:50,1:50),# 使用哪些维数
                       method = c("PCA","PCA","LSI"))

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           alpha = list(0,0,0),
                           supervised = T,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/rna_adt_atac/meta.rds')
cellname <- rownames(meta)
out_dir <- './Benchmarking/Reference-based_integration/output/Result/rna_adt_atac/'
d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'Palette_denovo_supervised.csv'),
          row.names = FALSE)

######################### RNA + ATAC ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:3){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta1.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}
path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:3){
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))
  meta2.list[[i]] <- readRDS(paste0(path2[i],'/meta.rds'))
}
ref <- c(rna1.list,rna2.list,adt.list,atac.list)
rm(rna1.list,rna2.list,adt.list,atac.list)

options(future.globals.maxSize = 500 * 1024^3)

setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s2d1_cite','./CITE/s2d4_cite','./CITE/s2d5_cite',
           './Multiome/s2d1_multi','./Multiome/s2d4_multi','./Multiome/s2d5_multi')

rna.list <- list()
# adt.list <- list()
atac.list <- list()
meta_q.list <- list()
for(i in 1:6){
  if(i %in% c(1,2,3)){
    rna.list[[i]] <-  readRDS(paste0(path1[i],'/rna.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[i-3]] <-  readRDS(paste0(path1[i],'/atac.rds'))
  }
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
  # meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$modality,'-',meta_q.list[[i]]$batch)
  meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$batch,'-',meta_q.list[[i]]$modality)
  
  if(i %in% c(1,2,3)){
    meta_q.list[[i]]$modality <- 'RNA'
  }
  if(i %in% c(4,5,6)){
    meta_q.list[[i]]$modality <- 'ATAC'
  }
}
q <- c(rna.list,atac.list)
rm(rna.list,atac.list)

musa_obj <- Create.Musa.Object(data.list = c(ref,q), 
                               samples = c('CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'RNA1','RNA2','RNA3',
                                           'ATAC1','ATAC2','ATAC3'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna","rna",
                                          "adt",'adt','adt',
                                          "atac","atac","atac",
                                          "rna","rna","rna",
                                          "atac","atac","atac"),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna","atac",'adt'), 
                               sparce = TRUE)

rm(ref,q)

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

meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/rna_atac/meta_supervised.rds')
meta <- meta[musa_obj@Meta$cell.name,]
musa_obj@Meta[['celltype']] <- meta$pred_l2_fastMNN
all.equal(musa_obj@Meta$cell.name,rownames(meta))

musa_obj <- DownSample(musa_obj,
                       modal = c('rna','adt','atac'),
                       supervised = T,
                       group = 'celltype',
                       dims = list(1:50,1:50,1:50),# 使用哪些维数
                       method = c("PCA","PCA","LSI"))

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           alpha = list(0,0,0),
                           supervised = T,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Result/rna_atac/meta.rds')
cellname <- rownames(meta)
out_dir <- './Benchmarking/Reference-based_integration/output/Result/rna_atac/'
d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'Palette_denovo_supervised.csv'),
          row.names = FALSE)

