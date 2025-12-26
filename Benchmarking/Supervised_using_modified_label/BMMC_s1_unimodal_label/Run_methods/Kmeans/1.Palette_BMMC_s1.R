######################### Random 1 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

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
path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:2){
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))
  meta2.list[[i]] <- readRDS(paste0(path2[i],'/meta.rds'))
}

clusters <- c('Seurat_cluster')
add_name <- c('seurat_')

for(i in 1){
  musa_obj <- Create.Musa.Object(data.list = c(rna1.list,
                                               rna2.list,
                                               adt.list,
                                               atac.list), 
                                 samples = c('CITE1','CITE2','CITE3',
                                             'Multiome1','Multiome2',
                                             'CITE1','CITE2','CITE3',
                                             'Multiome1','Multiome2'), 
                                 modals = c("rna","rna",
                                            "rna","rna","rna",
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

  mm <- readRDS('./Benchmarking/Supervised_using_modified_label/BMMC_s1_unimodal_label/output/Kmeans/meta_drop5.rds')
  mm <- mm[musa_obj@Meta$cell.name,]
  musa_obj@Meta[['celltype']] <- mm[,clusters[i]]
  
  all.equal(musa_obj@Meta$cell.name,mm[,1])
  
  musa_obj <- DownSample(musa_obj,
                         modal = c('rna','adt','atac'),
                         supervised = T,
                         group = 'celltype',
                         dims = list(1:50,1:50,1:50),# 使用哪些维数
                         method = c("PCA","PCA","LSI"))
  options(future.globals.maxSize = 50 * 1024^3)
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
  #a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
  mm <- readRDS('./Benchmarking/Supervised_using_modified_label/BMMC_s1_unimodal_label/output/Kmeans/meta_drop5.rds')
  out_dir <- './Benchmarking/Supervised_using_modified_label/BMMC_s1_unimodal_label/output/Kmeans'
  d <- a@Int.result[["bind"]]
  d <- d[rownames(mm),]
  write.csv(as.matrix(d),
            paste0(out_dir,'/',add_name[i],'Palette_drop5.csv'),
            row.names = FALSE)
  
}

######################### Random 2 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:2){
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

clusters <- c('Seurat_cluster')
add_name <- c('seurat_')

for(i in 1){
  musa_obj <- Create.Musa.Object(data.list = c(rna1.list,
                                               rna2.list,
                                               adt.list,
                                               atac.list), 
                                 samples = c('CITE1','CITE2',
                                             'Multiome1','Multiome2','Multiome3',
                                             'CITE1','CITE2',
                                             'Multiome1','Multiome2','Multiome3'), 
                                 modals = c("rna","rna",
                                            "rna","rna","rna",
                                            "adt",'adt',
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
  
  mm <- readRDS('./Benchmarking/Supervised_using_modified_label/BMMC_s1_unimodal_label/output/Kmeans/meta_drop1.rds')
  mm <- mm[musa_obj@Meta$cell.name,]
  musa_obj@Meta[['celltype']] <- mm[,clusters[i]]
  
  all.equal(musa_obj@Meta$cell.name,mm[,1])
  
  musa_obj <- DownSample(musa_obj,
                         modal = c('rna','adt','atac'),
                         supervised = T,
                         group = 'celltype',
                         dims = list(1:50,1:50,1:50),# 使用哪些维数
                         method = c("PCA","PCA","LSI"))
  options(future.globals.maxSize = 50 * 1024^3)
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
  # a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
  mm <- readRDS('./Benchmarking/Supervised_using_modified_label/BMMC_s1_unimodal_label/output/Kmeans/meta_drop1.rds')
  out_dir <- './Benchmarking/Supervised_using_modified_label/BMMC_s1_unimodal_label/output/Kmeans'
  d <- a@Int.result[["bind"]]
  d <- d[rownames(mm),]
  write.csv(as.matrix(d),
            paste0(out_dir,'/',add_name[i],'Palette_drop1.csv'),
            row.names = FALSE)
}

######################### Random 3 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

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
path2 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:2){
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))
  meta2.list[[i]] <- readRDS(paste0(path2[i],'/meta.rds'))
}
clusters <- c('Seurat_cluster')
add_name <- c('seurat_')

for(i in 1){
  musa_obj <- Create.Musa.Object(data.list = c(rna1.list,
                                               rna2.list,
                                               adt.list,
                                               atac.list), 
                                 samples = c('CITE1','CITE2','CITE3',
                                             'Multiome1','Multiome2',
                                             'CITE1','CITE2','CITE3',
                                             'Multiome1','Multiome2'), 
                                 modals = c("rna","rna",
                                            "rna","rna","rna",
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
  
  mm <- readRDS('./Benchmarking/Supervised_using_modified_label/BMMC_s1_unimodal_label/output/Kmeans/meta_drop6.rds')
  mm <- mm[musa_obj@Meta$cell.name,]
  musa_obj@Meta[['celltype']] <- mm[,clusters[i]]
  
  all.equal(musa_obj@Meta$cell.name,mm[,1])
  
  musa_obj <- DownSample(musa_obj,
                         modal = c('rna','adt','atac'),
                         supervised = T,
                         group = 'celltype',
                         dims = list(1:50,1:50,1:50),# 使用哪些维数
                         method = c("PCA","PCA","LSI"))
  options(future.globals.maxSize = 50 * 1024^3)
  musa_obj <- Find.Subspace3(musa_obj,
                             modal = c("rna",'adt','atac'),
                             lambda = list(0.8,0.8,0.5),
                             supervised = T,
                             L2.norm = F,
                             sub.dims = list(1:40,1:30,1:30),
                             cos.dims = list(1:50,1:50,1:50),
                             Angle.var = c(15,15,20),
                             max.Angle = c(50,50,50))
  a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  #a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
  mm <- readRDS('./Benchmarking/Supervised_using_modified_label/BMMC_s1_unimodal_label/output/Kmeans/meta_drop6.rds')
  out_dir <- './Benchmarking/Supervised_using_modified_label/BMMC_s1_unimodal_label/output/Kmeans'
  d <- a@Int.result[["bind"]]
  d <- d[rownames(mm),]
  write.csv(as.matrix(d),
            paste0(out_dir,'/',add_name[i],'Palette_drop6.csv'),
            row.names = FALSE)
  
}

######################### Random 4 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

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
path2 <- c('./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna2.list <- list()
atac.list <- list()
meta2.list <- list()
for(i in 1:2){
  rna2.list[[i]] <- readRDS(paste0(path2[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path2[i],'/atac.rds'))
  meta2.list[[i]] <- readRDS(paste0(path2[i],'/meta.rds'))
}

clusters <- c('Seurat_cluster')
add_name <- c('seurat_')

for(i in 1){
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
  
  
  mm <- readRDS('./Benchmarking/Supervised_using_modified_label/BMMC_s1_unimodal_label/output/Kmeans/meta_drop4.rds')
  mm <- mm[musa_obj@Meta$cell.name,]
  musa_obj@Meta[['celltype']] <- mm[,clusters[i]]
  all.equal(musa_obj@Meta$cell.name,mm[,1])
  
  musa_obj <- DownSample(musa_obj,
                         modal = c('rna','adt','atac'),
                         supervised = T,
                         group = 'celltype',
                         dims = list(1:50,1:50,1:50),# 使用哪些维数
                         method = c("PCA","PCA","LSI"))
  options(future.globals.maxSize = 50 * 1024^3)
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
  # a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
  mm <- readRDS('./Benchmarking/Supervised_using_modified_label/BMMC_s1_unimodal_label/output/Kmeans/meta_drop4.rds')
  out_dir <- './Benchmarking/Supervised_using_modified_label/BMMC_s1_unimodal_label/output/Kmeans'
  d <- a@Int.result[["bind"]]
  d <- d[rownames(mm),]
  write.csv(as.matrix(d),
            paste0(out_dir,'/',add_name[i],'Palette_drop4.csv'),
            row.names = FALSE)
}

######################### Random 5 ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta1.list <- list()
for(i in 1:2){
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

clusters <- c('Seurat_cluster')
add_name <- c('seurat_')

for(i in 1){
  mm <- readRDS('./Benchmarking/Supervised_using_modified_label/BMMC_s1_unimodal_label/output/Kmeans/meta_drop2.rds')
  mm <- mm[musa_obj@Meta$cell.name,]
  musa_obj@Meta[['celltype']] <- mm[,clusters[i]]
  
  all.equal(musa_obj@Meta$cell.name,mm[,1])
  
  musa_obj <- DownSample(musa_obj,
                         modal = c('rna','adt','atac'),
                         supervised = T,
                         group = 'celltype',
                         dims = list(1:50,1:50,1:50),# 使用哪些维数
                         method = c("PCA","PCA","LSI"))
  options(future.globals.maxSize = 50 * 1024^3)
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
  # a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
  mm <- readRDS('./Benchmarking/Supervised_using_modified_label/BMMC_s1_unimodal_label/output/Kmeans/meta_drop2.rds')
  out_dir <- './Benchmarking/Supervised_using_modified_label/BMMC_s1_unimodal_label/output/Kmeans'
  d <- a@Int.result[["bind"]]
  d <- d[rownames(mm),]
  write.csv(as.matrix(d),
            paste0(out_dir,'/',add_name[i],'Palette_drop2.csv'),
            row.names = FALSE)
}
