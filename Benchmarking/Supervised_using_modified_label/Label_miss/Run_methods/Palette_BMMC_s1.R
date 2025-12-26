######################### Random 1 ################################
source("./Benchmarking/Supervised_using_modified_label/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

out_dir <- './Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s1'

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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s1/meta_drop5.rds')
cellname <- rownames(meta)
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_miss5','label_miss10','label_miss20')
mod <- c('_miss5','_miss10','_miss20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,
                    modal = c('rna','adt','atac'),
                    supervised = T,
                    group = 'celltype',
                    dims = list(1:50,1:50,1:50),# 使用哪些维数
                    method = c("PCA","PCA","LSI"))
  options(future.globals.maxSize = 50 * 1024^3)
  obj <- Find.Subspace3(obj,
                        modal = c("rna",'adt','atac'),
                        lambda = list(0.8,0.8,0.5),
                        alpha = list(0,0,0),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:40,1:30,1:30),
                        cos.dims = list(1:50,1:50,1:50),
                        Angle.var = c(15,15,20),
                        max.Angle = c(50,50,50))
  a <- Run.Palette2(obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  d <- a@Int.result[["bind"]]
  d <- d[cellname,]
  write.csv(as.matrix(d),
            paste0(out_dir,'/Palette_drop5',paste0(mod[i]),'.csv'),
            row.names = FALSE)
}

######################### Random 2 ################################
source("./Benchmarking/Supervised_using_modified_label/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

out_dir <- './Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s1'

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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s1/meta_drop1.rds')
cellname <- rownames(meta)
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_miss5','label_miss10','label_miss20')
mod <- c('_miss5','_miss10','_miss20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,
                    modal = c('rna','adt','atac'),
                    supervised = T,
                    group = 'celltype',
                    dims = list(1:50,1:50,1:50),# 使用哪些维数
                    method = c("PCA","PCA","LSI"))
  options(future.globals.maxSize = 50 * 1024^3)
  obj <- Find.Subspace3(obj,
                        modal = c("rna",'adt','atac'),
                        lambda = list(0.8,0.8,0.5),
                        alpha = list(0,0,0),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:40,1:30,1:30),
                        cos.dims = list(1:50,1:50,1:50),
                        Angle.var = c(15,15,20),
                        max.Angle = c(50,50,50))
  a <- Run.Palette2(obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  d <- a@Int.result[["bind"]]
  d <- d[cellname,]
  write.csv(as.matrix(d),
            paste0(out_dir,'/Palette_drop1',paste0(mod[i]),'.csv'),
            row.names = FALSE)
}

######################### Random 3 ################################
source("./Benchmarking/Supervised_using_modified_label/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

out_dir <- './Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s1'

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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s1/meta_drop6.rds')
cellname <- rownames(meta)
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_miss5','label_miss10','label_miss20')
mod <- c('_miss5','_miss10','_miss20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,
                    modal = c('rna','adt','atac'),
                    supervised = T,
                    group = 'celltype',
                    dims = list(1:50,1:50,1:50),# 使用哪些维数
                    method = c("PCA","PCA","LSI"))
  options(future.globals.maxSize = 50 * 1024^3)
  obj <- Find.Subspace3(obj,
                        modal = c("rna",'adt','atac'),
                        lambda = list(0.8,0.8,0.5),
                        alpha = list(0,0,0),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:40,1:30,1:30),
                        cos.dims = list(1:50,1:50,1:50),
                        Angle.var = c(15,15,20),
                        max.Angle = c(50,50,50))
  a <- Run.Palette2(obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  d <- a@Int.result[["bind"]]
  d <- d[cellname,]
  write.csv(as.matrix(d),
            paste0(out_dir,'/Palette_drop6',paste0(mod[i]),'.csv'),
            row.names = FALSE)
}

######################### Random 4 ################################
source("./Benchmarking/Supervised_using_modified_label/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

out_dir <- './Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s1'

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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s1/meta_drop4.rds')
cellname <- rownames(meta)
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_miss5','label_miss10','label_miss20')
mod <- c('_miss5','_miss10','_miss20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,
                    modal = c('rna','adt','atac'),
                    supervised = T,
                    group = 'celltype',
                    dims = list(1:50,1:50,1:50),# 使用哪些维数
                    method = c("PCA","PCA","LSI"))
  options(future.globals.maxSize = 50 * 1024^3)
  obj <- Find.Subspace3(obj,
                        modal = c("rna",'adt','atac'),
                        lambda = list(0.8,0.8,0.5),
                        alpha = list(0,0,0),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:40,1:30,1:30),
                        cos.dims = list(1:50,1:50,1:50),
                        Angle.var = c(15,15,20),
                        max.Angle = c(50,50,50))
  a <- Run.Palette2(obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  d <- a@Int.result[["bind"]]
  d <- d[cellname,]
  write.csv(as.matrix(d),
            paste0(out_dir,'/Palette_drop4',paste0(mod[i]),'.csv'),
            row.names = FALSE)
}

######################### Random 5 ################################
source("./Benchmarking/Supervised_using_modified_label/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

out_dir <- './Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s1'
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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s1/meta_drop2.rds')
cellname <- rownames(meta)
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_miss5','label_miss10','label_miss20')
mod <- c('_miss5','_miss10','_miss20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,
                    modal = c('rna','adt','atac'),
                    supervised = T,
                    group = 'celltype',
                    dims = list(1:50,1:50,1:50),# 使用哪些维数
                    method = c("PCA","PCA","LSI"))
  options(future.globals.maxSize = 50 * 1024^3)
  obj <- Find.Subspace3(obj,
                        modal = c("rna",'adt','atac'),
                        lambda = list(0.8,0.8,0.5),
                        alpha = list(0,0,0),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:40,1:30,1:30),
                        cos.dims = list(1:50,1:50,1:50),
                        Angle.var = c(15,15,20),
                        max.Angle = c(50,50,50))
  a <- Run.Palette2(obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  d <- a@Int.result[["bind"]]
  d <- d[cellname,]
  write.csv(as.matrix(d),
            paste0(out_dir,'/Palette_drop2',paste0(mod[i]),'.csv'),
            row.names = FALSE)
}





