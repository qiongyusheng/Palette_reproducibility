######################### Random 1 ################################
source("./Benchmarking/Supervised_using_modified_label/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

out_dir <- './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/TEA_s1'
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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_shuffled/output/TEA_s1/meta.rds')
cellname <- rownames(meta)
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_shuffled5','label_shuffled10','label_shuffled20')
mod <- c('_shuffled5','_shuffled10','_shuffled20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,
                    modal = c('rna','adt','atac'),
                    supervised = T,
                    group = 'celltype',
                    dims = list(1:50,1:20,1:50),# 使用哪些维数
                    method = c("PCA","PCA","LSI"))
  options(future.globals.maxSize = 500 * 1024^3)
  obj <- Find.Subspace3(obj,
                        modal = c("rna",'adt','atac'),
                        lambda = list(0.8,0.8,0.5),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:40,1:20,1:30),
                        cos.dims = list(1:50,1:20,1:50),
                        Angle.var = c(15,15,20),
                        max.Angle = c(50,50,50))
  a <- Run.Palette2(obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  
  d <- a@Int.result[["bind"]]
  d <- d[cellname,]
  write.csv(as.matrix(d),
            paste0(out_dir,'/Palette_v1',paste0(mod[i]),'.csv'),
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

out_dir <- './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/TEA_s1'
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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_shuffled/output/TEA_s1/meta.rds')
cellname <- rownames(meta)
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_shuffled5','label_shuffled10','label_shuffled20')
mod <- c('_shuffled5','_shuffled10','_shuffled20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,
                    modal = c('rna','adt','atac'),
                    supervised = T,
                    group = 'celltype',
                    dims = list(1:50,1:20,1:50),# 使用哪些维数
                    method = c("PCA","PCA","LSI"))
  options(future.globals.maxSize = 500 * 1024^3)
  obj <- Find.Subspace3(obj,
                        modal = c("rna",'adt','atac'),
                        lambda = list(0.8,0.8,0.5),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:40,1:20,1:30),
                        cos.dims = list(1:50,1:20,1:50),
                        Angle.var = c(15,15,20),
                        max.Angle = c(50,50,50))
  a <- Run.Palette2(obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  
  d <- a@Int.result[["bind"]]
  d <- d[cellname,]
  write.csv(as.matrix(d),
            paste0(out_dir,'/Palette',paste0(mod[i]),'.csv'),
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

out_dir <- './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/TEA_s1'
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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_shuffled/output/TEA_s1/meta.rds')
cellname <- rownames(meta)
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_shuffled5','label_shuffled10','label_shuffled20')
mod <- c('_shuffled5','_shuffled10','_shuffled20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,
                    modal = c('rna','adt','atac'),
                    supervised = T,
                    group = 'celltype',
                    dims = list(1:50,1:20,1:50),# 使用哪些维数
                    method = c("PCA","PCA","LSI"))
  options(future.globals.maxSize = 500 * 1024^3)
  obj <- Find.Subspace3(obj,
                        modal = c("rna",'adt','atac'),
                        lambda = list(0.8,0.8,0.5),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:40,1:20,1:30),
                        cos.dims = list(1:50,1:20,1:50),
                        Angle.var = c(15,15,20),
                        max.Angle = c(50,50,50))
  a <- Run.Palette2(obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  
  d <- a@Int.result[["bind"]]
  d <- d[cellname,]
  write.csv(as.matrix(d),
            paste0(out_dir,'/Palette_v2',paste0(mod[i]),'.csv'),
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

out_dir <- './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/TEA_s1'
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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_shuffled/output/TEA_s1/meta.rds')
cellname <- rownames(meta)
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_shuffled5','label_shuffled10','label_shuffled20')
mod <- c('_shuffled5','_shuffled10','_shuffled20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,
                    modal = c('rna','adt','atac'),
                    supervised = T,
                    group = 'celltype',
                    dims = list(1:50,1:20,1:50),# 使用哪些维数
                    method = c("PCA","PCA","LSI"))
  options(future.globals.maxSize = 500 * 1024^3)
  obj <- Find.Subspace3(obj,
                        modal = c("rna",'adt','atac'),
                        lambda = list(0.8,0.8,0.5),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:40,1:20,1:30),
                        cos.dims = list(1:50,1:20,1:50),
                        Angle.var = c(15,15,20),
                        max.Angle = c(50,50,50))
  a <- Run.Palette2(obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  
  d <- a@Int.result[["bind"]]
  d <- d[cellname,]
  write.csv(as.matrix(d),
            paste0(out_dir,'/Palette_v3',paste0(mod[i]),'.csv'),
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

out_dir <- './Benchmarking/Supervised_using_modified_label/Label_shuffled/output/TEA_s1'
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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_shuffled/output/TEA_s1/meta.rds')
cellname <- rownames(meta)
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_shuffled5','label_shuffled10','label_shuffled20')
mod <- c('_shuffled5','_shuffled10','_shuffled20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,
                    modal = c('rna','adt','atac'),
                    supervised = T,
                    group = 'celltype',
                    dims = list(1:50,1:20,1:50),# 使用哪些维数
                    method = c("PCA","PCA","LSI"))
  options(future.globals.maxSize = 500 * 1024^3)
  obj <- Find.Subspace3(obj,
                        modal = c("rna",'adt','atac'),
                        lambda = list(0.8,0.8,0.5),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:40,1:20,1:30),
                        cos.dims = list(1:50,1:20,1:50),
                        Angle.var = c(15,15,20),
                        max.Angle = c(50,50,50))
  a <- Run.Palette2(obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  
  d <- a@Int.result[["bind"]]
  d <- d[cellname,]
  write.csv(as.matrix(d),
            paste0(out_dir,'/Palette_v4',paste0(mod[i]),'.csv'),
            row.names = FALSE)
}
