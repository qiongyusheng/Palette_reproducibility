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


meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_shuffled/output/Retina/meta_s1.rds')
cellname <- rownames(meta)
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_shuffled5','label_shuffled10','label_shuffled20')
mod <- c('_shuffled5','_shuffled10','_shuffled20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,
                         modal = c('rna','atac'),
                         supervised = T,
                         group = 'celltype',
                         dims = list(1:50,1:50),# 使用哪些维数
                         method = c("PCA","LSI"))
  
  obj <- Find.Subspace3(obj,
                             modal = c("rna",'atac'),
                             lambda = list(0.8,0.5),
                             supervised = T,
                             L2.norm = F,
                             sub.dims = list(1:40,1:30),
                             cos.dims = list(1:50,1:50),
                             Angle.var = c(15,20),
                             max.Angle = c(50,50))
  a <- Run.Palette2(obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  d <- a@Int.result[["bind"]]
  d <- d[cellname,]
  write.csv(as.matrix(d),
            paste0(out_dir,'/Palette_s1',paste0(mod[i]),'.csv'),
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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_shuffled/output/Retina/meta_s2.rds')
cellname <- rownames(meta)
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_shuffled5','label_shuffled10','label_shuffled20')
mod <- c('_shuffled5','_shuffled10','_shuffled20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,
                    modal = c('rna','atac'),
                    supervised = T,
                    group = 'celltype',
                    dims = list(1:50,1:50),# 使用哪些维数
                    method = c("PCA","LSI"))
  
  obj <- Find.Subspace3(obj,
                        modal = c("rna",'atac'),
                        lambda = list(0.8,0.5),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:40,1:30),
                        cos.dims = list(1:50,1:50),
                        Angle.var = c(15,20),
                        max.Angle = c(50,50))
  a <- Run.Palette2(obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  d <- a@Int.result[["bind"]]
  d <- d[cellname,]
  write.csv(as.matrix(d),
            paste0(out_dir,'/Palette_s2',paste0(mod[i]),'.csv'),
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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_shuffled/output/Retina/meta.rds')
cellname <- rownames(meta)
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_shuffled5','label_shuffled10','label_shuffled20')
mod <- c('_shuffled5','_shuffled10','_shuffled20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,
                    modal = c('rna','atac'),
                    supervised = T,
                    group = 'celltype',
                    dims = list(1:50,1:50),# 使用哪些维数
                    method = c("PCA","LSI"))
  
  obj <- Find.Subspace3(obj,
                        modal = c("rna",'atac'),
                        lambda = list(0.8,0.5),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:40,1:30),
                        cos.dims = list(1:50,1:50),
                        Angle.var = c(15,20),
                        max.Angle = c(50,50))
  a <- Run.Palette2(obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  d <- a@Int.result[["bind"]]
  d <- d[cellname,]
  write.csv(as.matrix(d),
            paste0(out_dir,'/Palette',paste0(mod[i]),'.csv'),
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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_shuffled/output/Retina/meta_s3.rds')
cellname <- rownames(meta)
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_shuffled5','label_shuffled10','label_shuffled20')
mod <- c('_shuffled5','_shuffled10','_shuffled20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,
                    modal = c('rna','atac'),
                    supervised = T,
                    group = 'celltype',
                    dims = list(1:50,1:50),# 使用哪些维数
                    method = c("PCA","LSI"))
  
  obj <- Find.Subspace3(obj,
                        modal = c("rna",'atac'),
                        lambda = list(0.8,0.5),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:40,1:30),
                        cos.dims = list(1:50,1:50),
                        Angle.var = c(15,20),
                        max.Angle = c(50,50))
  a <- Run.Palette2(obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  d <- a@Int.result[["bind"]]
  d <- d[cellname,]
  write.csv(as.matrix(d),
            paste0(out_dir,'/Palette_s3',paste0(mod[i]),'.csv'),
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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_shuffled/output/Retina/meta_s4.rds')
cellname <- rownames(meta)
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_shuffled5','label_shuffled10','label_shuffled20')
mod <- c('_shuffled5','_shuffled10','_shuffled20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,
                    modal = c('rna','atac'),
                    supervised = T,
                    group = 'celltype',
                    dims = list(1:50,1:50),# 使用哪些维数
                    method = c("PCA","LSI"))
  
  obj <- Find.Subspace3(obj,
                        modal = c("rna",'atac'),
                        lambda = list(0.8,0.5),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:40,1:30),
                        cos.dims = list(1:50,1:50),
                        Angle.var = c(15,20),
                        max.Angle = c(50,50))
  a <- Run.Palette2(obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  d <- a@Int.result[["bind"]]
  d <- d[cellname,]
  write.csv(as.matrix(d),
            paste0(out_dir,'/Palette_s4',paste0(mod[i]),'.csv'),
            row.names = FALSE)
}
