######################### Random 1 ################################
source("./Benchmarking/Supervised_using_modified_label/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

out_dir <- './Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s2'

setwd('./BMMC_CITE_Multiome/')
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

colnames(atac.list[[2]]) <- paste0(colnames(atac.list[[2]]),'-atac')
colnames(atac.list[[3]]) <- paste0(colnames(atac.list[[3]]),'-atac')

colnames(rna.list[[1]]) <- paste0(colnames(rna.list[[1]]),'-rna')
colnames(rna.list[[2]]) <- paste0(colnames(rna.list[[2]]),'-rna')
colnames(rna.list[[5]]) <- paste0(colnames(rna.list[[5]]),'-rna')
colnames(rna.list[[6]]) <- paste0(colnames(rna.list[[6]]),'-rna')

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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s2/meta_34.rds')
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
            paste0(out_dir,'/Palette_34',paste0(mod[i]),'.csv'),
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

out_dir <- './Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s2'

setwd('./BMMC_CITE_Multiome/')
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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s2/meta_35.rds')
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
            paste0(out_dir,'/Palette_35',paste0(mod[i]),'.csv'),
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

out_dir <- './Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s2'

setwd('./BMMC_CITE_Multiome/')
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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s2/meta_36.rds')
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
            paste0(out_dir,'/Palette_36',paste0(mod[i]),'.csv'),
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

out_dir <- './Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s2'

setwd('./BMMC_CITE_Multiome/')
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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s2/meta_24.rds')
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
            paste0(out_dir,'/Palette_24',paste0(mod[i]),'.csv'),
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

out_dir <- './Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s2'

setwd('./BMMC_CITE_Multiome/')
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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_miss/output/BMMC_s2/meta_26.rds')
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
            paste0(out_dir,'/Palette_26',paste0(mod[i]),'.csv'),
            row.names = FALSE)
}

