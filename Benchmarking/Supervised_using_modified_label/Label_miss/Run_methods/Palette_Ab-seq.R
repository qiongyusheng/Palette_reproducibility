######################### Random 1 ################################
source("./Benchmarking/Supervised_using_modified_label/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
setwd('./Young_Age')
out_dir <- './Benchmarking/Supervised_using_modified_label/Label_miss/output/Ab-seq'

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in paths){
  adt.list[[i]] <- readRDS(paste0('./',i,'/adt.rds'))
  rna.list[[i]] <- readRDS(paste0('./',i,'/rna.rds'))
  meta.list[[i]] <- readRDS(paste0('./',i,'/meta.rds'))
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
meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_miss/output/Ab-seq/meta_s1.rds')

cellname <- meta[,2]
rownames(meta) <- meta[,2]
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_miss5','label_miss10','label_miss20')
mod <- c('_miss5','_miss10','_miss20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,supervised = T,
                         group = 'celltype',
                         dims = list(1:50,1:50),# 使用哪些维数
                         method = c("PCA","PCA"))
  
  obj <- Find.Subspace3(obj,
                             modal = c("rna",'adt'),
                             lambda = list(0.8,0.8),
                             alpha = list(0,0),
                             supervised = T,
                             L2.norm = F,
                             sub.dims = list(1:40,1:30),
                             cos.dims = list(1:50,1:50),
                             Angle.var = c(15,15),
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
setwd('./Young_Age')
out_dir <- './Benchmarking/Supervised_using_modified_label/Label_miss/output/Ab-seq'

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in paths){
  adt.list[[i]] <- readRDS(paste0('./',i,'/adt.rds'))
  rna.list[[i]] <- readRDS(paste0('./',i,'/rna.rds'))
  meta.list[[i]] <- readRDS(paste0('./',i,'/meta.rds'))
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
meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_miss/output/Ab-seq/meta_s2.rds')

cellname <- meta[,2]
rownames(meta) <- meta[,2]
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_miss5','label_miss10','label_miss20')
mod <- c('_miss5','_miss10','_miss20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,supervised = T,
                    group = 'celltype',
                    dims = list(1:50,1:50),# 使用哪些维数
                    method = c("PCA","PCA"))
  
  obj <- Find.Subspace3(obj,
                        modal = c("rna",'adt'),
                        lambda = list(0.8,0.8),
                        alpha = list(0,0),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:40,1:30),
                        cos.dims = list(1:50,1:50),
                        Angle.var = c(15,15),
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
setwd('./Young_Age')
out_dir <- './Benchmarking/Supervised_using_modified_label/Label_miss/output/Ab-seq'
paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in paths){
  adt.list[[i]] <- readRDS(paste0('./',i,'/adt.rds'))
  rna.list[[i]] <- readRDS(paste0('./',i,'/rna.rds'))
  meta.list[[i]] <- readRDS(paste0('./',i,'/meta.rds'))
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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_miss/output/Ab-seq/meta_s3.rds')

cellname <- meta[,2]
rownames(meta) <- meta[,2]
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_miss5','label_miss10','label_miss20')
mod <- c('_miss5','_miss10','_miss20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,supervised = T,
                    group = 'celltype',
                    dims = list(1:50,1:50),# 使用哪些维数
                    method = c("PCA","PCA"))
  
  obj <- Find.Subspace3(obj,
                        modal = c("rna",'adt'),
                        lambda = list(0.8,0.8),
                        alpha = list(0,0),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:40,1:30),
                        cos.dims = list(1:50,1:50),
                        Angle.var = c(15,15),
                        max.Angle = c(50,50))
  
  a <- Run.Palette2(obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  d <- a@Int.result[["bind"]]
  d <- d[cellname,]
  write.csv(as.matrix(d),
            paste0(out_dir,'/Palette_s3',paste0(mod[i]),'.csv'),
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
setwd('./Young_Age')
out_dir <- './Benchmarking/Supervised_using_modified_label/Label_miss/output/Ab-seq'

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in paths){
  adt.list[[i]] <- readRDS(paste0('./',i,'/adt.rds'))
  rna.list[[i]] <- readRDS(paste0('./',i,'/rna.rds'))
  meta.list[[i]] <- readRDS(paste0('./',i,'/meta.rds'))
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

meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_miss/output/Ab-seq/meta.rds')

cellname <- meta[,2]
rownames(meta) <- meta[,2]
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_miss5','label_miss10','label_miss20')
mod <- c('_miss5','_miss10','_miss20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,supervised = T,
                    group = 'celltype',
                    dims = list(1:50,1:50),# 使用哪些维数
                    method = c("PCA","PCA"))
  
  obj <- Find.Subspace3(obj,
                        modal = c("rna",'adt'),
                        lambda = list(0.8,0.8),
                        alpha = list(0,0),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:40,1:30),
                        cos.dims = list(1:50,1:50),
                        Angle.var = c(15,15),
                        max.Angle = c(50,50))
  
  a <- Run.Palette2(obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  d <- a@Int.result[["bind"]]
  d <- d[cellname,]
  write.csv(as.matrix(d),
            paste0(out_dir,'/Palette',paste0(mod[i]),'.csv'),
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
setwd('./Young_Age')
out_dir <- './Benchmarking/Supervised_using_modified_label/Label_miss/output/Ab-seq'

paths <- c('Age1','Age2','Age3','BM1','BM2','BM3')
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in paths){
  adt.list[[i]] <- readRDS(paste0('./',i,'/adt.rds'))
  rna.list[[i]] <- readRDS(paste0('./',i,'/rna.rds'))
  meta.list[[i]] <- readRDS(paste0('./',i,'/meta.rds'))
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


meta <- readRDS('./Benchmarking/Supervised_using_modified_label/Label_miss/output/Ab-seq/meta_s4.rds')

cellname <- meta[,2]
rownames(meta) <- meta[,2]
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_miss5','label_miss10','label_miss20')
mod <- c('_miss5','_miss10','_miss20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,supervised = T,
                    group = 'celltype',
                    dims = list(1:50,1:50),# 使用哪些维数
                    method = c("PCA","PCA"))
  
  obj <- Find.Subspace3(obj,
                        modal = c("rna",'adt'),
                        lambda = list(0.8,0.8),
                        alpha = list(0,0),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:40,1:30),
                        cos.dims = list(1:50,1:50),
                        Angle.var = c(15,15),
                        max.Angle = c(50,50))
  
  a <- Run.Palette2(obj,nn=2,nDims = 20,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')
  d <- a@Int.result[["bind"]]
  d <- d[cellname,]
  write.csv(as.matrix(d),
            paste0(out_dir,'/Palette_s4',paste0(mod[i]),'.csv'),
            row.names = FALSE)
}