source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

setwd("./")
path <- c('./PBMC_ASAP/control_asap','./PBMC_ASAP/stim_asap',
          './PBMC_CITE/control_cite','./PBMC_CITE/stim_cite')
atac.list <- list()
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in 1:4){
    meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
    adt.list[[i]] <- readRDS(paste0(path[i],'/adt.rds'))
    if(i %in% c(1,2)){
        atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
    }
    if(i %in% c(3,4)){
        rna.list[[i-2]] <- readRDS(paste0(path[i],'/rna.rds'))
    }
}

musa_obj <- Create.Musa.Object(data.list = c(adt.list,atac.list,rna.list), 
                               samples = c('control_asap','stim_asap',
                                           'control_cite','stim_cite',
                                          'control_asap','stim_asap',
                                           'control_cite','stim_cite'), 
                               modals = c('adt','adt','adt','adt',
                                         'atac','atac','rna','rna'),
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

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'atac','adt'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:50,1:50),# 使用哪些维数
                         method = c("PCA","PCA",'LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

options(future.globals.maxSize = 50 * 1024^3)
musa_obj <- Find.Subspace3(musa_obj,
                          modal = c("rna",'atac','adt'),
                          lambda = list(0.8,0.5,0.8),
                          alpha = list(0,0,0),
                           joint = T,
                           L2.norm = F,
                         sub.dims = list(1:40,1:30,1:30),
                         cos.dims = list(1:50,1:50,1:50),
                          Angle.var = c(15,20,15),
                         max.Angle = c(50,50,50))

a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)


write.csv(as.matrix(a@Int.result[["bind"]]),
          "./Cross_condition_PBMC/output/Embedding/Palette.csv",
          row.names = FALSE)
