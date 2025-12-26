######################### unsupervised ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

setwd('./BMMC_CITE_Multiome')
options(future.globals.maxSize = 50 * 1024^3)
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


musa_obj <- Create.Musa.Object(data.list = c(rna1.list,
                                             rna2.list,
                                             adt.list,
                                             atac.list), 
                               samples = c('CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna","rna",
                                          "adt",'adt','adt',
                                          "atac","atac","atac"),
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
meta.list <- c(meta1.list,meta2.list)
cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))
mm <- data.frame(batch = c(rep("CITE1",5227),rep("CITE2",4978),rep("CITE3",6106),
                           rep("Multiome1",6224),rep("Multiome2",6740),rep("Multiome3",4279)),
                 celltype.l2 = c(as.character(meta1.list[[1]][,5]),
                                 as.character(meta1.list[[2]][,5]),
                                 as.character(meta1.list[[3]][,5]),
                                 as.character(meta2.list[[1]][,5]),
                                 as.character(meta2.list[[2]][,5]),
                                 as.character(meta2.list[[3]][,5])),
                 celltype.l1 = c(as.character(meta1.list[[1]][,4]),
                                 as.character(meta1.list[[2]][,4]),
                                 as.character(meta1.list[[3]][,4]),
                                 as.character(meta2.list[[1]][,4]),
                                 as.character(meta2.list[[2]][,4]),
                                 as.character(meta2.list[[3]][,4])),
                 modality = c(rep('CITE-seq',16311),rep('Multiome',17243))
)
out_dir <- './Benchmarking/Reference-based_integration/output/Reference'
d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_unsupervised.csv'),
          row.names = FALSE)


a <- Run.Palette2(musa_obj,nn=2,nDims = 30,lambda = 0.5,emb_L2 = T)
meta.list <- c(meta1.list,meta2.list)
cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))
mm <- data.frame(batch = c(rep("CITE1",5227),rep("CITE2",4978),rep("CITE3",6106),
                           rep("Multiome1",6224),rep("Multiome2",6740),rep("Multiome3",4279)),
                 celltype.l2 = c(as.character(meta1.list[[1]][,5]),
                                 as.character(meta1.list[[2]][,5]),
                                 as.character(meta1.list[[3]][,5]),
                                 as.character(meta2.list[[1]][,5]),
                                 as.character(meta2.list[[2]][,5]),
                                 as.character(meta2.list[[3]][,5])),
                 celltype.l1 = c(as.character(meta1.list[[1]][,4]),
                                 as.character(meta1.list[[2]][,4]),
                                 as.character(meta1.list[[3]][,4]),
                                 as.character(meta2.list[[1]][,4]),
                                 as.character(meta2.list[[2]][,4]),
                                 as.character(meta2.list[[3]][,4])),
                 modality = c(rep('CITE-seq',16311),rep('Multiome',17243))
)
out_dir <- './Benchmarking/Reference-based_integration/output/Reference'
d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_unsupervised_30.csv'),
          row.names = FALSE)
rownames(mm) <- cellname
saveRDS(mm,paste0(out_dir,'/meta.rds'))

meta.list <- c(meta1.list,meta2.list)
cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))
mm <- data.frame(batch = c(rep("CITE1",5227),rep("CITE2",4978),rep("CITE3",6106),
                           rep("Multiome1",6224),rep("Multiome2",6740),rep("Multiome3",4279)),
                 celltype.l2 = c(as.character(meta1.list[[1]][,5]),
                                 as.character(meta1.list[[2]][,5]),
                                 as.character(meta1.list[[3]][,5]),
                                 as.character(meta2.list[[1]][,5]),
                                 as.character(meta2.list[[2]][,5]),
                                 as.character(meta2.list[[3]][,5])),
                 celltype.l1 = c(as.character(meta1.list[[1]][,4]),
                                 as.character(meta1.list[[2]][,4]),
                                 as.character(meta1.list[[3]][,4]),
                                 as.character(meta2.list[[1]][,4]),
                                 as.character(meta2.list[[2]][,4]),
                                 as.character(meta2.list[[3]][,4])),
                 modality = c(rep('CITE-seq',16311),rep('Multiome',17243))
)
out_dir <- './Benchmarking/Reference-based_integration/output/Reference'
rownames(mm) <- cellname
saveRDS(mm,paste0(out_dir,'/meta.rds'))


######################### supervised ################################
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
                               samples = c('CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3',
                                           'CITE1','CITE2','CITE3',
                                           'Multiome1','Multiome2','Multiome3'), 
                               modals = c("rna","rna","rna",
                                          "rna","rna","rna",
                                          "adt",'adt','adt',
                                          "atac","atac","atac"),
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

mm <- data.frame(batch = c(rep("CITE1",5227),rep("CITE2",4978),rep("CITE3",6106),
                           rep("Multiome1",6224),rep("Multiome2",6740),rep("Multiome3",4279)),
                 ID = c(as.character(meta1.list[[1]][,1]),
                        as.character(meta1.list[[2]][,1]),
                        as.character(meta1.list[[3]][,1]),
                        as.character(meta2.list[[1]][,1]),
                        as.character(meta2.list[[2]][,1]),
                        as.character(meta2.list[[3]][,1])),
                 celltype.l2 = c(as.character(meta1.list[[1]][,5]),
                                 as.character(meta1.list[[2]][,5]),
                                 as.character(meta1.list[[3]][,5]),
                                 as.character(meta2.list[[1]][,5]),
                                 as.character(meta2.list[[2]][,5]),
                                 as.character(meta2.list[[3]][,5])),
                 celltype.l1 = c(as.character(meta1.list[[1]][,4]),
                                 as.character(meta1.list[[2]][,4]),
                                 as.character(meta1.list[[3]][,4]),
                                 as.character(meta2.list[[1]][,4]),
                                 as.character(meta2.list[[2]][,4]),
                                 as.character(meta2.list[[3]][,4])),
                 modality = c(rep('CITE-seq',16311),rep('Multiome',17243))
)

rownames(mm) <- mm[,2]
mm <- mm[musa_obj@Meta$cell.name,]
musa_obj@Meta[['celltype']] <- mm$celltype.l2

all.equal(musa_obj@Meta$cell.name,mm[,2])
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

meta.list <- c(meta1.list,meta2.list)
cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))
mm <- data.frame(batch = c(rep("CITE1",5227),rep("CITE2",4978),rep("CITE3",6106),
                           rep("Multiome1",6224),rep("Multiome2",6740),rep("Multiome3",4279)),
                 celltype.l2 = c(as.character(meta1.list[[1]][,5]),
                                 as.character(meta1.list[[2]][,5]),
                                 as.character(meta1.list[[3]][,5]),
                                 as.character(meta2.list[[1]][,5]),
                                 as.character(meta2.list[[2]][,5]),
                                 as.character(meta2.list[[3]][,5])),
                 celltype.l1 = c(as.character(meta1.list[[1]][,4]),
                                 as.character(meta1.list[[2]][,4]),
                                 as.character(meta1.list[[3]][,4]),
                                 as.character(meta2.list[[1]][,4]),
                                 as.character(meta2.list[[2]][,4]),
                                 as.character(meta2.list[[3]][,4])),
                 modality = c(rep('CITE-seq',16311),rep('Multiome',17243))
)
out_dir <- './Benchmarking/Reference-based_integration/output/Reference'
d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_supervised.csv'),
          row.names = FALSE)



a <- Run.Palette2(musa_obj,nn=2,nDims = 30,lambda = 0.5, supervised = T,emb_L2 = T, group = 'celltype')

meta.list <- c(meta1.list,meta2.list)
cellname <- Reduce(c,lapply(meta.list,function(x){
  x <- as.character(x[,1])
  x
}))
mm <- data.frame(batch = c(rep("CITE1",5227),rep("CITE2",4978),rep("CITE3",6106),
                           rep("Multiome1",6224),rep("Multiome2",6740),rep("Multiome3",4279)),
                 celltype.l2 = c(as.character(meta1.list[[1]][,5]),
                                 as.character(meta1.list[[2]][,5]),
                                 as.character(meta1.list[[3]][,5]),
                                 as.character(meta2.list[[1]][,5]),
                                 as.character(meta2.list[[2]][,5]),
                                 as.character(meta2.list[[3]][,5])),
                 celltype.l1 = c(as.character(meta1.list[[1]][,4]),
                                 as.character(meta1.list[[2]][,4]),
                                 as.character(meta1.list[[3]][,4]),
                                 as.character(meta2.list[[1]][,4]),
                                 as.character(meta2.list[[2]][,4]),
                                 as.character(meta2.list[[3]][,4])),
                 modality = c(rep('CITE-seq',16311),rep('Multiome',17243))
)
out_dir <- './Benchmarking/Reference-based_integration/output/Reference'
d <- a@Int.result[["bind"]]
d <- d[cellname,]
write.csv(as.matrix(d),
          paste0(out_dir,'/Palette_supervised_30.csv'),
          row.names = FALSE)
