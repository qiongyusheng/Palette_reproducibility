source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

setwd('./Young_Age')
out_dir <- './'

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

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'adt'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:50),# 使用哪些维数
                         method = c("PCA","PCA"),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt'),
                           lambda = list(0.8,0.8),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30),
                           cos.dims = list(1:50,1:50),
                           Angle.var = c(15,15),
                           max.Angle = c(50,50))

a <- Run.Palette2(musa_obj,nn=10,nDims = 20,lambda = 0.5, emb_L2 = T,origin.infer = T)

tmp <- a@Assays[[9]]@data

write.csv(as.matrix(tmp),
          "./ABseq/predict/palette/Age2rna.csv",
          row.names = FALSE)

tmp <- a@Assays[[10]]@data

write.csv(as.matrix(tmp),
          "./ABseq/predict/palette/BM3rna.csv",
          row.names = FALSE)

tmp <- a@Assays[[11]]@data

write.csv(as.matrix(tmp),
          "./ABseq/predict/palette/Age3adt.csv",
          row.names = FALSE)

tmp <- a@Assays[[12]]@data

write.csv(as.matrix(tmp),
          "./ABseq/predict/palette/BM2adt.csv",
          row.names = FALSE)
