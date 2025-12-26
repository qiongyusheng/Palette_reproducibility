source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
options(future.globals.maxSize = 50 * 1024^3)
setwd('./imputation/input/BMMC_s')

path <- c('CITE1rna','CITE1adt','CITE2rna','CITE2adt','CITE3',
          'Multiome1','Multiome2rna','Multiome2atac','Multiome3rna','Multiome3atac')
rna.list <- list()
adt.list <- list()
atac.list <- list()
meta.list <- list()
for(i in 1:10){
  
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
  if(i %in% c(1,3,5,6,7,9)){
    rna.list[[path[i]]] <- readRDS(paste0(path[i],'/rna.rds'))
  }
  
  if(i %in% c(2,4,5)){
    adt.list[[path[i]]] <- readRDS(paste0(path[i],'/adt.rds'))
  }
  if(i %in% c(6,8,10)){
    atac.list[[path[i]]] <- readRDS(paste0(path[i],'/atac.rds'))
  }
}

names(rna.list) <- NULL
names(adt.list) <- NULL
names(atac.list) <- NULL

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

options(future.globals.maxSize = 5000 * 1024^3)
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
                           alpha = list(0,0,0),
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))

a <- Run.Palette2(musa_obj,nn=10,nDims = 20,lambda = 0.5,origin.infer = T)

# cite1 rna
cite1_rna = a@Assays[[18]]@data
write.csv(as.matrix(cite1_rna),
          "./imputation/BMMC_s/Palette/cite1rna.csv",
          row.names = FALSE)

# cite2 rna
cite2_rna = a@Assays[[20]]@data
write.csv(as.matrix(cite2_rna),
          "./imputation/BMMC_s/Palette/cite2rna.csv",
          row.names = FALSE)

# cite1 adt
cite1_adt = a@Assays[[13]]@data
write.csv(as.matrix(cite1_adt),
          "./imputation/BMMC_s/Palette/cite1adt.csv",
          row.names = FALSE)

# cite2 adt
cite2_adt = a@Assays[[15]]@data
write.csv(as.matrix(cite2_adt),
          "./imputation/BMMC_s/Palette/cite2adt.csv",
          row.names = FALSE)

# mul2 atac
mul2_atac <- a@Assays[[24]]@data
mul2_atac <- as(mul2_atac,'dgCMatrix')
Matrix::writeMM(mul2_atac,
                "./imputation/BMMC_s/Palette/mul2atac.mtx")
write.csv(rownames(mul2_atac),"./imputation/BMMC_s/Palette/peak.csv")
saveRDS(rownames(mul2_atac),"./imputation/BMMC_s/Palette/peak.rds")

# mul2 rna
mul2_rna = a@Assays[[27]]@data
write.csv(as.matrix(mul2_rna),
          "./imputation/BMMC_s/Palette/mul2rna.csv",
          row.names = FALSE)

# mul3 atac
mul3_atac <- a@Assays[[26]]@data
mul3_atac <- as(mul3_atac,'dgCMatrix')
Matrix::writeMM(mul3_atac,
                "./imputation/BMMC_s/Palette/mul3atac.mtx")

# mul3 rna
mul3_rna = a@Assays[[29]]@data
write.csv(as.matrix(mul3_rna),
          "./imputation/BMMC_s/Palette/mul3rna.csv",
          row.names = FALSE)
