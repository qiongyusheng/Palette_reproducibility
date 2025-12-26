source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

setwd("./imputation/input/TEA")
rna.list <- list()
atac.list <- list()
adt.list <- list()
meta.list <- list()
for(i in 2:5){
    if(i %in% c(2,4)){
        rna.list[[i-1]] <- readRDS(paste0("./B",i,"/rna.rds"))
    }else{
        rna.list[[i-1]] <- NULL
    }
    if(i %in% c(2,3)){
        atac.list[[i-1]] <- readRDS(paste0("./B",i,"/atac.rds"))
    }else{
        atac.list[[i-1]] <- NULL
    }
    if(i %in% c(2,5)){
        adt.list[[i-1]] <- readRDS(paste0("./B",i,"/adt.rds"))
    }else{
        adt.list[[i-1]] <- NULL
    }
  # rna.list[[i-1]] <- readRDS(paste0("./B",i,"/rna.rds"))
  meta.list[[i-1]] <- readRDS(paste0("./B",i,"/meta.rds"))
  # atac.list[[i-1]] <- readRDS(paste0("./B",i,"/atac.rds"))
  # adt.list[[i-1]] <- readRDS(paste0("./B",i,"/adt.rds"))
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

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'atac','adt'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:50,1:20),# 使用哪些维数
                         method = c("PCA",'LSI',"PCA"),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

options(future.globals.maxSize = 5000 * 1024^3)
musa_obj <- Find.Subspace3(musa_obj,
                          modal = c("rna",'atac','adt'),
                          lambda = list(0.8,0.5,0.8),
                          alpha = list(0,0,0),
                         sub.dims = list(1:40,1:30,1:20),
                         cos.dims = list(1:50,1:50,1:20),
                          Angle.var = c(15,20,15),
                         max.Angle = c(50,50,50))

a <- Run.Palette2(musa_obj,nn=10,nDims = 20,lambda = 0.5,origin.infer = T)

# B3 rna
B3_rna <- a@Assays[[7]]@data
write.csv(as.matrix(B3_rna),
          "./imputation/TEA_s1/Palette/B3_rna.csv",
          row.names = FALSE)

# B3 adt
B3_adt <- a@Assays[[8]]@data
write.csv(as.matrix(B3_adt),
          "./imputation/TEA_s1/Palette/B3_adt.csv",
          row.names = FALSE)

# B4 atac
B4_atac <- a@Assays[[9]]@data
Matrix::writeMM(B4_atac,
                "./imputation/TEA_s1/Palette/B4_atac.mtx")
write.csv(rownames(B4_atac),"./imputation/TEA_s1/Palette/peak.csv")
saveRDS(rownames(B4_atac),"./imputation/TEA_s1/Palette/peak.rds")

# B4 adt
B4_adt <- a@Assays[[10]]@data
write.csv(as.matrix(B4_adt),
          "./imputation/TEA_s1/Palette/B4_adt.csv",
          row.names = FALSE)

# B5 atac
B5_atac <- a@Assays[[11]]@data
Matrix::writeMM(B5_atac,
                "./imputation/TEA_s1/Palette/B5_atac.mtx")

# B5 rna
B5_rna <- a@Assays[[12]]@data
write.csv(as.matrix(B5_rna),
          "./imputation/TEA_s1/Palette/B5_rna.csv",
          row.names = FALSE)
