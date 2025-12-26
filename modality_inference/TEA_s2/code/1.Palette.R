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
    if(i %in% c(3,4)){
        rna.list[[i-1]] <- readRDS(paste0("./B",i,"/rna.rds"))
    }else{
        rna.list[[i-1]] <- NULL
    }
    if(i %in% c(2,3)){
        atac.list[[i-1]] <- readRDS(paste0("./B",i,"/atac.rds"))
    }else{
        atac.list[[i-1]] <- NULL
    }
    if(i %in% c(4,5)){
        adt.list[[i-1]] <- readRDS(paste0("./B",i,"/adt.rds"))
    }else{
        adt.list[[i-1]] <- NULL
    }
  # rna.list[[i-1]] <- readRDS(paste0("./B",i,"/rna.rds"))
  meta.list[[i-1]] <- readRDS(paste0("./B",i,"/meta.rds"))
  # atac.list[[i-1]] <- readRDS(paste0("./B",i,"/atac.rds"))
  # adt.list[[i-1]] <- readRDS(paste0("./B",i,"/adt.rds"))
}
musa_obj <- Create.Musa.Object(data.list = c(atac.list[1:2],rna.list[2:3],adt.list[3:4]), 
                               samples = c('B1','B2','B2','B3','B3','B4'), 
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
                         modal = c("rna",'adt','atac'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:20,1:50),# 使用哪些维数
                         method = c("PCA",'PCA','LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)
options(future.globals.maxSize = 50 * 1024^3)
musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'adt','atac'),
                           lambda = list(0.8,0.8,0.5),
                           supervised = F,
                           L2.norm = F,
                           sub.dims = list(1:40,1:20,1:30),
                           cos.dims = list(1:50,1:20,1:50),
                           Angle.var = c(15,15,20),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=10,nDims = 20,lambda = 0.5,emb_L2 = T,origin.infer = T)
mm <- data.frame(batch = c(rep("B1",6194),rep("B2",6381),rep("B3",6412),rep("B4",6530)),
                 celltype = c(as.character(meta.list[[1]][,2]),
                              as.character(meta.list[[2]][,2]),as.character(meta.list[[3]][,2]),
                              as.character(meta.list[[4]][,2]))
)
meta <- Reduce(rbind,meta.list)
rownames(mm) <- rownames(meta)


# B2 rna
B2_rna <- a@Assays[[7]]@data
all.equal(colnames(B2_rna),rownames(meta.list[[1]]))
write.csv(as.matrix(B2_rna),
          "./imputation/TEA_s2/Palette/B2_rna.csv",
          row.names = FALSE)

# B2 adt
B2_adt <- a@Assays[[8]]@data
all.equal(colnames(B2_adt),rownames(meta.list[[1]]))
write.csv(as.matrix(B2_adt),
          "./imputation/TEA_s2/Palette/B2_adt.csv",
          row.names = FALSE)

# B3 adt
B3_adt <- a@Assays[[9]]@data
all.equal(colnames(B3_adt),rownames(meta.list[[2]]))
write.csv(as.matrix(B3_adt),
          "./imputation/TEA_s2/Palette/B3_adt.csv",
          row.names = FALSE)

# B4 atac
B4_atac <- a@Assays[[10]]@data
all.equal(colnames(B4_atac),rownames(meta.list[[3]]))
Matrix::writeMM(B4_atac,"./imputation/TEA_s2/Palette/B4_atac.mtx")
write.csv(rownames(B4_atac),"./imputation/TEA_s2/Palette/peak.csv")
saveRDS(rownames(B4_atac),"./imputation/TEA_s2/Palette/peak.rds")

# B5 atac
B5_atac <- a@Assays[[11]]@data
all.equal(colnames(B5_atac),rownames(meta.list[[4]]))
Matrix::writeMM(B5_atac,"./imputation/TEA_s2/Palette/B5_atac.mtx")

# B5 rna
B5_rna <- a@Assays[[12]]@data
all.equal(colnames(B5_rna),rownames(meta.list[[4]]))
write.csv(as.matrix(B5_rna),
          "./imputation/TEA_s2/Palette/B5_rna.csv",
          row.names = FALSE)
