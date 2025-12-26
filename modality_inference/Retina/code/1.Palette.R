source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

setwd('./imputation/input/Retina')

path <- c('./LGS1OD','./LGS1OS','./LGS2OD','./LGS2OS','./LGS3OD','./LGS3OS','./LVG1OD','./LVG1OS')

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
        atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
    }else{
       atac.list[[i]] <- NULL 
    }
    # rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds'))
    # atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
    meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
}

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

options(future.globals.maxSize = 5000 * 1024^3)
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- IDFLog_try(musa_obj,modal = "atac")
musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj,modal = "atac")

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'atac'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:50),# 使用哪些维数
                         method = c("PCA","LSI"),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

musa_obj <- Find.Subspace3(musa_obj,
                          modal = c("rna","atac"),
                          lambda = list(0.8,0.5),
                          alpha = list(0,0),
                         sub.dims = list(1:40,1:30),
                         cos.dims = list(1:50,1:50),
                          Angle.var = c(15,20),
                         max.Angle = c(50,50))

a <- Run.Palette2(musa_obj,nn=10,nDims = 20,lambda = 0.5,origin.infer = T)

# LGS2OD_atac
LGS2OD_atac <- a@Assays[[11]]@data
all.equal(rownames(meta.list[[3]]),colnames(LGS2OD_atac))
writeMM(LGS2OD_atac, 
        file="./imputation/Retina/Palette/LGS2OD_atac.mtx")

write.csv(data.frame(rownames=rownames(LGS2OD_atac)), 
          file="./imputation/Retina/Palette/peak.csv", 
          row.names=FALSE, 
          quote=FALSE)

# LGS2OS atac
LGS2OS_atac <- a@Assays[[12]]@data
all.equal(rownames(meta.list[[4]]),colnames(LGS2OS_atac))
writeMM(LGS2OS_atac, 
        file="./imputation/Retina/Palette/LGS2OS_atac.mtx")

# LGS3OD atac
LGS3OD_atac <- a@Assays[[13]]@data
all.equal(rownames(meta.list[[5]]),colnames(LGS3OD_atac))
writeMM(LGS3OD_atac, 
        file="./imputation/Retina/Palette/LGS3OD_atac.mtx")

# LGS3OS rna
LGS3OS_rna <- a@Assays[[14]]@data
all.equal(rownames(meta.list[[6]]),colnames(LGS3OS_rna))

write.csv(as.matrix(LGS3OS_rna),
          "./imputation/Retina/Palette/LGS3OS_rna.csv",
          row.names = FALSE)

# LVG1OD rna
LVG1OD_rna <- a@Assays[[15]]@data
all.equal(rownames(meta.list[[7]]),colnames(LVG1OD_rna))

write.csv(as.matrix(LVG1OD_rna),
          "./imputation/Retina/Palette/LVG1OD_rna.csv",
          row.names = FALSE)

# LVG1OS rna
LVG1OS_rna <- a@Assays[[16]]@data
all.equal(rownames(meta.list[[8]]),colnames(LVG1OS_rna))

write.csv(as.matrix(LVG1OS_rna),
          "./imputation/Retina/Palette/LVG1OS_rna.csv",
          row.names = FALSE)

