library(bindSC)
library(Seurat)
library(Signac)
library(data.table)
library(Matrix)
library(matrixStats)

setwd('./MERFISH_ATAC')

atac <- readRDS('./ATAC/atac.rds')
atac_co <- readRDS('./ATAC/co100_RF.rds')

rna_co <- readRDS('./MERFISH/co100_RF.rds')

# get cluster for bindsc x, using all x features
atac_n <- RunTFIDF(atac)
rm(atac)
svd_res <- RunSVD(atac_n)
xc_obj=CreateSeuratObject(counts=atac_n)
xc_obj[['lsi']] <- svd_res
xc_obj <- FindNeighbors(xc_obj, dims = 1:30,reduction = 'lsi')
xc_obj <- FindClusters(xc_obj)
x_cluster = as.factor(paste0('x_',as.character(Idents(xc_obj))))

# get cluster for bindsc x, using all x features
x_obj=CreateSeuratObject(counts=atac_co,assay="x")
x_obj <- NormalizeData(x_obj)

# get cluster for bindsc y, using all y features (variable)
y_obj <- CreateSeuratObject(counts=rna_co,assay="y")
y_obj <- ScaleData(y_obj, features = rownames(y_obj))
VariableFeatures(object = y_obj) <- rownames(rna_co)
y_obj <- RunPCA(y_obj, features = VariableFeatures(object = y_obj),npcs = 50)
y_obj <- FindNeighbors(y_obj, dims = 1:30)
y_obj <- FindClusters(y_obj)
y_cluster = as.factor(paste0('y_',as.character(Idents(y_obj))))

## for Z0
z_obj = y_obj

## now gather all the actual inputs
x_input = z_obj@assays[[1]]@data

## now gather all the actual inputs
x_input = as.matrix(fread('./MERFISH/bdsc30_100RF.csv',data.table = F))

# y_input = atac_n
y_input <- as.matrix(t(svd_res@cell.embeddings[,1:30]))
z0_input = as.matrix(fread('./ATAC/bdsc30_100RF.csv',data.table = F))

rm(xc_obj,x_obj,y_obj,z_obj,rna_co,atac_co,atac_n)
gc()

rm(svd_res)

# start bindsc
res <- BiCCA( X = x_input ,
              Y =  y_input, 
              Z0 = z0_input, 
              X.clst = y_cluster,
              Y.clst = x_cluster,
              alpha = 0.5,lambda = 0.5,
              K = 20,
              temp.path  = "./",
              num.iteration = 50)

write.csv(as.matrix(rbind(res$u,res$r)),
          './Reduced_feat/RF100/bindsc.csv',
          row.names=FALSE)