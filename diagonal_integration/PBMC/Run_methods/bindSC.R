library(bindSC)
library(Seurat)
library(data.table)
library(Matrix)
library(matrixStats)

setwd('./')

rna <- as(as.matrix(readRDS('./RNA/rna.rds')),'dgCMatrix')
cytof <- as(as.matrix(readRDS('./cytof/cytof.rds')),'dgCMatrix')
rna_co <- as.matrix(readRDS('./RNA/co.rds'))
cytof_co <- as.matrix(readRDS('./cytof/co.rds'))

rownames(rna) <- paste0('rna',rownames(rna))
rownames(cytof) <- paste0('cytof',rownames(cytof))

# get cluster for bindsc x, using all x features
xc_obj=CreateSeuratObject(counts=cytof,assay="x")
xc_obj <- NormalizeData(xc_obj)
xc_obj <- ScaleData(xc_obj, features = rownames(xc_obj))
xc_obj <- RunPCA(xc_obj, features = rownames(xc_obj))
xc_obj <- FindNeighbors(xc_obj, dims = 1:25)
xc_obj <- FindClusters(xc_obj)
x_cluster = as.factor(paste0('x_',as.character(Idents(xc_obj))))

# get cluster for bindsc x, using all x features
x_obj=CreateSeuratObject(counts=cytof_co,assay="x")
x_obj <- NormalizeData(x_obj)
x_obj <- ScaleData(x_obj, features = rownames(x_obj))# not used

# get cluster for bindsc y, using all y features (variable)
y_obj=CreateSeuratObject(counts=rna,assay="y")
y_obj <- NormalizeData(y_obj)
y_obj <- ScaleData(y_obj, features = rownames(y_obj))
VariableFeatures(object = y_obj) <- rownames(rna)
y_obj <- RunPCA(y_obj, features = VariableFeatures(object = y_obj))
y_obj <- FindNeighbors(y_obj, dims = 1:30)
y_obj <- FindClusters(y_obj)
y_cluster = as.factor(paste0('y_',as.character(Idents(y_obj))))

## for Z0
z_obj=CreateSeuratObject(counts=rna_co,assay="z")
z_obj <- NormalizeData(z_obj)
z_obj <- ScaleData(z_obj, features = rownames(z_obj))# not used

## now gather all the actual inputs
x_input = x_obj@assays$x@data
y_input = as.matrix(as.data.frame(y_obj@assays$y@data))
z0_input = z_obj@assays$z@data

rm(rna,rna_co,cytof,cytof_co,x_obj,y_obj,z_obj)
gc()

# start bindsc
res <- BiCCA( X = x_input ,
              Y =  y_input, 
              Z0 = z0_input, 
              X.clst = x_cluster,
              Y.clst = y_cluster,
              alpha = 0.5, 
              lambda = 0.5,
              K = 20,
              temp.path  = "./",
              num.iteration = 50)


write.csv(as.matrix(rbind(res$r,res$u)),
          './bindsc.csv',
          row.names=FALSE) 