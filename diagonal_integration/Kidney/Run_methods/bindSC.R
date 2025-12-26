library(bindSC)
library(Seurat)
library(Signac)
library(data.table)
library(Matrix)
library(matrixStats)

setwd('./')

rna <- readRDS('./RNA/X.rds')
rna_co <- readRDS('./RNA/co.rds')
atac <- readRDS('./ATAC/atac.rds')
atac_co <- readRDS('./ATAC/co.rds')

rownames(rna) <- paste0('rna',rownames(rna))

# get cluster for bindsc x, using all x features
atac_n <- RunTFIDF(atac)
xc_obj=CreateSeuratObject(counts=atac,assay="x")
xc_obj@assays[[1]]@scale.data <- as.matrix(atac_n)
xc_obj <- RunSVD(xc_obj, features = rownames(xc_obj))
xc_obj <- FindNeighbors(xc_obj, dims = 2:30,reduction = 'lsi')
xc_obj <- FindClusters(xc_obj)
x_cluster = as.factor(paste0('x_',as.character(Idents(xc_obj))))

# get cluster for bindsc x, using all x features
x_obj=CreateSeuratObject(counts=atac_co,assay="x")
x_obj <- NormalizeData(x_obj)

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
x_input = z_obj@assays[[1]]@data
# y_input = as.matrix(atac_n)
y_input <- as.matrix(t(xc_obj[['lsi']]@cell.embeddings[,1:30]))
z0_input = x_obj@assays[[1]]@data

rm(atac_n,xc_obj,x_obj,y_obj,z_obj,rna,atac,rna_co,atac_co)
gc()

# start bindsc
res <- BiCCA( X = x_input ,
              Y =  y_input, 
              Z0 = z0_input, 
              X.clst = y_cluster,
              Y.clst = x_cluster,
              alpha = 0.5,lambda = 0.5,
              K = 20,
              temp.path  = "./kidney",
              num.iteration = 50)

meta <- readRDS('./meta.rds')

write.csv(as.matrix(rbind(res$u,res$r)),
          './kidney/bindsc.csv',
          row.names=FALSE)
