library(Seurat)
library(cowplot)
library(dplyr)
library(data.table)
library(tibble)
library(Matrix)
library(SeuratData)
library(cowplot)

setwd('./immune_mtx/')
ref <- readMM('./matrix.mtx')
gene <- fread('./genes.tsv',data.table = F,head = F)
barcode <- fread('./barcodes.tsv',data.table = F,head = F)
rownames(ref) <- gene[,1]
colnames(ref) <- barcode[,1]

meta <- fread('./CensusImmune-BoneMarrow-10x_cell_type_2020-03-12.csv', # https://explore.data.humancellatlas.org/projects/cc95ff89-2e68-4a08-a234-480eca21ce79/project-matrices
          data.table = F)
head(meta)

meta_unique <- meta %>%
  group_by(barcode) %>%
  filter(n() == 1) %>%
  ungroup()
meta_unique <- as.data.frame(meta_unique)

nrow(meta_unique)
rownames(meta_unique) <- meta_unique$barcode
co_cell <- intersect(colnames(ref),rownames(meta_unique))
ref <- ref[,co_cell]
meta_unique <- meta_unique[co_cell,]
nrow(meta_unique)

meta <- meta_unique[,c(5,6)]
colnames(meta) <- c('celltype','barcode')
head(meta)

setwd('/home/server/sqy/data/mosaic')
rna_c <- readMM("./RAW/NIPS_CITE_Multiome_BMMC/CITE/rna.mtx")
adt <- readMM("./RAW/NIPS_CITE_Multiome_BMMC/CITE/adt.mtx")
gene_c <- fread("./RAW/NIPS_CITE_Multiome_BMMC/CITE/gene.csv",data.table = F)
p <- fread("./RAW/NIPS_CITE_Multiome_BMMC/CITE/protein.csv",data.table = F)
m_c <- fread("./RAW/NIPS_CITE_Multiome_BMMC/CITE/meta.csv",data.table = F)

rna_c <- t(rna_c)
adt <- t(adt)
gene_c <- gene_c[-1,]
p <- p[-1,]

rownames(rna_c) <- gene_c[,1]
rownames(adt) <- p[,1]
colnames(rna_c) <- colnames(adt) <- m_c[,1]

colnames(m_c)[1] <- 'barcode'
rownames(m_c) <- m_c[,1]

co_gene <- intersect(rownames(ref),rownames(rna_c))
ref_c <- ref[co_gene,]
rna_c <- rna_c[co_gene,]
querys <- CreateSeuratObject(rna_c,meta.data = m_c)
querys <- NormalizeData(querys)
query_list <- SplitObject(querys,split.by = 'batch')[1:3]
rm(querys)

obj <- CreateSeuratObject(ref_c,meta.data = meta)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)

for(i in 1:length(query_list)){
    query_list[[i]] <- FindVariableFeatures(query_list[[i]]) %>% ScaleData() %>% RunPCA()
    query_list[[i]] <- FindVariableFeatures(query_list[[i]],nfeatures = 5000)
    features <- SelectIntegrationFeatures(list(obj,query_list[[i]]))
    
    query_list[[i]] <- ScaleData(query_list[[i]], features = features, verbose = FALSE)
    query_list[[i]] <- RunPCA(query_list[[i]], features = features, verbose = FALSE)
    
    obj <- ScaleData(obj, features = features, verbose = FALSE)
    obj <- RunPCA(obj, features = features, verbose = FALSE)
    
    anchors = FindTransferAnchors(
      reference = obj,
      query = query_list[[i]],
      k.filter = NA,
      reduction = "pcaproject",
      reference.reduction = "pca",
      verbose = FALSE
    )
    predicted.labels = TransferData(anchorset = anchors,obj$celltype)
    query_list[[i]] <- AddMetaData(object = query_list[[i]], metadata = predicted.labels)
}

m_c_pred_list <- list()
for(i in 1:length(query_list)){
    m_c_pred_list[[i]] <- query_list[[i]]@meta.data[,c(colnames(m_c),'predicted.id')]
    }
m_c_pred <- Reduce(rbind,m_c_pred_list)

rna_m <- readMM("./RAW/NIPS_CITE_Multiome_BMMC/Multiome/rna.mtx")
# atac <- readMM("./RAW/NIPS_CITE_Multiome_BMMC/Multiome/atac.mtx")
gene_m <- fread("./RAW/NIPS_CITE_Multiome_BMMC/Multiome/gene.csv",data.table = F)
# peak <- fread("./RAW/NIPS_CITE_Multiome_BMMC/Multiome/peak.csv",data.table = F)
m_m <- fread("./RAW/NIPS_CITE_Multiome_BMMC/Multiome/meta.csv",data.table = F)

rna_m <- t(rna_m)
# atac <- t(atac)
gene_m <- gene_m[-1,]
# peak <- peak[-1,]

rownames(rna_m) <- gene_m[,1]
# rownames(atac) <- peak[,1]
colnames(rna_m) <- m_m[,1]
colnames(m_m)[1] <- 'barcode'
rownames(m_m) <- m_m[,1]



co_gene <- intersect(rownames(ref),rownames(rna_m))
ref_m <- ref[co_gene,]
rna_m <- rna_m[co_gene,]
querys <- CreateSeuratObject(rna_m,meta.data = m_m)
querys <- NormalizeData(querys)
query_list <- SplitObject(querys,split.by = 'batch')[1:3]
rm(querys,ref_c)

obj <- CreateSeuratObject(ref_m,meta.data = meta)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)

for(i in 1:length(query_list)){
    query_list[[i]] <- FindVariableFeatures(query_list[[i]]) %>% ScaleData() %>% RunPCA()
    query_list[[i]] <- FindVariableFeatures(query_list[[i]],nfeatures = 5000)
    features <- SelectIntegrationFeatures(list(obj,query_list[[i]]))
    
    query_list[[i]] <- ScaleData(query_list[[i]], features = features, verbose = FALSE)
    query_list[[i]] <- RunPCA(query_list[[i]], features = features, verbose = FALSE)
    
    obj <- ScaleData(obj, features = features, verbose = FALSE)
    obj <- RunPCA(obj, features = features, verbose = FALSE)
    
    anchors = FindTransferAnchors(
      reference = obj,
      query = query_list[[i]],
      k.filter = NA,
      reduction = "pcaproject",
      reference.reduction = "pca",
      verbose = FALSE
    )
    predicted.labels = TransferData(anchorset = anchors,obj$celltype)
    query_list[[i]] <- AddMetaData(object = query_list[[i]], metadata = predicted.labels)
}

m_m_pred_list <- list()
for(i in 1:length(query_list)){
    m_m_pred_list[[i]] <- query_list[[i]]@meta.data[,c(colnames(m_m),'predicted.id')]
    }
m_m_pred <- Reduce(rbind,m_m_pred_list)

meta_c <- m_c_pred[,c('barcode','batch',"cell_type",'l1_cell_type','l2_cell_type','predicted.id')]
meta_m<- m_m_pred[,c('barcode','batch',"cell_type",'l1_cell_type','l2_cell_type','predicted.id')]
colnames(meta_m) <- colnames(meta_c) <- c('barcode','batch',"celltype.origin",'celltype.l1','celltype.l2','predicted.id')
meta_m$modality <- 'Multiome'
meta_c$modality <- 'CITE'
rownames(meta_m) <- meta_m[,1]
rownames(meta_c) <- meta_c[,1]

setwd('./Benchmarking/Supervised_using_modified_label/BMMC_s1_unimodal_label/output/ReLabel')

b <- c('s1d1','s1d2','s1d3')
meta1.list <- list()
meta2.list <- list()
for(i in 1:3){
    meta1.list[[i]] <- dplyr::filter(meta_c,batch == b[i])
    meta2.list[[i]] <- dplyr::filter(meta_m,batch == b[i])
}

mm <- rbind(meta1.list[[2]],meta1.list[[3]],meta2.list[[1]],meta2.list[[2]],meta2.list[[3]])
mm$batch <- c(rep("CITE2",4978),rep("CITE3",6106),
              rep("Multiome1",6224),rep("Multiome2",6740),
              rep("Multiome3",4279))
rownames(mm) <- mm[,1]
saveRDS(mm ,'./BMMC/meta_drop1.rds')

mm <- rbind(meta1.list[[1]],meta1.list[[3]],
            meta2.list[[1]],meta2.list[[2]],meta2.list[[3]])
mm$batch <- c(rep("CITE1",5227),rep("CITE3",6106),
              rep("Multiome1",6224),
              rep("Multiome2",6740),rep("Multiome3",4279))
rownames(mm) <- mm[,1]
saveRDS(mm ,'./BMMC/meta_drop2.rds')

mm <- rbind(meta1.list[[1]],meta1.list[[2]],meta1.list[[3]],meta2.list[[2]],meta2.list[[3]])
mm$batch <- c(rep("CITE1",5227),rep("CITE2",4978),rep("CITE3",6106),
              rep("Multiome2",6740),rep("Multiome3",4279))
rownames(mm) <- mm[,1]
saveRDS(mm ,'./BMMC/meta_drop4.rds')

mm <- rbind(meta1.list[[1]],meta1.list[[2]],meta1.list[[3]],meta2.list[[1]],meta2.list[[3]])
mm$batch <- c(rep("CITE1",5227),rep("CITE2",4978),rep("CITE3",6106),
              rep("Multiome1",6224),
              rep("Multiome3",4279))
rownames(mm) <- mm[,1]
saveRDS(mm ,'./BMMC/meta_drop5.rds')

mm <- rbind(meta1.list[[1]],meta1.list[[2]],meta1.list[[3]],meta2.list[[1]],meta2.list[[2]])
mm$batch <- c(rep("CITE1",5227),rep("CITE2",4978),rep("CITE3",6106),
              rep("Multiome1",6224),
              rep("Multiome2",6740))
rownames(mm) <- mm[,1]
saveRDS(mm ,'./BMMC/meta_drop6.rds')
