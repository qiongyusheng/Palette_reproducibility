library(future)
library(Seurat)
library(data.table)
library(dplyr)
library(purrr)
library(Signac)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(Matrix)
library(matrixStats)
library(edgeR)
library(Rcpp)
sourceCpp("./Palette/utils.cpp")
Trans.Knowledge <- function(ref,
                            query,
                            ref.label = NULL,
                            k = 5){
  
  pred_label <- class::knn(ref,
                           query,
                           as.character(ref.label[,1]),
                           k = k, prob = TRUE) %>% as.character()
  
  out <- data.frame(predict_label = pred_label)
  
  return(out)
  
}

clust_labels_to_int = function(clust_labels) {
  uniq_labels = unique(clust_labels)
  clust_labels_int = rep(NA, length(clust_labels))
  
  for (ii in 1:length(clust_labels)) {
    clust_labels_int[ii] = match(clust_labels[ii], uniq_labels)
  }
  
  return(clust_labels_int)
}

meta <- readRDS('./meta.rds')

int <- as.matrix(fread('./Palette.csv',
                       data.table = F))
rownames(int) <- rownames(meta)
meta_merf <- dplyr::filter(meta,modality == 'merfish')
meta_atac <- dplyr::filter(meta,modality == 'atac')
asw <- silhouette_cpp(clust_labels_to_int(meta_merf$celltype),
                      as.matrix(t(int[rownames(meta_merf),])))
names(asw) <- rownames(meta_merf)
idx <- which(asw > 0)
asw_sub <- asw[idx]
asw_re <- asw[-idx]
int_ref <- int[names(asw_sub),]
int_query <- int[c(names(asw_re),rownames(meta_atac)),]
dim(int)
dim(int_ref)
dim(int_query)
meta_ref <- meta[rownames(int_ref),]
meta_query <- meta[rownames(int_query),]
pred <- Trans.Knowledge(int_ref,
                        int_query,
                        data.frame(meta_ref$celltype))
meta_new <- rbind(meta_ref,meta_query)
meta_new$Palette_anno <- c(meta_ref$celltype,pred[,1])
meta_new <- meta_new[rownames(meta),]
meta <- meta_new

int <- as.matrix(fread('./bindsc.csv',
                       data.table = F))
rownames(int) <- rownames(meta)
meta_merf <- dplyr::filter(meta,modality == 'merfish')
meta_atac <- dplyr::filter(meta,modality == 'atac')
asw <- silhouette_cpp(clust_labels_to_int(meta_merf$celltype),
                      as.matrix(t(int[rownames(meta_merf),])))
names(asw) <- rownames(meta_merf)
idx <- which(asw > 0)
asw_sub <- asw[idx]
asw_re <- asw[-idx]
int_ref <- int[names(asw_sub),]
int_query <- int[c(names(asw_re),rownames(meta_atac)),]
dim(int)
dim(int_ref)
dim(int_query)
meta_ref <- meta[rownames(int_ref),]
meta_query <- meta[rownames(int_query),]
pred <- Trans.Knowledge(int_ref,
                        int_query,
                        data.frame(meta_ref$celltype))
meta_new <- rbind(meta_ref,meta_query)
meta_new$bindsc_anno <- c(meta_ref$celltype,pred[,1])
meta_new <- meta_new[rownames(meta),]
meta <- meta_new


int <- as.matrix(fread('./fastMNN.csv',
                       data.table = F))
rownames(int) <- rownames(meta)
meta_merf <- dplyr::filter(meta,modality == 'merfish')
meta_atac <- dplyr::filter(meta,modality == 'atac')
asw <- silhouette_cpp(clust_labels_to_int(meta_merf$celltype),
                      as.matrix(t(int[rownames(meta_merf),])))
names(asw) <- rownames(meta_merf)
idx <- which(asw > 0)
asw_sub <- asw[idx]
asw_re <- asw[-idx]
int_ref <- int[names(asw_sub),]
int_query <- int[c(names(asw_re),rownames(meta_atac)),]
dim(int)
dim(int_ref)
dim(int_query)
meta_ref <- meta[rownames(int_ref),]
meta_query <- meta[rownames(int_query),]
pred <- Trans.Knowledge(int_ref,
                        int_query,
                        data.frame(meta_ref$celltype))
meta_new <- rbind(meta_ref,meta_query)
meta_new$fastMNN_anno <- c(meta_ref$celltype,pred[,1])
meta_new <- meta_new[rownames(meta),]
meta <- meta_new

int <- as.matrix(fread('./Harmony.csv',
                       data.table = F))
rownames(int) <- rownames(meta)
meta_merf <- dplyr::filter(meta,modality == 'merfish')
meta_atac <- dplyr::filter(meta,modality == 'atac')
asw <- silhouette_cpp(clust_labels_to_int(meta_merf$celltype),
                      as.matrix(t(int[rownames(meta_merf),])))
names(asw) <- rownames(meta_merf)
idx <- which(asw > 0)
asw_sub <- asw[idx]
asw_re <- asw[-idx]
int_ref <- int[names(asw_sub),]
int_query <- int[c(names(asw_re),rownames(meta_atac)),]
dim(int)
dim(int_ref)
dim(int_query)
meta_ref <- meta[rownames(int_ref),]
meta_query <- meta[rownames(int_query),]
pred <- Trans.Knowledge(int_ref,
                        int_query,
                        data.frame(meta_ref$celltype))
meta_new <- rbind(meta_ref,meta_query)
meta_new$Harmony_anno <- c(meta_ref$celltype,pred[,1])
meta_new <- meta_new[rownames(meta),]
meta <- meta_new

int <- as.matrix(fread('./midas.csv',
                       data.table = F))
rownames(int) <- rownames(meta)
meta_merf <- dplyr::filter(meta,modality == 'merfish')
meta_atac <- dplyr::filter(meta,modality == 'atac')
asw <- silhouette_cpp(clust_labels_to_int(meta_merf$celltype),
                      as.matrix(t(int[rownames(meta_merf),])))
names(asw) <- rownames(meta_merf)
idx <- which(asw > 0)
asw_sub <- asw[idx]
asw_re <- asw[-idx]
int_ref <- int[names(asw_sub),]
int_query <- int[c(names(asw_re),rownames(meta_atac)),]
dim(int)
dim(int_ref)
dim(int_query)
meta_ref <- meta[rownames(int_ref),]
meta_query <- meta[rownames(int_query),]
pred <- Trans.Knowledge(int_ref,
                        int_query,
                        data.frame(meta_ref$celltype))
meta_new <- rbind(meta_ref,meta_query)
meta_new$MIDAS_anno <- c(meta_ref$celltype,pred[,1])
meta_new <- meta_new[rownames(meta),]
meta <- meta_new

int <- as.matrix(fread('./scMoMaT.csv',
                       data.table = F,head = T))
rownames(int) <- rownames(meta)
meta_merf <- dplyr::filter(meta,modality == 'merfish')
meta_atac <- dplyr::filter(meta,modality == 'atac')
asw <- silhouette_cpp(clust_labels_to_int(meta_merf$celltype),
                      as.matrix(t(int[rownames(meta_merf),])))
names(asw) <- rownames(meta_merf)
idx <- which(asw > 0)
asw_sub <- asw[idx]
asw_re <- asw[-idx]
int_ref <- int[names(asw_sub),]
int_query <- int[c(names(asw_re),rownames(meta_atac)),]
dim(int)
dim(int_ref)
dim(int_query)
meta_ref <- meta[rownames(int_ref),]
meta_query <- meta[rownames(int_query),]
pred <- Trans.Knowledge(int_ref,
                        int_query,
                        data.frame(meta_ref$celltype))
meta_new <- rbind(meta_ref,meta_query)
meta_new$scMoMaT_anno <- c(meta_ref$celltype,pred[,1])
meta_new <- meta_new[rownames(meta),]
meta <- meta_new

int <- as.matrix(fread('./scVAEIT.csv',
                       data.table = F,head = T))
rownames(int) <- rownames(meta)
meta_merf <- dplyr::filter(meta,modality == 'merfish')
meta_atac <- dplyr::filter(meta,modality == 'atac')
asw <- silhouette_cpp(clust_labels_to_int(meta_merf$celltype),
                      as.matrix(t(int[rownames(meta_merf),])))
names(asw) <- rownames(meta_merf)
idx <- which(asw > 0)
asw_sub <- asw[idx]
asw_re <- asw[-idx]
int_ref <- int[names(asw_sub),]
int_query <- int[c(names(asw_re),rownames(meta_atac)),]
dim(int)
dim(int_ref)
dim(int_query)
meta_ref <- meta[rownames(int_ref),]
meta_query <- meta[rownames(int_query),]
pred <- Trans.Knowledge(int_ref,
                        int_query,
                        data.frame(meta_ref$celltype))
meta_new <- rbind(meta_ref,meta_query)
meta_new$scVAEIT_anno <- c(meta_ref$celltype,pred[,1])
meta_new <- meta_new[rownames(meta),]
meta <- meta_new

int <- as.matrix(fread('./stab.csv',
                       data.table = F))
rownames(int) <- rownames(meta)
meta_merf <- dplyr::filter(meta,modality == 'merfish')
meta_atac <- dplyr::filter(meta,modality == 'atac')
asw <- silhouette_cpp(clust_labels_to_int(meta_merf$celltype),
                      as.matrix(t(int[rownames(meta_merf),])))
names(asw) <- rownames(meta_merf)
idx <- which(asw > 0)
asw_sub <- asw[idx]
asw_re <- asw[-idx]
int_ref <- int[names(asw_sub),]
int_query <- int[c(names(asw_re),rownames(meta_atac)),]
dim(int)
dim(int_ref)
dim(int_query)
meta_ref <- meta[rownames(int_ref),]
meta_query <- meta[rownames(int_query),]
pred <- Trans.Knowledge(int_ref,
                        int_query,
                        data.frame(meta_ref$celltype))
meta_new <- rbind(meta_ref,meta_query)
meta_new$StabMap_anno <- c(meta_ref$celltype,pred[,1])
meta_new <- meta_new[rownames(meta),]
meta <- meta_new

int <- as.matrix(fread('./uniport.csv',
                       data.table = F,head = T))
rownames(int) <- rownames(meta)
meta_merf <- dplyr::filter(meta,modality == 'merfish')
meta_atac <- dplyr::filter(meta,modality == 'atac')
asw <- silhouette_cpp(clust_labels_to_int(meta_merf$celltype),
                      as.matrix(t(int[rownames(meta_merf),])))
names(asw) <- rownames(meta_merf)
idx <- which(asw > 0)
asw_sub <- asw[idx]
asw_re <- asw[-idx]
int_ref <- int[names(asw_sub),]
int_query <- int[c(names(asw_re),rownames(meta_atac)),]
dim(int)
dim(int_ref)
dim(int_query)
meta_ref <- meta[rownames(int_ref),]
meta_query <- meta[rownames(int_query),]
pred <- Trans.Knowledge(int_ref,
                        int_query,
                        data.frame(meta_ref$celltype))
meta_new <- rbind(meta_ref,meta_query)
meta_new$uniport_anno <- c(meta_ref$celltype,pred[,1])
meta_new <- meta_new[rownames(meta),]
meta <- meta_new

int <- as.matrix(fread('./Multigrate.csv',
                       data.table = F,head = T))
rownames(int) <- rownames(meta)
meta_merf <- dplyr::filter(meta,modality == 'merfish')
meta_atac <- dplyr::filter(meta,modality == 'atac')
asw <- silhouette_cpp(clust_labels_to_int(meta_merf$celltype),
                      as.matrix(t(int[rownames(meta_merf),])))
names(asw) <- rownames(meta_merf)
idx <- which(asw > 0)
asw_sub <- asw[idx]
asw_re <- asw[-idx]
int_ref <- int[names(asw_sub),]
int_query <- int[c(names(asw_re),rownames(meta_atac)),]
dim(int)
dim(int_ref)
dim(int_query)
meta_ref <- meta[rownames(int_ref),]
meta_query <- meta[rownames(int_query),]
pred <- Trans.Knowledge(int_ref,
                        int_query,
                        data.frame(meta_ref$celltype))
meta_new <- rbind(meta_ref,meta_query)
meta_new$Multigrate_anno <- c(meta_ref$celltype,pred[,1])
meta_new <- meta_new[rownames(meta),]
meta <- meta_new

int <- as.matrix(fread('./UINMF.csv',
                       data.table = F))
rownames(int) <- rownames(meta)
meta_merf <- dplyr::filter(meta,modality == 'merfish')
meta_atac <- dplyr::filter(meta,modality == 'atac')
asw <- silhouette_cpp(clust_labels_to_int(meta_merf$celltype),
                      as.matrix(t(int[rownames(meta_merf),])))
names(asw) <- rownames(meta_merf)
idx <- which(asw > 0)
asw_sub <- asw[idx]
asw_re <- asw[-idx]
int_ref <- int[names(asw_sub),]
int_query <- int[c(names(asw_re),rownames(meta_atac)),]
dim(int)
dim(int_ref)
dim(int_query)
meta_ref <- meta[rownames(int_ref),]
meta_query <- meta[rownames(int_query),]
pred <- Trans.Knowledge(int_ref,
                        int_query,
                        data.frame(meta_ref$celltype))
meta_new <- rbind(meta_ref,meta_query)
meta_new$UINMF_anno <- c(meta_ref$celltype,pred[,1])
meta_new <- meta_new[rownames(meta),]
meta <- meta_new

int <- as.matrix(fread('./Seurat.csv',
                       data.table = F))
rownames(int) <- rownames(meta)
meta_merf <- dplyr::filter(meta,modality == 'merfish')
meta_atac <- dplyr::filter(meta,modality == 'atac')
asw <- silhouette_cpp(clust_labels_to_int(meta_merf$celltype),
                      as.matrix(t(int[rownames(meta_merf),])))
names(asw) <- rownames(meta_merf)
idx <- which(asw > 0)
asw_sub <- asw[idx]
asw_re <- asw[-idx]
int_ref <- int[names(asw_sub),]
int_query <- int[c(names(asw_re),rownames(meta_atac)),]
dim(int)
dim(int_ref)
dim(int_query)
meta_ref <- meta[rownames(int_ref),]
meta_query <- meta[rownames(int_query),]
pred <- Trans.Knowledge(int_ref,
                        int_query,
                        data.frame(meta_ref$celltype))
meta_new <- rbind(meta_ref,meta_query)
meta_new$Seurat_anno <- c(meta_ref$celltype,pred[,1])
meta_new <- meta_new[rownames(meta),]
meta <- meta_new

int <- as.matrix(fread('./stab.csv',
                       data.table = F))
rownames(int) <- rownames(meta)
meta_merf <- dplyr::filter(meta,modality == 'merfish')
meta_atac <- dplyr::filter(meta,modality == 'atac')
asw <- silhouette_cpp(clust_labels_to_int(meta_merf$celltype),
                      as.matrix(t(int[rownames(meta_merf),])))
names(asw) <- rownames(meta_merf)
idx <- which(asw > 0)
asw_sub <- asw[idx]
asw_re <- asw[-idx]
int_ref <- int[names(asw_sub),]
int_query <- int[c(names(asw_re),rownames(meta_atac)),]
dim(int)
dim(int_ref)
dim(int_query)
meta_ref <- meta[rownames(int_ref),]
meta_query <- meta[rownames(int_query),]
pred <- Trans.Knowledge(int_ref,
                        int_query,
                        data.frame(meta_ref$celltype))
meta_new <- rbind(meta_ref,meta_query)
meta_new$Stab_anno <- c(meta_ref$celltype,pred[,1])
meta_new <- meta_new[rownames(meta),]
meta <- meta_new

int <- as.matrix(fread('/home/server/sqy/data/mosaic/result/NC_revise_1/weak_link/Mop_MERFISH_ATAC/maxfuse.csv',
                       data.table = F,head = T))
rownames(int) <- rownames(meta)
meta_merf <- dplyr::filter(meta,modality == 'merfish')
meta_atac <- dplyr::filter(meta,modality == 'atac')
asw <- silhouette_cpp(clust_labels_to_int(meta_merf$celltype),
                      as.matrix(t(int[rownames(meta_merf),])))
names(asw) <- rownames(meta_merf)
idx <- which(asw > 0)
asw_sub <- asw[idx]
asw_re <- asw[-idx]
int_ref <- int[names(asw_sub),]
int_query <- int[c(names(asw_re),rownames(meta_atac)),]
dim(int)
dim(int_ref)
dim(int_query)
meta_ref <- meta[rownames(int_ref),]
meta_query <- meta[rownames(int_query),]
pred <- Trans.Knowledge(int_ref,
                        int_query,
                        data.frame(meta_ref$celltype))
meta_new <- rbind(meta_ref,meta_query)
meta_new$maxfuse <- c(meta_ref$celltype,pred[,1])
meta_new <- meta_new[rownames(meta),]
meta <- meta_new

int <- as.matrix(fread('/home/server/sqy/data/mosaic/result/NC_revise_1/weak_link/Mop_MERFISH_ATAC/scCon.csv',
                       data.table = F,head = T))
rownames(int) <- rownames(meta)
meta_merf <- dplyr::filter(meta,modality == 'merfish')
meta_atac <- dplyr::filter(meta,modality == 'atac')
asw <- silhouette_cpp(clust_labels_to_int(meta_merf$celltype),
                      as.matrix(t(int[rownames(meta_merf),])))
names(asw) <- rownames(meta_merf)
idx <- which(asw > 0)
asw_sub <- asw[idx]
asw_re <- asw[-idx]
int_ref <- int[names(asw_sub),]
int_query <- int[c(names(asw_re),rownames(meta_atac)),]
dim(int)
dim(int_ref)
dim(int_query)
meta_ref <- meta[rownames(int_ref),]
meta_query <- meta[rownames(int_query),]
pred <- Trans.Knowledge(int_ref,
                        int_query,
                        data.frame(meta_ref$celltype))
meta_new <- rbind(meta_ref,meta_query)
meta_new$scConfluence <- c(meta_ref$celltype,pred[,1])
meta_new <- meta_new[rownames(meta),]
meta <- meta_new

saveRDS(meta,'./All_anno.rds')
