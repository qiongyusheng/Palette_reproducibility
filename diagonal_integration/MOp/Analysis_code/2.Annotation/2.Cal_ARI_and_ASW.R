library(data.table)
library(Matrix)
library(matrixStats)
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

meta <- readRDS('./All_anno.rds')
method = c('Palette','fastMNN_anno','Harmony_anno',
           'MIDAS_anno','scMoMaT_anno','scVAEIT_anno',
           'uniport_anno','Multigrate_anno',
           'UINMF_anno','Seurat_anno','Stab_anno','bindsc','maxfuse','scConfluence','celltype')

meta_atac <- dplyr::filter(meta,modality == "atac")
meta_merf <- dplyr::filter(meta,modality == "merfish")

require(tidyverse)

clust_labels_to_int = function(clust_labels) {
  uniq_labels = unique(clust_labels)
  clust_labels_int = rep(NA, length(clust_labels))
  
  for (ii in 1:length(clust_labels)) {
    clust_labels_int[ii] = match(clust_labels[ii], uniq_labels)
  }
  
  return(clust_labels_int)
}

ari = function(data = NULL, 
               meta = NULL, 
               label = 'celltype',
               methods = NULL,
               kmeans_nsim = 10) {
  ### The original x data (i.e., cells, features, and a column 'label' of cluster labels) 
  ### is under orig_fname_x.csv.
  ### The embeddings are in embed_fname_xi.csv where
  ### 0 <= i <= n_idx.
  ### The y data and embeddings are stored similarly.
  ### For each 0 <= i <= n_idx, calculate the following quantities:
  ###   1. ari_mix: one minus normalized adjusted random index between the k-means clustering result (with two clusters)
  ###               and the dataset index (x and y). This is a measure of mixing.
  ###   2. ari_clust: normalized adjusted random index between the k-means clustering result (with ground truth number of clusters) 
  ###                 and the ground truth cell type cluster labels.
  ###                 This is a measure of how well the embeddings are if used to do clustering 
  ###   3. ari_f1: 2 * ari * slt_clust / (ari_mix + ari_clust) --- F1 score of ari_mix and ari_clust
  ###
  ### Since the k-means implementation in R is highly dependent on random initializations, we do k-means for kmeans_nsim times
  ### and average the result
  ### Each one of ari_mix, ari_clust, ari_f1 is a length (n_idx+1) vector.
  ### Return cbind(ari_mix, ari_clust, ari_f1)
  
  ARI_MIN = -1
  ARI_MAX = 1
  
  clust_label <- clust_labels_to_int(meta[,label])
  n_clusters <- length(unique(clust_label))
  ll <- list()
   for(i in 1:length(methods)){
       ll[[i]] <- clust_labels_to_int(meta[,methods[i]])
   } 
    
  avg_ari <- 0
  max_ari <- 0
  res <- data.frame(ARI = rep(0,length(method)))
  for(i in 1:kmeans_nsim){
    est_clust_label = kmeans(data, centers=n_clusters)$cluster
      for(j in 1:length(methods)){
         ari_clust = (mclust::adjustedRandIndex(est_clust_label, ll[[j]]) - ARI_MIN)  / (ARI_MAX - ARI_MIN)
         res[j,1] <- res[j,1]+ ari_clust / kmeans_nsim
      }
  }
  
  return(res)
}


merf = readRDS('./MERFISH_ATAC/MERFISH/pca30.rds')
atac = readRDS("./MERFISH_ATAC/ATAC/SVD30.rds")

out <- data.frame(merf_ARI = rep(0,length(method)),atac_ARI = rep(0,length(method)))

for(i in 1:2){
    ari_tmp = ari(merf,meta_merf,'Palette_anno',method)
    out[,'merf_ARI'] <- ari_tmp[,1]
    ari_tmp = ari(atac,meta_atac,'Palette_anno',method)
    out[,'atac_ARI'] <- ari_tmp[,1]
}

rownames(out) <- method
out$avg <- (out[,1]+out[,2])/2

rownames(out) <- c('Palette','fastMNN','Harmony',
           'MIDAS','scMoMaT','scVAEIT',
           'uniPort','Multigrate',
           'UINMF','Seurat','StabMap','BindSC','MaxFuse','scConfluence','RAW')

saveRDS(out,'./anno_ARI.rds')

out <- data.frame(merf_SW = rep(0,length(method)),atac_SW = rep(0,length(method)))

for(i in 1:length(method)){
    asw <- silhouette_cpp(clust_labels_to_int(meta_merf[,method[i]]),
                      as.matrix(t(merf)))
    out[i,'merf_SW'] <- mean(asw)
    asw <- silhouette_cpp(clust_labels_to_int(meta_atac[,method[i]]),
                      as.matrix(t(atac)))
    out[i,'atac_SW'] <- mean(asw)
}

rownames(out) <- c('Palette','fastMNN','Harmony',
           'MIDAS','scMoMaT','scVAEIT',
           'uniPort','Multigrate',
           'UINMF','Seurat','StabMap','BindSC','MaxFuse','scConfluence','RAW')

rownames(out) <- method
out$avg <- (out[,1]+out[,2])/2

saveRDS(out,'./anno_ASW.rds')
