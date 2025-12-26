library(batchelor)
library(harmony)
library(Seurat)
library(scater)
library(irlba)
library(BiocNeighbors)
library(SingleCellExperiment)
library(Matrix)
library(umap)
library(dplyr)
library(Rcpp)
library(dplyr)
library(cowplot)

setwd('./')
rna <- as.matrix(readRDS('./RNA/co.rds'))
atac <- as.matrix(readRDS('./ATAC/co.rds'))
meta <- readRDS('./meta.rds')
out_dir <- './kidney/'

run_Seurat = function(batches, meta, is.normalize = TRUE, nfeatures = 2000, vargs = NULL, reduction = c("cca", "rpca", "rlsi"),
                      k.filter = 200, out.npcs = 30)
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  # modified rownames
  if (!is.null(vargs)) vargs = rownames(batch_seurat)[which(rownames(batches[[1]]) %in% vargs)]
  batch_list = SplitObject(batch_seurat, split.by = "modality")
  for(i in 1:length(batch_list)) {
    if(is.normalize == TRUE) {
      batch_list[[i]] <- NormalizeData(batch_list[[i]],normalization.method = "LogNormalize", scale.factor = 10000)
    }
    batch_list[[i]] <- FindVariableFeatures(batch_list[[i]], selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  }
  
  # run Seurat 
  t1 = Sys.time()
  if (is.null(vargs)) {
    features <- SelectIntegrationFeatures(batch_list)
    if (reduction == "rpca") {
      for(i in 1:length(batch_list)) {
        batch_list[[i]] <- ScaleData(batch_list[[i]], features = features, verbose = FALSE)
        batch_list[[i]] <- RunPCA(batch_list[[i]], features = features, verbose = FALSE)
      }
    }
    cell_anchors = FindIntegrationAnchors(object.list = batch_list, k.filter = k.filter, reduction = reduction)
  } else {
    if (reduction == "rpca") {
      for(i in 1:length(batch_list)) {
        batch_list[[i]] <- ScaleData(batch_list[[i]], features = vargs, verbose = FALSE)
        batch_list[[i]] <- RunPCA(batch_list[[i]], features = vargs, verbose = FALSE)
      }
    }
    cell_anchors = FindIntegrationAnchors(object.list = batch_list, k.filter = k.filter, reduction = reduction, anchor.features = vargs)
  }
  rm(batch_list, batch_seurat)
  gc()
  batch_correct = IntegrateData(anchorset = cell_anchors)
  t2 = Sys.time()
  print(t2-t1)
  
  DefaultAssay(batch_correct) = "integrated"
  batch_correct = ScaleData(object = batch_correct)
  batch_correct = RunPCA(object = batch_correct, npcs = out.npcs, verbose = FALSE)
  seurat_res = t(as.data.frame(batch_correct@reductions$pca@cell.embeddings))
  colnames(seurat_res) = rownames(meta)
  return(seurat_res)
}

seurat_res <- run_Seurat(list(rna,atac),
                         meta = meta,
                         is.normalize = T,
                         nfeatures = 1565,
                         reduction = 'cca',
                         out.npcs = 20)

write.csv(as.matrix(t(seurat_res)),
          paste0(out_dir,'Seurat.csv'),
          row.names = FALSE)

run_Harmony = function(batches, meta, is.normalize = TRUE, nfeatures = 2000, group.by.vars = "modality", theta = 2, vargs = NULL, 
                        out.npcs = 30)
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  # modified rownames
  if (!is.null(vargs)) vargs = rownames(batch_seurat)[which(rownames(batches[[1]]) %in% vargs)]
  if (is.normalize == TRUE) {
    batch_seurat <- NormalizeData(object = batch_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  batch_seurat <- FindVariableFeatures(object = batch_seurat, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  if (!is.null(vargs)) VariableFeatures(batch_seurat) = vargs
  batch_seurat <- ScaleData(object = batch_seurat)
  batch_seurat = RunPCA(object = batch_seurat, npcs = out.npcs, features = VariableFeatures(object = batch_seurat))
  
  #run Harmony
  t1 = Sys.time()
  batch_seurat = RunHarmony(object = batch_seurat, group.by.vars = group.by.vars, theta = theta, plot_convergence = TRUE, 
                             nclust = 50, max.iter.cluster = 100)
  t2 = Sys.time()
  print(t2-t1)
  
  harmony_res = t(as.data.frame(batch_seurat@reductions[["harmony"]]@cell.embeddings))
  colnames(harmony_res) = rownames(meta)
  return(harmony_res)
}

harmony_res <- run_Harmony(list(rna,atac),
                           meta = meta,
                           nfeatures = 1565,
                           out.npcs = 20)

write.csv(as.matrix(t(harmony_res)),
          paste0(out_dir,'Harmony.csv'),
          row.names = FALSE)

run_fastMNN = function(batches, meta, is.normalize = TRUE, nfeatures = 2000, vargs = NULL,
                    k = 20, out.npcs = 30)
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  # modified rownames
  if (!is.null(vargs)) vargs = rownames(batch_seurat)[which(rownames(batches[[1]]) %in% vargs)]
  if (is.normalize == TRUE) {
    batch_seurat <- NormalizeData(batch_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  batch_seurat <- FindVariableFeatures(object = batch_seurat, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  if (!is.null(vargs)) VariableFeatures(batch_seurat) = vargs
  vargenes <- batch_seurat@assays[["RNA"]]@var.features

  # log-normalized data matrices as input
  data_list = lapply(1:length(batches), function(i) batch_seurat@assays[["RNA"]]@data[vargenes, colnames(batches[[i]])])
  ##########################################################
  # run fastMNN
  t1 = Sys.time()
  out_mnn_total = do.call(batchelor::fastMNN, c(data_list, k = k, d = out.npcs))
  t2 = Sys.time()
  print(t2-t1)
  
  fastmnn_res = t(out_mnn_total@assays@data@listData[["reconstructed"]]@seed@components)
  colnames(fastmnn_res) = rownames(meta)
  return(fastmnn_res)
}

fastMNN_res <- run_fastMNN(list(rna,atac),
                           meta = meta,
                           nfeatures = 1565,
                           out.npcs = 20)

write.csv(as.matrix(t(fastMNN_res)),
          paste0(out_dir,'fastMNN.csv'),
          row.names = FALSE)
