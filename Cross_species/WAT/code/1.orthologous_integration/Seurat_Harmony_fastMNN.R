library(batchelor)
library(conos)
library(harmony)
library(rliger)
library(Seurat)
library(scInt)
library(scater)
library(irlba)
library(BiocNeighbors)
library(SingleCellExperiment)
library(Matrix)
library(umap)
library(dplyr)
library(Rcpp)

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
  out_mnn_total = do.call(batchelor::fastMNN, c(data_list, k = k, d = out.npcs, auto.merge = T))
  t2 = Sys.time()
  print(t2-t1)
  
  fastmnn_res = t(out_mnn_total@assays@data@listData[["reconstructed"]]@seed@components)
  colnames(fastmnn_res) = rownames(meta)
  return(fastmnn_res)
}
run_Seurat = function(batches, meta, is.normalize = TRUE, nfeatures = 2000, vargs = NULL, reduction = c("cca", "rpca", "rlsi"),
                      k.filter = 200, out.npcs = 30)
{
  # preprocessing
  batch_seurat = CreateSeuratObject(do.call(cbind, batches), meta.data = meta, min.cells = 0, min.features = 0)
  # modified rownames
  if (!is.null(vargs)) vargs = rownames(batch_seurat)[which(rownames(batches[[1]]) %in% vargs)]
  batch_list = SplitObject(batch_seurat, split.by = "Batch")
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
run_Harmony = function(batches, meta, is.normalize = TRUE, nfeatures = 2000, group.by.vars = "Batch", theta = 2, vargs = NULL, 
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

is.normalize = TRUE
npcs = 20
harmony_var = c("Batch")
harmony_theta = c(2)

setwd('./human_mouse_WAT')
meta <- readRDS('./meta.rds')
human_homo <- as(as.matrix(readRDS('./human/rna_homo.rds')),'dgCMatrix')
mouse_homo <- as(as.matrix(readRDS('./mouse/rna_homo.rds')),'dgCMatrix')
X <- cbind(human_homo,mouse_homo)
rm(human_homo,mouse_homo)
dataset = lapply(unique(meta$Batch), function(x) X[, which(meta$Batch == x)])
names(dataset) = unique(meta$Batch)
print(identical(rownames(meta), Reduce('c', lapply(dataset, colnames))))

head(rownames(dataset[[1]]))
head(meta)

out_dir <- './'

res_fastMNN = run_fastMNN(dataset, meta, is.normalize = is.normalize,
                          vargs = rownames(dataset[[1]]), out.npcs = npcs)
write.csv(as.matrix(t(res_fastMNN)),
          paste0(out_dir,'fastMNN.csv'),
          row.names = FALSE)
res_Seurat = run_Seurat(dataset, meta, is.normalize = is.normalize, 
                        vargs = rownames(dataset[[1]]), out.npcs = npcs, reduction = "rpca")
write.csv(as.matrix(t(res_Seurat)),
          paste0(out_dir,'Seurat.csv'),
          row.names = FALSE)
res_harmony = run_Harmony(dataset, meta, is.normalize = is.normalize,
                          vargs = rownames(dataset[[1]]), group.by.vars = harmony_var,
                          theta = harmony_theta, out.npcs = npcs)
write.csv(as.matrix(t(res_harmony)),
          paste0(out_dir,'Harmony.csv'),
          row.names = FALSE)
