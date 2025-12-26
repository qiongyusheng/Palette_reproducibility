library(batchelor)
library(harmony)
library(rliger)
library(Seurat)
library(scater)
library(irlba)
library(BiocNeighbors)
library(SingleCellExperiment)
library(Matrix)
library(dplyr)
library(stats)

###################################################################################
### Seurat ###
###################################################################################
# k.filter (200 by default)
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

setwd('./BMMC_CITE_Multiome')
path <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite',
           './Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna.list <- list()
meta.list <- list()
for(i in 1:6){
  rna.list[[i]] <- readRDS(paste0(path[i],'/rna.rds'))
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
}

meta <- Reduce(rbind,meta.list)
head(meta)
meta$Batch <- paste0(meta$batch,'_',meta$modality)

srt_res <- run_Seurat(rna.list,meta,nfeatures = nrow(rna.list[[1]]),reduction = 'cca')

srt_clust <- kmeans(as.matrix(t(srt_res)),
                centers=length(unique(meta$celltype.l2)),nstart = 25)$cluster

meta$Seurat_cluster <- paste0('Seurat_',as.character(srt_clust))

meta_list <- split(meta, meta$Batch)
meta1.list <- meta_list[c(1,3,5)]
meta2.list <- meta_list[c(2,4,6)]

setwd('./Benchmarking/Supervised_using_modified_label/BMMC_s1_unimodal_label/output/Kmeans')

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
