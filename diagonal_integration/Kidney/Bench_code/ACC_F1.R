library(Seurat)
library(data.table)
library(Matrix)
library(matrixStats)

setwd('/home/server/sqy/data/mosaic/result/weak_link/kidney/')

method = c("Palette","stab","Multigrate","scMoMaT",
           "midas",'scVAEIT',"bindsc",'uniport',"UINMF",
            "fastMNN","Seurat","Harmony",'glue','scCon','SIMBA')

res <- list()
for(i in 1:length(method)){
    res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}
res[[13]] <- res[[13]][,-1]

meta <- readRDS('./meta.rds')

require(tidyverse)
evaluate <- function(true_lab, pred_lab) {
  "
    Conf: confusion matrix
    macro : macro F1-score
    F1 : F1-score per class
    Acc : accuracy
    PercUnl : percentage of unlabeled cells
    PopSize : number of cells per cell type
    "
  unique_true <- unlist(unique(true_lab))
  unique_pred <- unlist(unique(pred_lab))
  
  unique_all <- unique(c(unique_true, unique_pred))
  conf <- table(true_lab, pred_lab)
  pop_size <- rowSums(conf)
  
  conf_F1 <- table(true_lab, pred_lab)
  F1 <- vector()
  sum_acc <- 0
  
  for (i in c(1:length(unique_true))) {
    findLabel <- colnames(conf_F1) == row.names(conf_F1)[i]
    if (sum(findLabel)) {
      prec <- conf_F1[i, findLabel] / colSums(conf_F1)[findLabel]
      rec <- conf_F1[i, findLabel] / rowSums(conf_F1)[i]
      if (prec == 0 || rec == 0) {
        F1[i] <- 0
      } else {
        F1[i] <- (2 * prec * rec) / (prec + rec)
      }
      sum_acc <- sum_acc + conf_F1[i, findLabel]
    } else {
      F1[i] <- 0
    }
  }
  
  pop_size <- pop_size[pop_size > 0]
  names(F1) <- names(pop_size)
  macro_F1 <- mean(F1)
  total <- length(pred_lab)
  num_unlab <- sum(pred_lab == "unassigned") + sum(pred_lab == "Unassigned") + sum(pred_lab == "rand") + sum(pred_lab == "Unknown") + sum(pred_lab == "unknown") + sum(pred_lab == "Node") + sum(pred_lab == "ambiguous") + sum(pred_lab == "unassign")
  per_unlab <- num_unlab / total
  acc <- sum_acc / sum(conf_F1)
  
  result <- list(Acc = acc, Macro_F1 = macro_F1, Conf = conf, F1 = F1, True_labels = true_lab, Pred_labels = pred_lab)
  return(result)
}

meta_rna <- dplyr::filter(meta,modality == "RNA")
meta_codex <- dplyr::filter(meta,modality == "ATAC")

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

out1 <- data.frame(ACC = rep(0,length(method)),macro_F1 = rep(0,length(method)))
rownames(out1) <- method

pred_label_df1 <- as.data.frame(matrix(0,nrow = nrow(meta_rna),ncol = length(method)+1))
colnames(pred_label_df1) <- c('ground_truth',method)
rownames(pred_label_df1) <- rownames(meta_rna)
pred_label_df1$ground_truth <- meta_rna$label

for(i in 1:length(method)){
    
    if(method[i] != 'Seurat'){
       pred_tmp <- Trans.Knowledge(res[[i]][19986:nrow(meta),],
        res[[i]][1:19985,],
                      data.frame(meta_codex[,2]))
    
    pred_label_df1[,method[i]] <- pred_tmp[,1]
    score_tmp <- evaluate(meta_rna[,2],pred_tmp[,1])
    out1[method[i],'ACC'] <- score_tmp[[1]]
    out1[method[i],'macro_F1'] <- score_tmp[[2]] 
    }
}

out2 <- data.frame(ACC = rep(0,length(method)),macro_F1 = rep(0,length(method)))
rownames(out2) <- method

pred_label_df2 <- as.data.frame(matrix(0,nrow = nrow(meta_codex),ncol = length(method)+1))
colnames(pred_label_df2) <- c('ground_truth',method)
rownames(pred_label_df2) <- rownames(meta_codex)
pred_label_df2$ground_truth <- meta_codex$label

for(i in 1:length(method)){
    
    if(method[i] != 'Seurat'){
         pred_tmp <- Trans.Knowledge(res[[i]][1:19985,],res[[i]][19986:nrow(meta),],
                      data.frame(meta_rna[,2]))
    
    pred_label_df2[,method[i]] <- pred_tmp[,1]
    score_tmp <- evaluate(meta_codex[,2],pred_tmp[,1])
    out2[method[i],'ACC'] <- score_tmp[[1]]
    out2[method[i],'macro_F1'] <- score_tmp[[2]]   
    }
}

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

rna <- as.matrix(readRDS('./RNA/co.rds'))
atac <- as.matrix(readRDS('./ATAC/co.rds'))

run_Seurat = function(batches, meta, 
                      is.normalize = TRUE, 
                      nfeatures = 2000, 
                      vargs = NULL, 
                      reduction = c("cca", "rpca", "rlsi"),
                      k.filter = 200, 
                      out.npcs = 30)
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
    features <- SelectIntegrationFeatures(batch_list)
    for(i in 1:length(batch_list)) {
        batch_list[[i]] <- ScaleData(batch_list[[i]], features = features, verbose = FALSE)
        batch_list[[i]] <- RunPCA(batch_list[[i]], features = features, verbose = FALSE)
      }
  # run Seurat 
  t1 = Sys.time()
    anchors = FindTransferAnchors(
      reference = batch_list[[1]],
      query = batch_list[[2]],
      k.filter = NA,
      reduction = "pcaproject",
      reference.reduction = "pca",
      verbose = FALSE
    )
    
  gc()
  p1 = TransferData(anchorset = anchors,batch_list[[1]][['label']][,1])
   anchors = FindTransferAnchors(
      reference = batch_list[[2]],
      query = batch_list[[1]],
      k.filter = NA,
      reduction = "pcaproject",
      reference.reduction = "pca",
      verbose = FALSE
    ) 
  p2 = TransferData(anchorset = anchors,batch_list[[2]][['label']][,1])
  t2 = Sys.time()
  print(t2-t1)
  
  return(list(p1,p2))
}

seurat_res <- run_Seurat(list(rna,atac),
                         meta = meta,
                         is.normalize = T,
                         nfeatures = 1565,
                         reduction = 'cca',
                         out.npcs = 20)

score_tmp <- evaluate(meta_rna[,2],seurat_res[[2]][,1])
out1['Seurat','ACC'] <- score_tmp[[1]]
out1['Seurat','macro_F1'] <- score_tmp[[1]]

score_tmp <- evaluate(meta_codex[,2],seurat_res[[1]][,1])
out2['Seurat','ACC'] <- score_tmp[[1]]
out2['Seurat','macro_F1'] <- score_tmp[[1]]

rownames(out1) <- c("Palette","StabMap","Multigrate","scMoMaT",
                    "midas",'scVAEIT',"bindSC",'uniPort',"UINMF",
                    "fastMNN","Seurat","Harmony",'GLUE','scConfluence','SIMBA')

colnames(pred_label_df1) <- c('ground_truth',"Palette","StabMap","Multigrate","scMoMaT",
                              "midas",'scVAEIT',"bindSC",'uniPort',"UINMF",
                               "fastMNN","Seurat","Harmony",'GLUE','scConfluence','SIMBA')

rownames(out2) <- c("Palette","StabMap","Multigrate","scMoMaT",
                    "midas",'scVAEIT',"bindSC",'uniPort',"UINMF",
                    "fastMNN","Seurat","Harmony",'GLUE','scConfluence','SIMBA')

colnames(pred_label_df2) <- c('ground_truth',"Palette","StabMap","Multigrate","scMoMaT",
                              "midas",'scVAEIT',"bindSC",'uniPort',"UINMF",
                               "fastMNN","Seurat","Harmony",'GLUE','scConfluence','SIMBA')

saveRDS(out1,'./ACC_macro_f1_ATAC2RNA.rds')
saveRDS(out2,'./ACC_macro_f1_RNA2ATAC.rds')
