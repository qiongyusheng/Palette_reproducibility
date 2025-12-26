######################### ADT ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised.csv',data.table = F)))
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
colnames(latent_data) <- rownames(meta)

setwd('/home/server/sqy/data/mosaic/input/BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite')
ref <- list()
meta_ref.list <- list()
for(i in 1:3){
  ref[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta_ref.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}

setwd('/home/server/sqy/data/mosaic/input/BMMC_CITE_Multiome')
path1 <- c('./CITE/s2d1_cite','./CITE/s2d4_cite','./CITE/s2d5_cite')
q <- list()
meta_q.list <- list()
for(i in 1:3){
  q[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
  # meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$modality,'-',meta_q.list[[i]]$batch)
  meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$batch,'-',meta_q.list[[i]]$modality)
  meta_q.list[[i]]$modality <- 'ADT'
}

res <- Palette_map(c(ref,q),
                   modal = c('adt','adt','adt',
                             'adt','adt','adt'),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's2d1_cite','s2d4_cite','s2d5_cite'),
                   q_or_batch = c(rep('ref',3),rep('query',3)),
                   modal_order = c('adt'),
                   method = 'Harmony',
                   method_out.npcs = 40,
                   normalize_method = c('CLR'),
                   latent = latent_data,Bi_sPCA = T)

result_harmony <- run_Harmony(list(res[[1]],latent_data),
                              cell_name = c(colnames(res[[1]]),colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)

res <- Palette_map(c(ref,q),
                   modal = c('adt','adt','adt',
                             'adt','adt','adt'),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's2d1_cite','s2d4_cite','s2d5_cite'),
                   q_or_batch = c(rep('ref',3),rep('query',3)),
                   modal_order = c('adt'),
                   method = 'fastMNN',
                   method_out.npcs = 40,
                   normalize_method = c('CLR'),
                   latent = latent_data,Bi_sPCA = T)

result_fastMNN <- run_fastMNN(list(res[[1]],latent_data),
                              cell_name = c(colnames(res[[1]]),colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)


latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised_30.csv',
                                 data.table = F)))
colnames(latent_data) <- rownames(meta)

res <- Palette_map(c(ref,q),
                   modal = c('adt','adt','adt',
                             'adt','adt','adt'),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's2d1_cite','s2d4_cite','s2d5_cite'),
                   q_or_batch = c(rep('ref',3),rep('query',3)),
                   modal_order = c('adt'),
                   method = 'Seurat',
                   method_out.npcs = 40,
                   normalize_method = c('CLR'),
                   latent = latent_data,Bi_sPCA = T)

result_Seurat <- run_Seurat(list(res[[1]],latent_data),
                            cell_name = c(colnames(res[[1]]),colnames(latent_data)),
                            dims = 1:25,
                            is.normalize = F,out.npcs = NA)

meta_ref <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
meta_ref$RQ <- 'reference'
meta_q <- Reduce(rbind,meta_q.list)
meta_q$RQ <- 'query'
meta <- rbind(meta_ref,meta_q)

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

pred_label_l1_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_harmony <- evaluate(meta[rownames(meta_q),4],pred_label_l1_harmony[,1])
meta$pred_l1_harmony <- c(meta_ref[,4],pred_label_l1_harmony[,1])

pred_label_l2_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_harmony <- evaluate(meta[rownames(meta_q),5],pred_label_l2_harmony[,1])
meta$pred_l2_harmony <- c(meta_ref[,5],pred_label_l2_harmony[,1])

pred_label_l1_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_fastMNN <- evaluate(meta[rownames(meta_q),4],pred_label_l1_fastMNN[,1])
meta$pred_l1_fastMNN <- c(meta_ref[,4],pred_label_l1_fastMNN[,1])

pred_label_l2_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_fastMNN <- evaluate(meta[rownames(meta_q),5],pred_label_l2_fastMNN[,1])
meta$pred_l2_fastMNN <- c(meta_ref[,5],pred_label_l2_fastMNN[,1])

pred_label_l1_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),4]))
scores_l1_Seurat <- evaluate(meta[rownames(meta_q),4],pred_label_l1_Seurat[,1])
meta$pred_l1_Seurat <- c(meta_ref[,4],pred_label_l1_Seurat[,1])

pred_label_l2_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),5]))
scores_l2_Seurat <- evaluate(meta[rownames(meta_q),5],pred_label_l2_Seurat[,1])
meta$pred_l2_Seurat <- c(meta_ref[,5],pred_label_l2_Seurat[,1])

meta <- meta[c(rownames(meta_ref),rownames(meta_q)),]

score <- data.frame(Acc = c(scores_l1_harmony[[1]],scores_l2_harmony[[1]],
                            scores_l1_fastMNN[[1]],scores_l2_fastMNN[[1]],
                            scores_l1_Seurat[[1]],scores_l2_Seurat[[1]]),
                    Macro_F1 = c(scores_l1_harmony[[2]],scores_l2_harmony[[2]],
                                 scores_l1_fastMNN[[2]],scores_l2_fastMNN[[2]],
                                 scores_l1_Seurat[[2]],scores_l2_Seurat[[2]]))
rownames(score) <- c('l1_harmony','l2_harmony',
                     'l1_fastMNN','l2_fastMNN',
                     'l1_Seurat','l2_Seurat')

out_dir <- './Benchmarking/Reference-based_integration/output/Result/adt/'
write.csv(as.matrix(t(result_harmony[,rownames(meta)])),
          paste0(out_dir,'Palette_harmony_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_fastMNN[,rownames(meta)])),
          paste0(out_dir,'Palette_fastMNN_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_Seurat[,rownames(meta)])),
          paste0(out_dir,'Palette_Seurat_supervised.csv'),
          row.names = FALSE)

saveRDS(meta,paste0(out_dir,'meta_supervised.rds'))
saveRDS(meta_q,paste0(out_dir,'meta_query_supervised.rds'))
saveRDS(meta_ref,paste0(out_dir,'meta_ref_supervised.rds'))
saveRDS(score,paste0(out_dir,'score_supervised.rds'))

######################### ADT + ATAC ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised.csv',data.table = F)))
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
colnames(latent_data) <- rownames(meta)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite',
           './Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
# rna.list <- list()
adt.list <- list()
atac.list <- list()
meta_ref.list <- list()

for(i in 1:6){
  if(i %in% c(1,2,3)){
    adt.list[[i]] <-  readRDS(paste0(path1[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[i-3]] <-  readRDS(paste0(path1[i],'/atac.rds'))
  }
  meta_ref.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}

ref <- c(adt.list,atac.list)

setwd('/home/server/sqy/data/mosaic/input/BMMC_CITE_Multiome')
path1 <- c('./CITE/s2d1_cite','./CITE/s2d4_cite','./CITE/s2d5_cite',
           './Multiome/s2d1_multi','./Multiome/s2d4_multi','./Multiome/s2d5_multi')

# rna.list <- list()
adt.list <- list()
atac.list <- list()
meta_q.list <- list()
for(i in 1:6){
  if(i %in% c(1,2,3)){
    adt.list[[i]] <-  readRDS(paste0(path1[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[i-3]] <-  readRDS(paste0(path1[i],'/atac.rds'))
  }
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
  # meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$modality,'-',meta_q.list[[i]]$batch)
  meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$batch,'-',meta_q.list[[i]]$modality)
  
  if(i %in% c(1,2,3)){
    meta_q.list[[i]]$modality <- 'ADT'
  }else{
    meta_q.list[[i]]$modality <- 'ATAC'
  }
  
}
q <- c(adt.list,atac.list)
rm(adt.list,atac.list)
res <- Palette_map(c(ref,q),
                   modal = c(rep('adt',3),rep('atac',3),
                             rep('adt',3),rep('atac',3)),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',6),rep('query',6)),
                   modal_order = c('adt','atac'),
                   method = 'Harmony',
                   method_out.npcs = 40,
                   normalize_method = c('CLR','TF-IDF'),
                   latent = latent_data,Bi_sPCA = T)

result_harmony <- run_Harmony(list(res[[1]],
                                   res[[2]],
                                   latent_data),
                              cell_name = c(colnames(res[[1]]),
                                            colnames(res[[2]]),
                                            colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)

res <- Palette_map(c(ref,q),
                   modal = c(rep('adt',3),rep('atac',3),
                             rep('adt',3),rep('atac',3)),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',6),rep('query',6)),
                   modal_order = c('adt','atac'),
                   method = 'fastMNN',
                   method_out.npcs = 40,
                   normalize_method = c('CLR','TF-IDF'),
                   latent = latent_data,Bi_sPCA = T)

result_fastMNN <- run_fastMNN(list(res[[1]],
                                   res[[2]],
                                   latent_data),
                              cell_name = c(colnames(res[[1]]),
                                            colnames(res[[2]]),
                                            colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)
latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised_30.csv',
                                 data.table = F)))
colnames(latent_data) <- rownames(meta)

res <- Palette_map(c(ref,q),
                   modal = c(rep('adt',3),rep('atac',3),
                             rep('adt',3),rep('atac',3)),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',6),rep('query',6)),
                   modal_order = c('adt','atac'),
                   method = 'Seurat',
                   method_out.npcs = 40,
                   normalize_method = c('CLR','TF-IDF'),
                   latent = latent_data,Bi_sPCA = T)

result_Seurat <- run_Seurat(list(res[[1]],
                                 res[[2]],
                                 latent_data),
                            cell_name = c(colnames(res[[1]]),
                                          colnames(res[[2]]),
                                          colnames(latent_data)),
                            dims = 1:25,
                            is.normalize = F,out.npcs = NA)

meta_ref <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
meta_ref$RQ <- 'reference'
meta_q <- Reduce(rbind,meta_q.list)
meta_q$RQ <- 'query'
meta <- rbind(meta_ref,meta_q)

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

pred_label_l1_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_harmony <- evaluate(meta[rownames(meta_q),4],pred_label_l1_harmony[,1])
meta$pred_l1_harmony <- c(meta_ref[,4],pred_label_l1_harmony[,1])

pred_label_l2_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_harmony <- evaluate(meta[rownames(meta_q),5],pred_label_l2_harmony[,1])
meta$pred_l2_harmony <- c(meta_ref[,5],pred_label_l2_harmony[,1])

pred_label_l1_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_fastMNN <- evaluate(meta[rownames(meta_q),4],pred_label_l1_fastMNN[,1])
meta$pred_l1_fastMNN <- c(meta_ref[,4],pred_label_l1_fastMNN[,1])

pred_label_l2_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_fastMNN <- evaluate(meta[rownames(meta_q),5],pred_label_l2_fastMNN[,1])
meta$pred_l2_fastMNN <- c(meta_ref[,5],pred_label_l2_fastMNN[,1])

pred_label_l1_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),4]))
scores_l1_Seurat <- evaluate(meta[rownames(meta_q),4],pred_label_l1_Seurat[,1])
meta$pred_l1_Seurat <- c(meta_ref[,4],pred_label_l1_Seurat[,1])

pred_label_l2_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),5]))
scores_l2_Seurat <- evaluate(meta[rownames(meta_q),5],pred_label_l2_Seurat[,1])
meta$pred_l2_Seurat <- c(meta_ref[,5],pred_label_l2_Seurat[,1])

meta <- meta[c(rownames(meta_ref),rownames(meta_q)),]

score <- data.frame(Acc = c(scores_l1_harmony[[1]],scores_l2_harmony[[1]],
                            scores_l1_fastMNN[[1]],scores_l2_fastMNN[[1]],
                            scores_l1_Seurat[[1]],scores_l2_Seurat[[1]]),
                    Macro_F1 = c(scores_l1_harmony[[2]],scores_l2_harmony[[2]],
                                 scores_l1_fastMNN[[2]],scores_l2_fastMNN[[2]],
                                 scores_l1_Seurat[[2]],scores_l2_Seurat[[2]]))
rownames(score) <- c('l1_harmony','l2_harmony',
                     'l1_fastMNN','l2_fastMNN',
                     'l1_Seurat','l2_Seurat')

out_dir <- './Benchmarking/Reference-based_integration/output/Result/adt_atac/'
write.csv(as.matrix(t(result_harmony[,rownames(meta)])),
          paste0(out_dir,'Palette_harmony_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_fastMNN[,rownames(meta)])),
          paste0(out_dir,'Palette_fastMNN_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_Seurat[,rownames(meta)])),
          paste0(out_dir,'Palette_Seurat_supervised.csv'),
          row.names = FALSE)
saveRDS(meta,paste0(out_dir,'meta_supervised.rds'))
saveRDS(meta_q,paste0(out_dir,'meta_query_supervised.rds'))
saveRDS(meta_ref,paste0(out_dir,'meta_ref_supervised.rds'))
saveRDS(score,paste0(out_dir,'score_supervised.rds'))


######################### ATAC ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised.csv',data.table = F)))
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
colnames(latent_data) <- rownames(meta)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
ref <- list()
meta_ref.list <- list()
for(i in 1:3){
  ref[[i]] <- readRDS(paste0(path1[i],'/atac.rds'))
  meta_ref.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}

setwd('/home/server/sqy/data/mosaic/input/BMMC_CITE_Multiome')
path1 <- c('./Multiome/s2d1_multi','./Multiome/s2d4_multi','./Multiome/s2d5_multi')
q <- list()
meta_q.list <- list()
for(i in 1:3){
  q[[i]] <- readRDS(paste0(path1[i],'/atac.rds'))
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
  # meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$modality,'-',meta_q.list[[i]]$batch)
  meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$batch,'-',meta_q.list[[i]]$modality)
  meta_q.list[[i]]$modality <- 'ATAC'
}

res <- Palette_map(c(ref,q),
                   modal = c('atac','atac','atac',
                             'atac','atac','atac'),
                   batch = c('s1d1_multi','s1d2_multi','s1d3_multi',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',3),rep('query',3)),
                   modal_order = c('atac'),
                   method = 'Harmony',
                   method_out.npcs = 40,
                   normalize_method = c('TF-IDF'),
                   latent = latent_data,Bi_sPCA = T)

result_harmony <- run_Harmony(list(res[[1]],latent_data),
                              cell_name = c(colnames(res[[1]]),colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)

res <- Palette_map(c(ref,q),
                   modal = c('atac','atac','atac',
                             'atac','atac','atac'),
                   batch = c('s1d1_multi','s1d2_multi','s1d3_multi',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',3),rep('query',3)),
                   modal_order = c('atac'),
                   method = 'fastMNN',
                   method_out.npcs = 40,
                   normalize_method = c('TF-IDF'),
                   latent = latent_data,Bi_sPCA = T)

result_fastMNN <- run_fastMNN(list(res[[1]],latent_data),
                              cell_name = c(colnames(res[[1]]),colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)

latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised_30.csv',
                                 data.table = F)))
colnames(latent_data) <- rownames(meta)

res <- Palette_map(c(ref,q),
                   modal = c('atac','atac','atac',
                             'atac','atac','atac'),
                   batch = c('s1d1_multi','s1d2_multi','s1d3_multi',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',3),rep('query',3)),
                   modal_order = c('atac'),
                   method = 'Seurat',
                   method_out.npcs = 40,
                   normalize_method = c('TF-IDF'),
                   latent = latent_data,Bi_sPCA = T)

result_Seurat <- run_Seurat(list(res[[1]],latent_data),
                            cell_name = c(colnames(res[[1]]),colnames(latent_data)),
                            dims = 1:25,
                            is.normalize = F,out.npcs = NA)

meta_ref <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
meta_ref$RQ <- 'reference'
meta_q <- Reduce(rbind,meta_q.list)
meta_q$RQ <- 'query'
meta <- rbind(meta_ref,meta_q)

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
pred_label_l1_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_harmony <- evaluate(meta[rownames(meta_q),4],pred_label_l1_harmony[,1])
meta$pred_l1_harmony <- c(meta_ref[,4],pred_label_l1_harmony[,1])

pred_label_l2_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_harmony <- evaluate(meta[rownames(meta_q),5],pred_label_l2_harmony[,1])
meta$pred_l2_harmony <- c(meta_ref[,5],pred_label_l2_harmony[,1])

pred_label_l1_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_fastMNN <- evaluate(meta[rownames(meta_q),4],pred_label_l1_fastMNN[,1])
meta$pred_l1_fastMNN <- c(meta_ref[,4],pred_label_l1_fastMNN[,1])

pred_label_l2_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_fastMNN <- evaluate(meta[rownames(meta_q),5],pred_label_l2_fastMNN[,1])
meta$pred_l2_fastMNN <- c(meta_ref[,5],pred_label_l2_fastMNN[,1])


pred_label_l1_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),4]))
scores_l1_Seurat <- evaluate(meta[rownames(meta_q),4],pred_label_l1_Seurat[,1])
meta$pred_l1_Seurat <- c(meta_ref[,4],pred_label_l1_Seurat[,1])

pred_label_l2_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),5]))
scores_l2_Seurat <- evaluate(meta[rownames(meta_q),5],pred_label_l2_Seurat[,1])
meta$pred_l2_Seurat <- c(meta_ref[,5],pred_label_l2_Seurat[,1])

meta <- meta[c(rownames(meta_ref),rownames(meta_q)),]

score <- data.frame(Acc = c(scores_l1_harmony[[1]],scores_l2_harmony[[1]],
                            scores_l1_fastMNN[[1]],scores_l2_fastMNN[[1]],
                            scores_l1_Seurat[[1]],scores_l2_Seurat[[1]]),
                    Macro_F1 = c(scores_l1_harmony[[2]],scores_l2_harmony[[2]],
                                 scores_l1_fastMNN[[2]],scores_l2_fastMNN[[2]],
                                 scores_l1_Seurat[[2]],scores_l2_Seurat[[2]]))
rownames(score) <- c('l1_harmony','l2_harmony',
                     'l1_fastMNN','l2_fastMNN',
                     'l1_Seurat','l2_Seurat')

out_dir <- './Benchmarking/Reference-based_integration/output/Result/atac/'
write.csv(as.matrix(t(result_harmony[,rownames(meta)])),
          paste0(out_dir,'Palette_harmony_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_fastMNN[,rownames(meta)])),
          paste0(out_dir,'Palette_fastMNN_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_Seurat[,rownames(meta)])),
          paste0(out_dir,'Palette_Seurat_supervised.csv'),
          row.names = FALSE)
saveRDS(meta,paste0(out_dir,'meta_supervised.rds'))
saveRDS(meta_q,paste0(out_dir,'meta_query_supervised.rds'))
saveRDS(meta_ref,paste0(out_dir,'meta_ref_supervised.rds'))
saveRDS(score,paste0(out_dir,'score_supervised.rds'))

######################### CITE ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised.csv',data.table = F)))
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
colnames(latent_data) <- rownames(meta)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite')
rna1.list <- list()
adt.list <- list()
meta_ref.list <- list()
for(i in 1:3){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta_ref.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}
ref <- c(rna1.list,adt.list)

setwd('/home/server/sqy/data/mosaic/input/BMMC_CITE_Multiome')
path1 <- c('./CITE/s2d1_cite','./CITE/s2d4_cite','./CITE/s2d5_cite')
rna1.list <- list()
adt.list <- list()
meta_q.list <- list()
for(i in 1:3){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  adt.list[[i]] <- readRDS(paste0(path1[i],'/adt.rds'))
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}
q <- c(rna1.list,adt.list)
rm(rna1.list,adt.list)

res <- Palette_map(c(ref,q),
                   modal = c('rna','rna','rna',
                             'adt','adt','adt',
                             'rna','rna','rna',
                             'adt','adt','adt'),
                   batch = c('s1d1','s1d2','s1d3',
                             's1d1','s1d2','s1d3',
                             's2d1','s2d4','s2d5',
                             's2d1','s2d4','s2d5'),
                   q_or_batch = c(rep('ref',6),rep('query',6)),
                   modal_order = c('rna','adt'),
                   method = 'Harmony',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize','CLR'),
                   latent = latent_data,Bi_sPCA = T)

result_harmony <- run_Harmony(list(res[[1]],latent_data),
                              cell_name = c(colnames(res[[1]]),colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)
res <- Palette_map(c(ref,q),
                   modal = c('rna','rna','rna',
                             'adt','adt','adt',
                             'rna','rna','rna',
                             'adt','adt','adt'),
                   batch = c('s1d1','s1d2','s1d3',
                             's1d1','s1d2','s1d3',
                             's2d1','s2d4','s2d5',
                             's2d1','s2d4','s2d5'),
                   q_or_batch = c(rep('ref',6),rep('query',6)),
                   modal_order = c('rna','adt'),
                   method = 'fastMNN',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize','CLR'),
                   latent = latent_data,Bi_sPCA = T)

result_fastMNN <- run_fastMNN(list(res[[1]],latent_data),
                              cell_name = c(colnames(res[[1]]),colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)

latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised_30.csv',
                                 data.table = F)))
colnames(latent_data) <- rownames(meta)

res <- Palette_map(c(ref,q),
                   modal = c('rna','rna','rna',
                             'adt','adt','adt',
                             'rna','rna','rna',
                             'adt','adt','adt'),
                   batch = c('s1d1','s1d2','s1d3',
                             's1d1','s1d2','s1d3',
                             's2d1','s2d4','s2d5',
                             's2d1','s2d4','s2d5'),
                   q_or_batch = c(rep('ref',6),rep('query',6)),
                   modal_order = c('rna','adt'),
                   method = 'Seurat',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize','CLR'),
                   latent = latent_data,Bi_sPCA = T)

result_Seurat <- run_Seurat(list(res[[1]],latent_data),
                            cell_name = c(colnames(res[[1]]),colnames(latent_data)),
                            dims = 1:25,
                            is.normalize = F,out.npcs = NA)

meta_ref <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
meta_ref$RQ <- 'reference'
meta_q <- Reduce(rbind,meta_q.list)
meta_q$RQ <- 'query'
meta <- rbind(meta_ref,meta_q)

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

pred_label_l1_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_harmony <- evaluate(meta[rownames(meta_q),4],pred_label_l1_harmony[,1])
meta$pred_l1_harmony <- c(meta_ref[,4],pred_label_l1_harmony[,1])

pred_label_l2_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_harmony <- evaluate(meta[rownames(meta_q),5],pred_label_l2_harmony[,1])
meta$pred_l2_harmony <- c(meta_ref[,5],pred_label_l2_harmony[,1])

pred_label_l1_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_fastMNN <- evaluate(meta[rownames(meta_q),4],pred_label_l1_fastMNN[,1])
meta$pred_l1_fastMNN <- c(meta_ref[,4],pred_label_l1_fastMNN[,1])

pred_label_l2_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_fastMNN <- evaluate(meta[rownames(meta_q),5],pred_label_l2_fastMNN[,1])
meta$pred_l2_fastMNN <- c(meta_ref[,5],pred_label_l2_fastMNN[,1])

pred_label_l1_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),4]))
scores_l1_Seurat <- evaluate(meta[rownames(meta_q),4],pred_label_l1_Seurat[,1])
meta$pred_l1_Seurat <- c(meta_ref[,4],pred_label_l1_Seurat[,1])

pred_label_l2_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),5]))
scores_l2_Seurat <- evaluate(meta[rownames(meta_q),5],pred_label_l2_Seurat[,1])
meta$pred_l2_Seurat <- c(meta_ref[,5],pred_label_l2_Seurat[,1])

meta <- meta[c(rownames(meta_ref),rownames(meta_q)),]

score <- data.frame(Acc = c(scores_l1_harmony[[1]],scores_l2_harmony[[1]],
                            scores_l1_fastMNN[[1]],scores_l2_fastMNN[[1]],
                            scores_l1_Seurat[[1]],scores_l2_Seurat[[1]]),
                    Macro_F1 = c(scores_l1_harmony[[2]],scores_l2_harmony[[2]],
                                 scores_l1_fastMNN[[2]],scores_l2_fastMNN[[2]],
                                 scores_l1_Seurat[[2]],scores_l2_Seurat[[2]]))
rownames(score) <- c('l1_harmony','l2_harmony',
                     'l1_fastMNN','l2_fastMNN',
                     'l1_Seurat','l2_Seurat')


out_dir <- './Benchmarking/Reference-based_integration/output/Result/cite/'
write.csv(as.matrix(t(result_harmony[,rownames(meta)])),
          paste0(out_dir,'Palette_harmony_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_fastMNN[,rownames(meta)])),
          paste0(out_dir,'Palette_fastMNN_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_Seurat[,rownames(meta)])),
          paste0(out_dir,'Palette_Seurat_supervised.csv'),
          row.names = FALSE)
saveRDS(meta,paste0(out_dir,'meta_supervised.rds'))
saveRDS(meta_q,paste0(out_dir,'meta_query_supervised.rds'))
saveRDS(meta_ref,paste0(out_dir,'meta_ref_supervised.rds'))
saveRDS(score,paste0(out_dir,'score_supervised.rds'))

######################### CITE + Multiome ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised.csv',data.table = F)))
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
colnames(latent_data) <- rownames(meta)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite',
           './Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna.list <- list()
adt.list <- list()
atac.list <- list()
meta_ref.list <- list()

for(i in 1:6){
  rna.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  if(i %in% c(1,2,3)){
    adt.list[[i]] <-  readRDS(paste0(path1[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[i-3]] <-  readRDS(paste0(path1[i],'/atac.rds'))
  }
  meta_ref.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}

ref <- c(rna.list,adt.list,atac.list)

setwd('/home/server/sqy/data/mosaic/input/BMMC_CITE_Multiome')
path1 <- c('./CITE/s2d1_cite','./CITE/s2d4_cite','./CITE/s2d5_cite',
           './Multiome/s2d1_multi','./Multiome/s2d4_multi','./Multiome/s2d5_multi')

rna.list <- list()
adt.list <- list()
atac.list <- list()
meta_q.list <- list()
for(i in 1:6){
  rna.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  if(i %in% c(1,2,3)){
    adt.list[[i]] <-  readRDS(paste0(path1[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[i-3]] <-  readRDS(paste0(path1[i],'/atac.rds'))
  }
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
  # meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$modality,'-',meta_q.list[[i]]$batch)
  meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$batch,'-',meta_q.list[[i]]$modality)
}
q <- c(rna.list,adt.list,atac.list)
rm(rna.list,adt.list,atac.list)

res <- Palette_map(c(ref,q),
                   modal = c(rep('rna',6),rep('adt',3),rep('atac',3),
                             rep('rna',6),rep('adt',3),rep('atac',3)),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',12),rep('query',12)),
                   modal_order = c('rna','adt','atac'),
                   method = 'Harmony',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize','CLR','TF-IDF'),
                   latent = latent_data,Bi_sPCA = T)
result_harmony <- run_Harmony(list(res[[1]],
                                   res[[2]],
                                   latent_data),
                              cell_name = c(colnames(res[[1]]),
                                            colnames(res[[2]]),
                                            colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)

res <- Palette_map(c(ref,q),
                   modal = c(rep('rna',6),rep('adt',3),rep('atac',3),
                             rep('rna',6),rep('adt',3),rep('atac',3)),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',12),rep('query',12)),
                   modal_order = c('rna','adt','atac'),
                   method = 'fastMNN',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize','CLR','TF-IDF'),
                   latent = latent_data,Bi_sPCA = T)

result_fastMNN <- run_fastMNN(list(res[[1]],
                                   res[[2]],
                                   latent_data),
                              cell_name = c(colnames(res[[1]]),
                                            colnames(res[[2]]),
                                            colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)


latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised_30.csv',
                                 data.table = F)))
colnames(latent_data) <- rownames(meta)

res <- Palette_map(c(ref,q),
                   modal = c(rep('rna',6),rep('adt',3),rep('atac',3),
                             rep('rna',6),rep('adt',3),rep('atac',3)),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',12),rep('query',12)),
                   modal_order = c('rna','adt','atac'),
                   method = 'Seurat',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize','CLR','TF-IDF'),
                   latent = latent_data,Bi_sPCA = T)

result_Seurat <- run_Seurat(list(res[[1]],
                                 res[[2]],
                                 latent_data),
                            cell_name = c(colnames(res[[1]]),
                                          colnames(res[[2]]),
                                          colnames(latent_data)),
                            dims = 1:25,
                            is.normalize = F,out.npcs = NA)

meta_ref <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
meta_ref$RQ <- 'reference'
meta_q <- Reduce(rbind,meta_q.list)
meta_q$RQ <- 'query'
meta <- rbind(meta_ref,meta_q)
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

pred_label_l1_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_harmony <- evaluate(meta[rownames(meta_q),4],pred_label_l1_harmony[,1])
meta$pred_l1_harmony <- c(meta_ref[,4],pred_label_l1_harmony[,1])

pred_label_l2_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_harmony <- evaluate(meta[rownames(meta_q),5],pred_label_l2_harmony[,1])
meta$pred_l2_harmony <- c(meta_ref[,5],pred_label_l2_harmony[,1])

pred_label_l1_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_fastMNN <- evaluate(meta[rownames(meta_q),4],pred_label_l1_fastMNN[,1])
meta$pred_l1_fastMNN <- c(meta_ref[,4],pred_label_l1_fastMNN[,1])

pred_label_l2_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_fastMNN <- evaluate(meta[rownames(meta_q),5],pred_label_l2_fastMNN[,1])
meta$pred_l2_fastMNN <- c(meta_ref[,5],pred_label_l2_fastMNN[,1])

pred_label_l1_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),4]))
scores_l1_Seurat <- evaluate(meta[rownames(meta_q),4],pred_label_l1_Seurat[,1])
meta$pred_l1_Seurat <- c(meta_ref[,4],pred_label_l1_Seurat[,1])

pred_label_l2_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),5]))
scores_l2_Seurat <- evaluate(meta[rownames(meta_q),5],pred_label_l2_Seurat[,1])
meta$pred_l2_Seurat <- c(meta_ref[,5],pred_label_l2_Seurat[,1])

meta <- meta[c(rownames(meta_ref),rownames(meta_q)),]

score <- data.frame(Acc = c(scores_l1_harmony[[1]],scores_l2_harmony[[1]],
                            scores_l1_fastMNN[[1]],scores_l2_fastMNN[[1]],
                            scores_l1_Seurat[[1]],scores_l2_Seurat[[1]]),
                    Macro_F1 = c(scores_l1_harmony[[2]],scores_l2_harmony[[2]],
                                 scores_l1_fastMNN[[2]],scores_l2_fastMNN[[2]],
                                 scores_l1_Seurat[[2]],scores_l2_Seurat[[2]]))
rownames(score) <- c('l1_harmony','l2_harmony',
                     'l1_fastMNN','l2_fastMNN',
                     'l1_Seurat','l2_Seurat')

out_dir <- './Benchmarking/Reference-based_integration/output/Result/cite_multi/'
write.csv(as.matrix(t(result_harmony[,rownames(meta)])),
          paste0(out_dir,'Palette_harmony_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_fastMNN[,rownames(meta)])),
          paste0(out_dir,'Palette_fastMNN_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_Seurat[,rownames(meta)])),
          paste0(out_dir,'Palette_Seurat_supervised.csv'),
          row.names = FALSE)
saveRDS(meta,paste0(out_dir,'meta_supervised.rds'))
saveRDS(meta_q,paste0(out_dir,'meta_query_supervised.rds'))
saveRDS(meta_ref,paste0(out_dir,'meta_ref_supervised.rds'))
saveRDS(score,paste0(out_dir,'score_supervised.rds'))


######################### Multiome ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised.csv',data.table = F)))
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
colnames(latent_data) <- rownames(meta)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna1.list <- list()
atac.list <- list()
meta_ref.list <- list()
for(i in 1:3){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path1[i],'/atac.rds'))
  meta_ref.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}
ref <- c(rna1.list,atac.list)
rm(rna1.list,atac.list)

setwd('/home/server/sqy/data/mosaic/input/BMMC_CITE_Multiome')
path1 <- c('./Multiome/s2d1_multi','./Multiome/s2d4_multi','./Multiome/s2d5_multi')
rna1.list <- list()
atac.list <- list()
meta_q.list <- list()
for(i in 1:3){
  rna1.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  atac.list[[i]] <- readRDS(paste0(path1[i],'/atac.rds'))
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}
q <- c(rna1.list,atac.list)
rm(rna1.list,atac.list)
gc()

res <- Palette_map(c(ref,q),
                   modal = c('rna','rna','rna',
                             'atac','atac','atac',
                             'rna','rna','rna',
                             'atac','atac','atac'),
                   batch = c('s1d1','s1d2','s1d3',
                             's1d1','s1d2','s1d3',
                             's2d1','s2d4','s2d5',
                             's2d1','s2d4','s2d5'),
                   q_or_batch = c(rep('ref',6),rep('query',6)),
                   modal_order = c('rna','atac'),
                   method = 'Harmony',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize','TF-IDF'),
                   latent = latent_data,Bi_sPCA = T)
result_harmony <- run_Harmony(list(res[[1]],latent_data),
                              cell_name = c(colnames(res[[1]]),colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)

res <- Palette_map(c(ref,q),
                   modal = c('rna','rna','rna',
                             'atac','atac','atac',
                             'rna','rna','rna',
                             'atac','atac','atac'),
                   batch = c('s1d1','s1d2','s1d3',
                             's1d1','s1d2','s1d3',
                             's2d1','s2d4','s2d5',
                             's2d1','s2d4','s2d5'),
                   q_or_batch = c(rep('ref',6),rep('query',6)),
                   modal_order = c('rna','atac'),
                   method = 'fastMNN',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize','TF-IDF'),
                   latent = latent_data,Bi_sPCA = T)
result_fastMNN <- run_fastMNN(list(res[[1]],latent_data),
                              cell_name = c(colnames(res[[1]]),colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)

latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised_30.csv',
                                 data.table = F)))
colnames(latent_data) <- rownames(meta)

res <- Palette_map(c(ref,q),
                   modal = c('rna','rna','rna',
                             'atac','atac','atac',
                             'rna','rna','rna',
                             'atac','atac','atac'),
                   batch = c('s1d1','s1d2','s1d3',
                             's1d1','s1d2','s1d3',
                             's2d1','s2d4','s2d5',
                             's2d1','s2d4','s2d5'),
                   q_or_batch = c(rep('ref',6),rep('query',6)),
                   modal_order = c('rna','atac'),
                   method = 'Seurat',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize','TF-IDF'),
                   latent = latent_data,Bi_sPCA = T)

result_Seurat <- run_Seurat(list(res[[1]],latent_data),
                            cell_name = c(colnames(res[[1]]),colnames(latent_data)),
                            dims = 1:25,
                            is.normalize = F,out.npcs = NA)

meta_ref <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
meta_ref$RQ <- 'reference'
meta_q <- Reduce(rbind,meta_q.list)
meta_q$RQ <- 'query'
meta <- rbind(meta_ref,meta_q)

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

pred_label_l1_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_harmony <- evaluate(meta[rownames(meta_q),4],pred_label_l1_harmony[,1])
meta$pred_l1_harmony <- c(meta_ref[,4],pred_label_l1_harmony[,1])

pred_label_l2_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_harmony <- evaluate(meta[rownames(meta_q),5],pred_label_l2_harmony[,1])
meta$pred_l2_harmony <- c(meta_ref[,5],pred_label_l2_harmony[,1])

pred_label_l1_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_fastMNN <- evaluate(meta[rownames(meta_q),4],pred_label_l1_fastMNN[,1])
meta$pred_l1_fastMNN <- c(meta_ref[,4],pred_label_l1_fastMNN[,1])

pred_label_l2_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_fastMNN <- evaluate(meta[rownames(meta_q),5],pred_label_l2_fastMNN[,1])
meta$pred_l2_fastMNN <- c(meta_ref[,5],pred_label_l2_fastMNN[,1])

pred_label_l1_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),4]))
scores_l1_Seurat <- evaluate(meta[rownames(meta_q),4],pred_label_l1_Seurat[,1])
meta$pred_l1_Seurat <- c(meta_ref[,4],pred_label_l1_Seurat[,1])

pred_label_l2_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),5]))
scores_l2_Seurat <- evaluate(meta[rownames(meta_q),5],pred_label_l2_Seurat[,1])
meta$pred_l2_Seurat <- c(meta_ref[,5],pred_label_l2_Seurat[,1])

meta <- meta[c(rownames(meta_ref),rownames(meta_q)),]

score <- data.frame(Acc = c(scores_l1_harmony[[1]],scores_l2_harmony[[1]],
                            scores_l1_fastMNN[[1]],scores_l2_fastMNN[[1]],
                            scores_l1_Seurat[[1]],scores_l2_Seurat[[1]]),
                    Macro_F1 = c(scores_l1_harmony[[2]],scores_l2_harmony[[2]],
                                 scores_l1_fastMNN[[2]],scores_l2_fastMNN[[2]],
                                 scores_l1_Seurat[[2]],scores_l2_Seurat[[2]]))
rownames(score) <- c('l1_harmony','l2_harmony',
                     'l1_fastMNN','l2_fastMNN',
                     'l1_Seurat','l2_Seurat')

out_dir <- './Benchmarking/Reference-based_integration/output/Result/multiome/'
write.csv(as.matrix(t(result_harmony[,rownames(meta)])),
          paste0(out_dir,'Palette_harmony_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_fastMNN[,rownames(meta)])),
          paste0(out_dir,'Palette_fastMNN_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_Seurat[,rownames(meta)])),
          paste0(out_dir,'Palette_Seurat_supervised.csv'),
          row.names = FALSE)
saveRDS(meta,paste0(out_dir,'meta_supervised.rds'))
saveRDS(meta_q,paste0(out_dir,'meta_query_supervised.rds'))
saveRDS(meta_ref,paste0(out_dir,'meta_ref_supervised.rds'))
saveRDS(score,paste0(out_dir,'score_supervised.rds'))

######################### RNA ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised.csv',data.table = F)))
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
colnames(latent_data) <- rownames(meta)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite',
           './Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
ref <- list()
meta_ref.list <- list()
for(i in 1:6){
  ref[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  meta_ref.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}

setwd('/home/server/sqy/data/mosaic/input/BMMC_CITE_Multiome')
path1 <- c('./CITE/s2d1_cite','./CITE/s2d4_cite','./CITE/s2d5_cite',
           './Multiome/s2d1_multi','./Multiome/s2d4_multi','./Multiome/s2d5_multi')
q <- list()
meta_q.list <- list()
for(i in 1:6){
  q[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
  # meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$modality,'-',meta_q.list[[i]]$batch)
  meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$batch,'-',meta_q.list[[i]]$modality)
  meta_q.list[[i]]$modality <- 'RNA'
}
res <- Palette_map(c(ref,q),
                   modal = c('rna','rna','rna',
                             'rna','rna','rna',
                             'rna','rna','rna',
                             'rna','rna','rna'),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',6),rep('query',6)),
                   modal_order = c('rna'),
                   method = 'Harmony',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize'),
                   latent = latent_data,Bi_sPCA = T)

result_harmony <- run_Harmony(list(res[[1]],latent_data),
                              cell_name = c(colnames(res[[1]]),colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)

res <- Palette_map(c(ref,q),
                   modal = c('rna','rna','rna',
                             'rna','rna','rna',
                             'rna','rna','rna',
                             'rna','rna','rna'),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',6),rep('query',6)),
                   modal_order = c('rna'),
                   method = 'fastMNN',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize'),
                   latent = latent_data,Bi_sPCA = T)

result_fastMNN <- run_fastMNN(list(res[[1]],latent_data),
                              cell_name = c(colnames(res[[1]]),colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)

latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised_30.csv',
                                 data.table = F)))
colnames(latent_data) <- rownames(meta)
res <- Palette_map(c(ref,q),
                   modal = c('rna','rna','rna',
                             'rna','rna','rna',
                             'rna','rna','rna',
                             'rna','rna','rna'),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',6),rep('query',6)),
                   modal_order = c('rna'),
                   method = 'Seurat',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize'),
                   latent = latent_data,Bi_sPCA = T)

result_Seurat <- run_Seurat(list(res[[1]],latent_data),
                            cell_name = c(colnames(res[[1]]),colnames(latent_data)),
                            dims = 1:25,
                            is.normalize = F,out.npcs = NA)

meta_ref <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
meta_ref$RQ <- 'reference'
meta_q <- Reduce(rbind,meta_q.list)
meta_q$RQ <- 'query'
meta <- rbind(meta_ref,meta_q)

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

pred_label_l1_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_harmony <- evaluate(meta[rownames(meta_q),4],pred_label_l1_harmony[,1])
meta$pred_l1_harmony <- c(meta_ref[,4],pred_label_l1_harmony[,1])

pred_label_l2_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_harmony <- evaluate(meta[rownames(meta_q),5],pred_label_l2_harmony[,1])
meta$pred_l2_harmony <- c(meta_ref[,5],pred_label_l2_harmony[,1])

pred_label_l1_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_fastMNN <- evaluate(meta[rownames(meta_q),4],pred_label_l1_fastMNN[,1])
meta$pred_l1_fastMNN <- c(meta_ref[,4],pred_label_l1_fastMNN[,1])

pred_label_l2_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_fastMNN <- evaluate(meta[rownames(meta_q),5],pred_label_l2_fastMNN[,1])
meta$pred_l2_fastMNN <- c(meta_ref[,5],pred_label_l2_fastMNN[,1])

pred_label_l1_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),4]))
scores_l1_Seurat <- evaluate(meta[rownames(meta_q),4],pred_label_l1_Seurat[,1])
meta$pred_l1_Seurat <- c(meta_ref[,4],pred_label_l1_Seurat[,1])

pred_label_l2_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),5]))
scores_l2_Seurat <- evaluate(meta[rownames(meta_q),5],pred_label_l2_Seurat[,1])
meta$pred_l2_Seurat <- c(meta_ref[,5],pred_label_l2_Seurat[,1])

meta <- meta[c(rownames(meta_ref),rownames(meta_q)),]

score <- data.frame(Acc = c(scores_l1_harmony[[1]],scores_l2_harmony[[1]],
                            scores_l1_fastMNN[[1]],scores_l2_fastMNN[[1]],
                            scores_l1_Seurat[[1]],scores_l2_Seurat[[1]]),
                    Macro_F1 = c(scores_l1_harmony[[2]],scores_l2_harmony[[2]],
                                 scores_l1_fastMNN[[2]],scores_l2_fastMNN[[2]],
                                 scores_l1_Seurat[[2]],scores_l2_Seurat[[2]]))
rownames(score) <- c('l1_harmony','l2_harmony',
                     'l1_fastMNN','l2_fastMNN',
                     'l1_Seurat','l2_Seurat')
out_dir <- './Benchmarking/Reference-based_integration/output/Result/rna/'
write.csv(as.matrix(t(result_harmony[,rownames(meta)])),
          paste0(out_dir,'Palette_harmony_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_fastMNN[,rownames(meta)])),
          paste0(out_dir,'Palette_fastMNN_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_Seurat[,rownames(meta)])),
          paste0(out_dir,'Palette_Seurat_supervised.csv'),
          row.names = FALSE)
saveRDS(meta,paste0(out_dir,'meta_supervised.rds'))
saveRDS(meta_q,paste0(out_dir,'meta_query_supervised.rds'))
saveRDS(meta_ref,paste0(out_dir,'meta_ref_supervised.rds'))
saveRDS(score,paste0(out_dir,'score_supervised.rds'))

######################### RNA + ADT ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised.csv',data.table = F)))
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
colnames(latent_data) <- rownames(meta)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite',
           './Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna.list <- list()
adt.list <- list()
# atac.list <- list()
meta_ref.list <- list()

for(i in 1:6){
  rna.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  if(i %in% c(1,2,3)){
    adt.list[[i]] <-  readRDS(paste0(path1[i],'/adt.rds'))
  }
  meta_ref.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}

ref <- c(rna.list,adt.list)

setwd('/home/server/sqy/data/mosaic/input/BMMC_CITE_Multiome')
path1 <- c('./CITE/s2d1_cite','./CITE/s2d4_cite','./CITE/s2d5_cite',
           './Multiome/s2d1_multi','./Multiome/s2d4_multi','./Multiome/s2d5_multi')

rna.list <- list()
adt.list <- list()
# atac.list <- list()
meta_q.list <- list()
for(i in 1:6){
  if(i %in% c(1,2,3)){
    adt.list[[i]] <-  readRDS(paste0(path1[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    rna.list[[i-3]] <-  readRDS(paste0(path1[i],'/rna.rds'))
  }
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
  # meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$modality,'-',meta_q.list[[i]]$batch)
  meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$batch,'-',meta_q.list[[i]]$modality)
  if(i %in% c(1,2,3)){
    meta_q.list[[i]]$modality <- 'ADT'
  }else{
    meta_q.list[[i]]$modality <- 'RNA'
  }
}
q <- c(adt.list,rna.list)
rm(rna.list,adt.list)

res <- Palette_map(c(ref,q),
                   modal = c(rep('rna',6),rep('adt',3),
                             rep('adt',3),rep('rna',3)),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's1d1_cite','s1d2_cite','s1d3_cite',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',9),rep('query',6)),
                   modal_order = c('rna','adt'),
                   method = 'Harmony',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize','CLR'),
                   latent = latent_data,Bi_sPCA = T)

result_harmony <- run_Harmony(list(res[[1]],
                                   res[[2]],
                                   latent_data),
                              cell_name = c(colnames(res[[1]]),
                                            colnames(res[[2]]),
                                            colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)

res <- Palette_map(c(ref,q),
                   modal = c(rep('rna',6),rep('adt',3),
                             rep('adt',3),rep('rna',3)),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's1d1_cite','s1d2_cite','s1d3_cite',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',9),rep('query',6)),
                   modal_order = c('rna','adt'),
                   method = 'fastMNN',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize','CLR'),
                   latent = latent_data,Bi_sPCA = T)

result_fastMNN <- run_fastMNN(list(res[[1]],
                                   res[[2]],
                                   latent_data),
                              cell_name = c(colnames(res[[1]]),
                                            colnames(res[[2]]),
                                            colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)

latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised_30.csv',
                                 data.table = F)))
colnames(latent_data) <- rownames(meta)

res <- Palette_map(c(ref,q),
                   modal = c(rep('rna',6),rep('adt',3),
                             rep('adt',3),rep('rna',3)),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's1d1_cite','s1d2_cite','s1d3_cite',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',9),rep('query',6)),
                   modal_order = c('rna','adt'),
                   method = 'Seurat',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize','CLR'),
                   latent = latent_data,Bi_sPCA = T)

result_Seurat <- run_Seurat(list(res[[1]],
                                 res[[2]],
                                 latent_data),
                            cell_name = c(colnames(res[[1]]),
                                          colnames(res[[2]]),
                                          colnames(latent_data)),
                            dims = 1:25,
                            is.normalize = F,out.npcs = NA)

meta_ref <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
meta_ref$RQ <- 'reference'
meta_q <- Reduce(rbind,meta_q.list)
meta_q$RQ <- 'query'
meta <- rbind(meta_ref,meta_q)

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

pred_label_l1_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_harmony <- evaluate(meta[rownames(meta_q),4],pred_label_l1_harmony[,1])
meta$pred_l1_harmony <- c(meta_ref[,4],pred_label_l1_harmony[,1])

pred_label_l2_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_harmony <- evaluate(meta[rownames(meta_q),5],pred_label_l2_harmony[,1])
meta$pred_l2_harmony <- c(meta_ref[,5],pred_label_l2_harmony[,1])

pred_label_l1_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_fastMNN <- evaluate(meta[rownames(meta_q),4],pred_label_l1_fastMNN[,1])
meta$pred_l1_fastMNN <- c(meta_ref[,4],pred_label_l1_fastMNN[,1])

pred_label_l2_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_fastMNN <- evaluate(meta[rownames(meta_q),5],pred_label_l2_fastMNN[,1])
meta$pred_l2_fastMNN <- c(meta_ref[,5],pred_label_l2_fastMNN[,1])

pred_label_l1_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),4]))
scores_l1_Seurat <- evaluate(meta[rownames(meta_q),4],pred_label_l1_Seurat[,1])
meta$pred_l1_Seurat <- c(meta_ref[,4],pred_label_l1_Seurat[,1])

pred_label_l2_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),5]))
scores_l2_Seurat <- evaluate(meta[rownames(meta_q),5],pred_label_l2_Seurat[,1])
meta$pred_l2_Seurat <- c(meta_ref[,5],pred_label_l2_Seurat[,1])

meta <- meta[c(rownames(meta_ref),rownames(meta_q)),]

score <- data.frame(Acc = c(scores_l1_harmony[[1]],scores_l2_harmony[[1]],
                            scores_l1_fastMNN[[1]],scores_l2_fastMNN[[1]],
                            scores_l1_Seurat[[1]],scores_l2_Seurat[[1]]),
                    Macro_F1 = c(scores_l1_harmony[[2]],scores_l2_harmony[[2]],
                                 scores_l1_fastMNN[[2]],scores_l2_fastMNN[[2]],
                                 scores_l1_Seurat[[2]],scores_l2_Seurat[[2]]))
rownames(score) <- c('l1_harmony','l2_harmony',
                     'l1_fastMNN','l2_fastMNN',
                     'l1_Seurat','l2_Seurat')

out_dir <- './Benchmarking/Reference-based_integration/output/Result/rna_adt/'
write.csv(as.matrix(t(result_harmony[,rownames(meta)])),
          paste0(out_dir,'Palette_harmony_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_fastMNN[,rownames(meta)])),
          paste0(out_dir,'Palette_fastMNN_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_Seurat[,rownames(meta)])),
          paste0(out_dir,'Palette_Seurat_supervised.csv'),
          row.names = FALSE)
saveRDS(meta,paste0(out_dir,'meta_supervised.rds'))
saveRDS(meta_q,paste0(out_dir,'meta_query_supervised.rds'))
saveRDS(meta_ref,paste0(out_dir,'meta_ref_supervised.rds'))
saveRDS(score,paste0(out_dir,'score_supervised.rds'))


######################### RNA + ADT + ATAC ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised.csv',data.table = F)))
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
colnames(latent_data) <- rownames(meta)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite',
           './Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna.list <- list()
adt.list <- list()
atac.list <- list()
meta_ref.list <- list()

for(i in 1:6){
  rna.list[[i]] <-  readRDS(paste0(path1[i],'/rna.rds'))
  if(i %in% c(1,2,3)){
    adt.list[[i]] <-  readRDS(paste0(path1[i],'/adt.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[i-3]] <-  readRDS(paste0(path1[i],'/atac.rds'))
  }
  meta_ref.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}

ref <- c(rna.list,adt.list,atac.list)

setwd('/home/server/sqy/data/mosaic/input/BMMC_CITE_Multiome')
path1 <- c('./CITE/s2d1_cite','./CITE/s2d4_cite','./CITE/s2d5_cite',
           './Multiome/s2d1_multi','./Multiome/s2d4_multi','./Multiome/s2d5_multi')

rna.list <- list()
adt.list <- list()
atac.list <- list()
meta_q.list <- list()
for(i in 1:6){
  if(i %in% c(1,2)){
    adt.list[[i]] <-  readRDS(paste0(path1[i],'/adt.rds'))
  }
  if(i %in% c(3,4)){
    rna.list[[i-2]] <-  readRDS(paste0(path1[i],'/rna.rds'))
  }
  if(i %in% c(5,6)){
    atac.list[[i-4]] <-  readRDS(paste0(path1[i],'/atac.rds'))
  }
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
  # meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$modality,'-',meta_q.list[[i]]$batch)
  meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$batch,'-',meta_q.list[[i]]$modality)
  if(i %in% c(1,2)){
    meta_q.list[[i]]$modality <- 'ADT'
  }
  if(i %in% c(3,4)){
    meta_q.list[[i]]$modality <- 'RNA'
  }
  if(i %in% c(5,6)){
    meta_q.list[[i]]$modality <- 'ATAC'
  }
  
}
q <- c(adt.list,rna.list,atac.list)
rm(adt.list,atac.list,rna.list)

res <- Palette_map(c(ref,q),
                   modal = c(rep('rna',6),rep('adt',3),rep('atac',3),
                             rep('adt',2),rep('rna',2),rep('atac',2)),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',12),rep('query',6)),
                   modal_order = c('rna','adt','atac'),
                   method = 'Harmony',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize','CLR','TF-IDF'),
                   latent = latent_data,Bi_sPCA = T)

result_harmony <- run_Harmony(list(res[[1]],
                                   res[[2]],
                                   res[[3]],
                                   latent_data),
                              cell_name = c(colnames(res[[1]]),
                                            colnames(res[[2]]),
                                            colnames(res[[3]]),
                                            colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)

res <- Palette_map(c(ref,q),
                   modal = c(rep('rna',6),rep('adt',3),rep('atac',3),
                             rep('adt',2),rep('rna',2),rep('atac',2)),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',12),rep('query',6)),
                   modal_order = c('rna','adt','atac'),
                   method = 'fastMNN',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize','CLR','TF-IDF'),
                   latent = latent_data,Bi_sPCA = T)

result_fastMNN <- run_fastMNN(list(res[[1]],
                                   res[[2]],
                                   res[[3]],
                                   latent_data),
                              cell_name = c(colnames(res[[1]]),
                                            colnames(res[[2]]),
                                            colnames(res[[3]]),
                                            colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)

latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised_30.csv',
                                 data.table = F)))
colnames(latent_data) <- rownames(meta)

res <- Palette_map(c(ref,q),
                   modal = c(rep('rna',6),rep('adt',3),rep('atac',3),
                             rep('adt',2),rep('rna',2),rep('atac',2)),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',12),rep('query',6)),
                   modal_order = c('rna','adt','atac'),
                   method = 'Seurat',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize','CLR','TF-IDF'),
                   latent = latent_data,Bi_sPCA = T)

result_Seurat <- run_Seurat(list(res[[1]],
                                 res[[2]],
                                 res[[3]],
                                 latent_data),
                            cell_name = c(colnames(res[[1]]),
                                          colnames(res[[2]]),
                                          colnames(res[[3]]),
                                          colnames(latent_data)),
                            dims = 1:25,
                            is.normalize = F,out.npcs = NA)

meta_ref <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
meta_ref$RQ <- 'reference'
meta_q <- Reduce(rbind,meta_q.list)
meta_q$RQ <- 'query'
meta <- rbind(meta_ref,meta_q)
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

pred_label_l1_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_harmony <- evaluate(meta[rownames(meta_q),4],pred_label_l1_harmony[,1])
meta$pred_l1_harmony <- c(meta_ref[,4],pred_label_l1_harmony[,1])

pred_label_l2_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_harmony <- evaluate(meta[rownames(meta_q),5],pred_label_l2_harmony[,1])
meta$pred_l2_harmony <- c(meta_ref[,5],pred_label_l2_harmony[,1])

pred_label_l1_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_fastMNN <- evaluate(meta[rownames(meta_q),4],pred_label_l1_fastMNN[,1])
meta$pred_l1_fastMNN <- c(meta_ref[,4],pred_label_l1_fastMNN[,1])

pred_label_l2_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_fastMNN <- evaluate(meta[rownames(meta_q),5],pred_label_l2_fastMNN[,1])
meta$pred_l2_fastMNN <- c(meta_ref[,5],pred_label_l2_fastMNN[,1])

pred_label_l1_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),4]))
scores_l1_Seurat <- evaluate(meta[rownames(meta_q),4],pred_label_l1_Seurat[,1])
meta$pred_l1_Seurat <- c(meta_ref[,4],pred_label_l1_Seurat[,1])

pred_label_l2_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),5]))
scores_l2_Seurat <- evaluate(meta[rownames(meta_q),5],pred_label_l2_Seurat[,1])
meta$pred_l2_Seurat <- c(meta_ref[,5],pred_label_l2_Seurat[,1])

meta <- meta[c(rownames(meta_ref),rownames(meta_q)),]

score <- data.frame(Acc = c(scores_l1_harmony[[1]],scores_l2_harmony[[1]],
                            scores_l1_fastMNN[[1]],scores_l2_fastMNN[[1]],
                            scores_l1_Seurat[[1]],scores_l2_Seurat[[1]]),
                    Macro_F1 = c(scores_l1_harmony[[2]],scores_l2_harmony[[2]],
                                 scores_l1_fastMNN[[2]],scores_l2_fastMNN[[2]],
                                 scores_l1_Seurat[[2]],scores_l2_Seurat[[2]]))
rownames(score) <- c('l1_harmony','l2_harmony',
                     'l1_fastMNN','l2_fastMNN',
                     'l1_Seurat','l2_Seurat')

out_dir <- './Benchmarking/Reference-based_integration/output/Result/rna_adt_atac/'
write.csv(as.matrix(t(result_harmony[,rownames(meta)])),
          paste0(out_dir,'Palette_harmony_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_fastMNN[,rownames(meta)])),
          paste0(out_dir,'Palette_fastMNN_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_Seurat[,rownames(meta)])),
          paste0(out_dir,'Palette_Seurat_supervised.csv'),
          row.names = FALSE)
saveRDS(meta,paste0(out_dir,'meta_supervised.rds'))
saveRDS(meta_q,paste0(out_dir,'meta_query_supervised.rds'))
saveRDS(meta_ref,paste0(out_dir,'meta_ref_supervised.rds'))
saveRDS(score,paste0(out_dir,'score_supervised.rds'))


######################### RNA + ATAC ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
library(data.table)
options(future.globals.maxSize = 5000 * 1024^3)
latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised.csv',data.table = F)))
meta <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
colnames(latent_data) <- rownames(meta)
setwd('./BMMC_CITE_Multiome')
path1 <- c('./CITE/s1d1_cite','./CITE/s1d2_cite','./CITE/s1d3_cite',
           './Multiome/s1d1_multi','./Multiome/s1d2_multi','./Multiome/s1d3_multi')
rna.list <- list()
# adt.list <- list()
atac.list <- list()
meta_ref.list <- list()

for(i in 1:6){
  rna.list[[i]] <- readRDS(paste0(path1[i],'/rna.rds'))
  if(i %in% c(4,5,6)){
    atac.list[[i-3]] <-  readRDS(paste0(path1[i],'/atac.rds'))
  }
  meta_ref.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
}

ref <- c(rna.list,atac.list)

setwd('/home/server/sqy/data/mosaic/input/BMMC_CITE_Multiome')
path1 <- c('./CITE/s2d1_cite','./CITE/s2d4_cite','./CITE/s2d5_cite',
           './Multiome/s2d1_multi','./Multiome/s2d4_multi','./Multiome/s2d5_multi')

rna.list <- list()
# adt.list <- list()
atac.list <- list()
meta_q.list <- list()
for(i in 1:6){
  if(i %in% c(1,2,3)){
    rna.list[[i]] <-  readRDS(paste0(path1[i],'/rna.rds'))
  }
  if(i %in% c(4,5,6)){
    atac.list[[i-3]] <-  readRDS(paste0(path1[i],'/atac.rds'))
  }
  meta_q.list[[i]] <- readRDS(paste0(path1[i],'/meta.rds'))
  # meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$modality,'-',meta_q.list[[i]]$batch)
  meta_q.list[[i]]$batch <- paste0(meta_q.list[[i]]$batch,'-',meta_q.list[[i]]$modality)
  
  if(i %in% c(1,2,3)){
    meta_q.list[[i]]$modality <- 'RNA'
  }
  if(i %in% c(4,5,6)){
    meta_q.list[[i]]$modality <- 'ATAC'
  }
}
q <- c(rna.list,atac.list)
rm(rna.list,atac.list)

res <- Palette_map(c(ref,q),
                   modal = c(rep('rna',6),rep('atac',3),
                             rep('rna',3),rep('atac',3)),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',9),rep('query',6)),
                   modal_order = c('rna','atac'),
                   method = 'Harmony',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize','TF-IDF'),
                   latent = latent_data,Bi_sPCA = T)
result_harmony <- run_Harmony(list(res[[1]],
                                   res[[2]],
                                   latent_data),
                              cell_name = c(colnames(res[[1]]),
                                            colnames(res[[2]]),
                                            colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)

res <- Palette_map(c(ref,q),
                   modal = c(rep('rna',6),rep('atac',3),
                             rep('rna',3),rep('atac',3)),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',9),rep('query',6)),
                   modal_order = c('rna','atac'),
                   method = 'fastMNN',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize','TF-IDF'),
                   latent = latent_data,Bi_sPCA = T)

result_fastMNN <- run_fastMNN(list(res[[1]],
                                   res[[2]],
                                   latent_data),
                              cell_name = c(colnames(res[[1]]),
                                            colnames(res[[2]]),
                                            colnames(latent_data)),
                              is.normalize = F,out.npcs = NA)

latent_data <- t(as.matrix(fread('./Benchmarking/Reference-based_integration/output/Palette_Reference/Palette_supervised_30.csv',
                                 data.table = F)))
colnames(latent_data) <- rownames(meta)

res <- Palette_map(c(ref,q),
                   modal = c(rep('rna',6),rep('atac',3),
                             rep('rna',3),rep('atac',3)),
                   batch = c('s1d1_cite','s1d2_cite','s1d3_cite',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's1d1_multi','s1d2_multi','s1d3_multi',
                             's2d1_cite','s2d4_cite','s2d5_cite',
                             's2d1_multi','s2d4_multi','s2d5_multi'),
                   q_or_batch = c(rep('ref',9),rep('query',6)),
                   modal_order = c('rna','atac'),
                   method = 'Seurat',
                   method_out.npcs = 40,
                   normalize_method = c('LogNormalize','TF-IDF'),
                   latent = latent_data,Bi_sPCA = T)

result_Seurat <- run_Seurat(list(res[[1]],
                                 res[[2]],
                                 latent_data),
                            cell_name = c(colnames(res[[1]]),
                                          colnames(res[[2]]),
                                          colnames(latent_data)),
                            dims = 1:25,
                            is.normalize = F,out.npcs = NA)

meta_ref <- readRDS('./Benchmarking/Reference-based_integration/output/Palette_Reference/meta.rds')
meta_ref$RQ <- 'reference'
meta_q <- Reduce(rbind,meta_q.list)
meta_q$RQ <- 'query'
meta <- rbind(meta_ref,meta_q)

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
pred_label_l1_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_harmony <- evaluate(meta[rownames(meta_q),4],pred_label_l1_harmony[,1])
meta$pred_l1_harmony <- c(meta_ref[,4],pred_label_l1_harmony[,1])

pred_label_l2_harmony <- Trans.Knowledge(t(result_harmony[,rownames(meta_ref)]),
                                         t(result_harmony[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_harmony <- evaluate(meta[rownames(meta_q),5],pred_label_l2_harmony[,1])
meta$pred_l2_harmony <- c(meta_ref[,5],pred_label_l2_harmony[,1])

pred_label_l1_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),4]))
scores_l1_fastMNN <- evaluate(meta[rownames(meta_q),4],pred_label_l1_fastMNN[,1])
meta$pred_l1_fastMNN <- c(meta_ref[,4],pred_label_l1_fastMNN[,1])

pred_label_l2_fastMNN <- Trans.Knowledge(t(result_fastMNN[,rownames(meta_ref)]),
                                         t(result_fastMNN[,rownames(meta_q)]),
                                         data.frame(meta[rownames(meta_ref),5]))
scores_l2_fastMNN <- evaluate(meta[rownames(meta_q),5],pred_label_l2_fastMNN[,1])
meta$pred_l2_fastMNN <- c(meta_ref[,5],pred_label_l2_fastMNN[,1])

pred_label_l1_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),4]))
scores_l1_Seurat <- evaluate(meta[rownames(meta_q),4],pred_label_l1_Seurat[,1])
meta$pred_l1_Seurat <- c(meta_ref[,4],pred_label_l1_Seurat[,1])

pred_label_l2_Seurat <- Trans.Knowledge(t(result_Seurat[,rownames(meta_ref)]),
                                        t(result_Seurat[,rownames(meta_q)]),
                                        data.frame(meta[rownames(meta_ref),5]))
scores_l2_Seurat <- evaluate(meta[rownames(meta_q),5],pred_label_l2_Seurat[,1])
meta$pred_l2_Seurat <- c(meta_ref[,5],pred_label_l2_Seurat[,1])
meta <- meta[c(rownames(meta_ref),rownames(meta_q)),]
score <- data.frame(Acc = c(scores_l1_harmony[[1]],scores_l2_harmony[[1]],
                            scores_l1_fastMNN[[1]],scores_l2_fastMNN[[1]],
                            scores_l1_Seurat[[1]],scores_l2_Seurat[[1]]),
                    Macro_F1 = c(scores_l1_harmony[[2]],scores_l2_harmony[[2]],
                                 scores_l1_fastMNN[[2]],scores_l2_fastMNN[[2]],
                                 scores_l1_Seurat[[2]],scores_l2_Seurat[[2]]))
rownames(score) <- c('l1_harmony','l2_harmony',
                     'l1_fastMNN','l2_fastMNN',
                     'l1_Seurat','l2_Seurat')

out_dir <- './Benchmarking/Reference-based_integration/output/Result/rna_atac/'
write.csv(as.matrix(t(result_harmony[,rownames(meta)])),
          paste0(out_dir,'Palette_harmony_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_fastMNN[,rownames(meta)])),
          paste0(out_dir,'Palette_fastMNN_supervised.csv'),
          row.names = FALSE)
write.csv(as.matrix(t(result_Seurat[,rownames(meta)])),
          paste0(out_dir,'Palette_Seurat_supervised.csv'),
          row.names = FALSE)
saveRDS(meta,paste0(out_dir,'meta_supervised.rds'))
saveRDS(meta_q,paste0(out_dir,'meta_query_supervised.rds'))
saveRDS(meta_ref,paste0(out_dir,'meta_ref_supervised.rds'))
saveRDS(score,paste0(out_dir,'score_supervised.rds'))
