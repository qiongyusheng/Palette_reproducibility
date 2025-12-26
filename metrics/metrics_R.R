library(data.table)
library(Matrix)
library(matrixStats)
library(data.table)
library(Matrix)
library(matrixStats)
require(tidyverse)

clust_labels_to_int = function(clust_labels) {
  uniq_labels = unique(clust_labels)
  clust_labels_int = rep(NA, length(clust_labels))
  
  for (ii in 1:length(clust_labels)) {
    clust_labels_int[ii] = match(clust_labels[ii], uniq_labels)
  }
  
  return(clust_labels_int)
}


lisi = function(data = NULL,
                meta = NULL, 
                condition = "batch",
                label = 'celltype') {
  
  cLISI <- 0
  iLISI <- 0
  
  clust_label <- clust_labels_to_int(meta[,label])
  data_label <- clust_labels_to_int(meta[,condition])
  
  labels = as.data.frame(cbind(data_label, clust_label))
  colnames(labels) = c('data_label', 'clust_label')
  
  lisi_res <- lisi::compute_lisi(data, labels, c('clust_label','data_label'), perplexity = 40)
  Bi <- length(unique(data_label))
  Bc <-  length(unique(clust_label))
  cLISI_vec <- (Bc-lisi_res$clust_label)/(Bc-1)
  iLISI_vec <- (lisi_res$data_label-1)/(Bi-1)
  # iLISI_vec <- lisi_res$data_label
  # cLISI_vec <- lisi_res$clust_label
  iLISI <- median(iLISI_vec)
  cLISI <- median(cLISI_vec)
  
  LISI_f1_vec <- (2*iLISI_vec*cLISI_vec)/(iLISI_vec+cLISI_vec)
  
  # LISI_f1 <- (2*cLISI*iLISI)/(cLISI+iLISI)
  LISI_f1 <- median(LISI_f1_vec)
  return(cbind(iLISI, cLISI, LISI_f1))
  # return(list(cLISI_vec,iLISI_vec))
}

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

iLISI = function(data = NULL,
                 meta = NULL, 
                 condition = "batch") {
  
  iLISI <- 0
  
  data_label <- clust_labels_to_int(meta[,condition])
  
  labels = as.data.frame(data_label,)
  colnames(labels) = c('data_label')
  
  n_samples <- nrow(data)  
  safe_perplexity <- max(1, floor(min(40, (n_samples - 1) / 3)))
  lisi_res <- lisi::compute_lisi(data, labels, 'data_label', perplexity = safe_perplexity)
  Bi <- length(unique(data_label))
  iLISI_vec <- (lisi_res$data_label-1)/(Bi-1)
  # iLISI <- median(iLISI_vec)
  
  return(iLISI_vec)
}

CiLISI_F1 = function(data = NULL,
                     meta = NULL, 
                     condition = "batch",
                     label = 'celltype',
                     min.cell = 20){
  
  data.list <- split.data.frame(data,meta[,label])
  meta.list <- split.data.frame(meta,meta[,label])
  
  unique_label <- unique(meta[,label])
  
  idx <- sapply(unique_label, function(x) {
    m <- meta.list[[x]]
    celltype_count <- nrow(m)
    sample_count <- length(unique(m[[condition]])) 
    return(celltype_count >= min.cell &&
             sample_count >= 2)
  })
  
  unique_label <- unique_label[idx]
  
  data.list <- data.list[unique_label]
  meta.list <- meta.list[unique_label]
  
  iLISI_list <- lapply(seq_along(data.list), function(i) {
    x <- data.list[[i]]
    m <- meta.list[[i]]
    
    tmp_iLISI <- iLISI(data = x, meta = m, condition = condition)
    return(tmp_iLISI)
  })
  
  CiLISI_vec <- unlist(iLISI_list)
  CiLISI <- median(CiLISI_vec)
  # return(CiLISI)
  
  data <- Reduce(rbind,data.list)
  meta <- Reduce(rbind,meta.list)
  clust_label <- clust_labels_to_int(meta[,label])
  labels = as.data.frame(data_label,)
  colnames(labels) = c('data_label')
  
  lisi_res <- lisi::compute_lisi(data, labels, 'data_label', perplexity = 40)
  
  Bc <-  length(unique(clust_label))
  cLISI_vec <- (Bc-lisi_res$clust_label)/(Bc-1)
  cLISI <- median(cLISI_vec)
  
  LISI_f1_vec <- (2*CiLISI_vec*cLISI_vec)/(CiLISI_vec+cLISI_vec)
  LISI_f1 <- median(LISI_f1_vec)
  
  return(cbind(CiLISI, cLISI, LISI_f1))
}


check_data = function(data,meta,cut.off = 5,condition){
  data.list <- split.data.frame(data,meta[,condition])
  data_size <- unlist(lapply(1:length(data.list),function(x) nrow(data.list[[x]])))
  if(length(which(data_size > 5)) > 1){
    return(TRUE)
  }else{
    return(FALSE)
  }
  
}

sample_equal = function(data,meta,cut.off = 5,condition){
  
  data.list <- split.data.frame(data,meta[,condition])
  
  data_size <- unlist(lapply(1:length(data.list),function(x) nrow(data.list[[x]])))
  
  idx <- which(data_size > 5)
  data_size <- data_size[idx]
  data.list <- data.list[idx]
  meta.list <- split.data.frame(meta,meta[,condition])
  meta.list <- meta.list[idx]
  
  sample_size <- min(data_size)
  data_out <- list()
  meta_out <- list()
  for(i in 1:length(data.list)){
    n <- nrow(data.list[[i]])
    idx_tmp <- sample(1:n,size = sample_size,replace = FALSE)
    data_out[[i]] <- data.list[[i]][idx_tmp,]
    meta_out[[i]] <- meta.list[[i]][idx_tmp,]
  }
  
  return(list(data_out,
              meta_out))
}

CiLISI_sample = function(data = NULL,
                         meta = NULL, 
                         condition = "batch",
                         label = 'celltype',
                         min.cell = 20,
                         n_repeat = 10){
  
  data.list <- split.data.frame(data,meta[,label])
  meta.list <- split.data.frame(meta,meta[,label])
  
  unique_label <- unique(meta[,label])
  
  idx <- sapply(unique_label, function(x) {
    m <- meta.list[[x]]
    celltype_count <- nrow(m)
    sample_count <- length(unique(m[[condition]])) 
    return(celltype_count >= min.cell &&
             sample_count >= 2)
  })
  
  unique_label <- unique_label[idx]
  
  data.list <- data.list[unique_label]
  meta.list <- meta.list[unique_label]
  
  iLISI_list <- list()
  for(i in seq_along(data.list)){
    if(check_data(data = data.list[[i]],
                  meta = meta.list[[i]],
                  condition = condition)){
      
      iLISI_vec <- unlist(lapply(1:n_repeat,function(x){
        res <- sample_equal(data = data.list[[i]],
                            meta = meta.list[[i]],
                            cut.off = 5,
                            condition = condition)
        tmp_iLISI <- median(iLISI(data = Reduce(rbind,res[[1]]), 
                                  meta = Reduce(rbind,res[[2]]), 
                                  condition = condition))
        
        return(tmp_iLISI)
      }))
      
      iLISI_list[[i]] <- iLISI_vec
      
    }
  }
  
  CiLISI_vec <- unlist(iLISI_list)
  CiLISI <- median(CiLISI_vec)
  return(CiLISI)
}