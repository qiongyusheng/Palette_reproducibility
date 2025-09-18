library(data.table)
library(Matrix)
library(matrixStats)

data_dir <- '.../'
meta_dir <- '.../'
out_dir <- '.../'

method = c("Palette",
           'Multigrate',
           'scMoMaT',
           'midas',
           'scVAEIT',
           "StabMap")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0(data_dir,method[i],'.csv'),data.table = F,header = T)[,1:20]
}
meta <- readRDS(paste0(meta_dir,'meta.rds'))

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
  iLISI <- median(iLISI_vec)
  cLISI <- median(cLISI_vec)
  
  LISI_f1_vec <- (2*iLISI_vec*cLISI_vec)/(iLISI_vec+cLISI_vec)
  
  LISI_f1 <- median(LISI_f1_vec)
  return(cbind(iLISI, cLISI, LISI_f1))
}

out <- data.frame(F1_LISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  f1_lisi <- lisi(res[[i]],meta,'modality','celltype')
  out[method[i],'F1_LISI'] <- f1_lisi[1,3]
}

saveRDS(out,paste0(out_dir,'.../f1_lisi_score.rds'))


#################

# predict label evaluation
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

meta_modal1 <- dplyr::filter(meta,modality == "modal1")
meta_modal2 <- dplyr::filter(meta,modality == "modal2")

out <- data.frame(ACC = rep(0,length(method)),macro_F1 = rep(0,length(method)))
rownames(out) <- method

pred_label_df <- as.data.frame(matrix(0,nrow = nrow(meta_modal2),ncol = length(method)+1))
colnames(pred_label_df) <- c('ground_truth',method)
rownames(pred_label_df) <- rownames(meta_modal2)
pred_label_df$ground_truth <- meta_modal2$celltype

for(i in 1:length(method)){
  
  pred_tmp <- Trans.Knowledge(res[[i]][rownames(meta_modal1),],res[[i]][rownames(meta_modal2),],
                              data.frame(meta_modal1[,2]))
  
  pred_label_df[,method[i]] <- pred_tmp[,1]
  score_tmp <- evaluate(meta_modal2[,2],pred_tmp[,1])
  out[method[i],'ACC'] <- score_tmp[[1]]
  out[method[i],'macro_F1'] <- score_tmp[[2]]
}
saveRDS(out,'./ACC_macro_f1.rds')