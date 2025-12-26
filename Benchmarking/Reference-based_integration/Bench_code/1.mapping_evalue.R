######################### ADT ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Reference-based_integration/output/Result/adt')
meta <- readRDS('./meta.rds')
method <- c('MIDAS','Multigrate',
            'Palette_harmony','Palette_harmony_supervised',
            'Palette_fastMNN','Palette_fastMNN_supervised',
            'Palette_Seurat','Palette_Seurat_supervised',
            'Palette_denovo','Palette_denovo_supervised')

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}


out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_modality.rds')


############

out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'RQ','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_RQ.rds')

# Multigrate label transfer
mtg <- as.matrix(res[[2]])
rownames(mtg) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),4]))
scores_l1_mtg<- evaluate(meta[rownames(meta_q),4],pred_label_l1_mtg[,1])
meta$pred_l1_mtg <- c(meta_ref[,4],pred_label_l1_mtg[,1])

pred_label_l2_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),5]))
scores_l2_mtg<- evaluate(meta[rownames(meta_q),5],pred_label_l2_mtg[,1])
meta$pred_l2_mtg <- c(meta_ref[,5],pred_label_l2_mtg[,1])

score <- data.frame(Acc = c(scores_l1_mtg[[1]],scores_l2_mtg[[1]]),
                    Macro_F1 = c(scores_l1_mtg[[2]],scores_l2_mtg[[2]]))
rownames(score) <- c('l1_mtg','l2_mtg')
saveRDS(score,'./score_mtg.rds')

# MIDAS label transfer
MIDAS <- as.matrix(res[[1]])
rownames(MIDAS) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),4]))
scores_l1_midas<- evaluate(meta[rownames(meta_q),4],pred_label_l1_midas[,1])
meta$pred_l1_midas <- c(meta_ref[,4],pred_label_l1_midas[,1])

pred_label_l2_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),5]))
scores_l2_midas<- evaluate(meta[rownames(meta_q),5],pred_label_l2_midas[,1])
meta$pred_l2_midas <- c(meta_ref[,5],pred_label_l2_midas[,1])

score <- data.frame(Acc = c(scores_l1_midas[[1]],scores_l2_midas[[1]]),
                    Macro_F1 = c(scores_l1_midas[[2]],scores_l2_midas[[2]]))
rownames(score) <- c('l1_midas','l2_midas')
saveRDS(score,'./score_midas.rds')


######################### ADT + ATAC ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Reference-based_integration/output/Result/adt_atac')
meta <- readRDS('./meta.rds')
method <- c('MIDAS','Multigrate',
            'Palette_harmony','Palette_harmony_supervised',
            'Palette_fastMNN','Palette_fastMNN_supervised',
            'Palette_Seurat','Palette_Seurat_supervised',
            'Palette_denovo','Palette_denovo_supervised')

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}


out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_modality.rds')


############

out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'RQ','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_RQ.rds')

# Multigrate label transfer
mtg <- as.matrix(res[[2]])
rownames(mtg) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),4]))
scores_l1_mtg<- evaluate(meta[rownames(meta_q),4],pred_label_l1_mtg[,1])
meta$pred_l1_mtg <- c(meta_ref[,4],pred_label_l1_mtg[,1])

pred_label_l2_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),5]))
scores_l2_mtg<- evaluate(meta[rownames(meta_q),5],pred_label_l2_mtg[,1])
meta$pred_l2_mtg <- c(meta_ref[,5],pred_label_l2_mtg[,1])

score <- data.frame(Acc = c(scores_l1_mtg[[1]],scores_l2_mtg[[1]]),
                    Macro_F1 = c(scores_l1_mtg[[2]],scores_l2_mtg[[2]]))
rownames(score) <- c('l1_mtg','l2_mtg')
saveRDS(score,'./score_mtg.rds')

# MIDAS label transfer
MIDAS <- as.matrix(res[[1]])
rownames(MIDAS) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),4]))
scores_l1_midas<- evaluate(meta[rownames(meta_q),4],pred_label_l1_midas[,1])
meta$pred_l1_midas <- c(meta_ref[,4],pred_label_l1_midas[,1])

pred_label_l2_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),5]))
scores_l2_midas<- evaluate(meta[rownames(meta_q),5],pred_label_l2_midas[,1])
meta$pred_l2_midas <- c(meta_ref[,5],pred_label_l2_midas[,1])

score <- data.frame(Acc = c(scores_l1_midas[[1]],scores_l2_midas[[1]]),
                    Macro_F1 = c(scores_l1_midas[[2]],scores_l2_midas[[2]]))
rownames(score) <- c('l1_midas','l2_midas')
saveRDS(score,'./score_midas.rds')

######################### ATAC ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Reference-based_integration/output/Result/atac')
meta <- readRDS('./meta.rds')
method <- c('MIDAS','Multigrate',
            'Palette_harmony','Palette_harmony_supervised',
            'Palette_fastMNN','Palette_fastMNN_supervised',
            'Palette_Seurat','Palette_Seurat_supervised',
            'Palette_denovo','Palette_denovo_supervised')

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}


out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_modality.rds')


############

out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'RQ','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_RQ.rds')

# Multigrate label transfer
mtg <- as.matrix(res[[2]])
rownames(mtg) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),4]))
scores_l1_mtg<- evaluate(meta[rownames(meta_q),4],pred_label_l1_mtg[,1])
meta$pred_l1_mtg <- c(meta_ref[,4],pred_label_l1_mtg[,1])

pred_label_l2_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),5]))
scores_l2_mtg<- evaluate(meta[rownames(meta_q),5],pred_label_l2_mtg[,1])
meta$pred_l2_mtg <- c(meta_ref[,5],pred_label_l2_mtg[,1])

score <- data.frame(Acc = c(scores_l1_mtg[[1]],scores_l2_mtg[[1]]),
                    Macro_F1 = c(scores_l1_mtg[[2]],scores_l2_mtg[[2]]))
rownames(score) <- c('l1_mtg','l2_mtg')
saveRDS(score,'./score_mtg.rds')

# MIDAS label transfer
MIDAS <- as.matrix(res[[1]])
rownames(MIDAS) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),4]))
scores_l1_midas<- evaluate(meta[rownames(meta_q),4],pred_label_l1_midas[,1])
meta$pred_l1_midas <- c(meta_ref[,4],pred_label_l1_midas[,1])

pred_label_l2_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),5]))
scores_l2_midas<- evaluate(meta[rownames(meta_q),5],pred_label_l2_midas[,1])
meta$pred_l2_midas <- c(meta_ref[,5],pred_label_l2_midas[,1])

score <- data.frame(Acc = c(scores_l1_midas[[1]],scores_l2_midas[[1]]),
                    Macro_F1 = c(scores_l1_midas[[2]],scores_l2_midas[[2]]))
rownames(score) <- c('l1_midas','l2_midas')
saveRDS(score,'./score_midas.rds')

######################### CITE ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Reference-based_integration/output/Result/cite')
meta <- readRDS('./meta.rds')
method <- c('MIDAS','Multigrate',
            'Palette_harmony','Palette_harmony_supervised',
            'Palette_fastMNN','Palette_fastMNN_supervised',
            'Palette_Seurat','Palette_Seurat_supervised',
            'Palette_denovo','Palette_denovo_supervised')

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}


out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_modality.rds')


############

out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'RQ','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_RQ.rds')

# Multigrate label transfer
mtg <- as.matrix(res[[2]])
rownames(mtg) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),4]))
scores_l1_mtg<- evaluate(meta[rownames(meta_q),4],pred_label_l1_mtg[,1])
meta$pred_l1_mtg <- c(meta_ref[,4],pred_label_l1_mtg[,1])

pred_label_l2_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),5]))
scores_l2_mtg<- evaluate(meta[rownames(meta_q),5],pred_label_l2_mtg[,1])
meta$pred_l2_mtg <- c(meta_ref[,5],pred_label_l2_mtg[,1])

score <- data.frame(Acc = c(scores_l1_mtg[[1]],scores_l2_mtg[[1]]),
                    Macro_F1 = c(scores_l1_mtg[[2]],scores_l2_mtg[[2]]))
rownames(score) <- c('l1_mtg','l2_mtg')
saveRDS(score,'./score_mtg.rds')

# MIDAS label transfer
MIDAS <- as.matrix(res[[1]])
rownames(MIDAS) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),4]))
scores_l1_midas<- evaluate(meta[rownames(meta_q),4],pred_label_l1_midas[,1])
meta$pred_l1_midas <- c(meta_ref[,4],pred_label_l1_midas[,1])

pred_label_l2_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),5]))
scores_l2_midas<- evaluate(meta[rownames(meta_q),5],pred_label_l2_midas[,1])
meta$pred_l2_midas <- c(meta_ref[,5],pred_label_l2_midas[,1])

score <- data.frame(Acc = c(scores_l1_midas[[1]],scores_l2_midas[[1]]),
                    Macro_F1 = c(scores_l1_midas[[2]],scores_l2_midas[[2]]))
rownames(score) <- c('l1_midas','l2_midas')
saveRDS(score,'./score_midas.rds')


######################### CITE + Multiome ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Reference-based_integration/output/Result/cite_multi')
meta <- readRDS('./meta.rds')
method <- c('MIDAS','Multigrate',
            'Palette_harmony','Palette_harmony_supervised',
            'Palette_fastMNN','Palette_fastMNN_supervised',
            'Palette_Seurat','Palette_Seurat_supervised',
            'Palette_denovo','Palette_denovo_supervised')

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}


out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_modality.rds')


############

out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'RQ','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_RQ.rds')

# Multigrate label transfer
mtg <- as.matrix(res[[2]])
rownames(mtg) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),4]))
scores_l1_mtg<- evaluate(meta[rownames(meta_q),4],pred_label_l1_mtg[,1])
meta$pred_l1_mtg <- c(meta_ref[,4],pred_label_l1_mtg[,1])

pred_label_l2_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),5]))
scores_l2_mtg<- evaluate(meta[rownames(meta_q),5],pred_label_l2_mtg[,1])
meta$pred_l2_mtg <- c(meta_ref[,5],pred_label_l2_mtg[,1])

score <- data.frame(Acc = c(scores_l1_mtg[[1]],scores_l2_mtg[[1]]),
                    Macro_F1 = c(scores_l1_mtg[[2]],scores_l2_mtg[[2]]))
rownames(score) <- c('l1_mtg','l2_mtg')
saveRDS(score,'./score_mtg.rds')

# MIDAS label transfer
MIDAS <- as.matrix(res[[1]])
rownames(MIDAS) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),4]))
scores_l1_midas<- evaluate(meta[rownames(meta_q),4],pred_label_l1_midas[,1])
meta$pred_l1_midas <- c(meta_ref[,4],pred_label_l1_midas[,1])

pred_label_l2_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),5]))
scores_l2_midas<- evaluate(meta[rownames(meta_q),5],pred_label_l2_midas[,1])
meta$pred_l2_midas <- c(meta_ref[,5],pred_label_l2_midas[,1])

score <- data.frame(Acc = c(scores_l1_midas[[1]],scores_l2_midas[[1]]),
                    Macro_F1 = c(scores_l1_midas[[2]],scores_l2_midas[[2]]))
rownames(score) <- c('l1_midas','l2_midas')
saveRDS(score,'./score_midas.rds')

######################### Multiome ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Reference-based_integration/output/Result/multiome')
meta <- readRDS('./meta.rds')
method <- c('MIDAS','Multigrate',
            'Palette_harmony','Palette_harmony_supervised',
            'Palette_fastMNN','Palette_fastMNN_supervised',
            'Palette_Seurat','Palette_Seurat_supervised',
            'Palette_denovo','Palette_denovo_supervised')

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}


out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_modality.rds')


############

out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'RQ','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_RQ.rds')

# Multigrate label transfer
mtg <- as.matrix(res[[2]])
rownames(mtg) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),4]))
scores_l1_mtg<- evaluate(meta[rownames(meta_q),4],pred_label_l1_mtg[,1])
meta$pred_l1_mtg <- c(meta_ref[,4],pred_label_l1_mtg[,1])

pred_label_l2_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),5]))
scores_l2_mtg<- evaluate(meta[rownames(meta_q),5],pred_label_l2_mtg[,1])
meta$pred_l2_mtg <- c(meta_ref[,5],pred_label_l2_mtg[,1])

score <- data.frame(Acc = c(scores_l1_mtg[[1]],scores_l2_mtg[[1]]),
                    Macro_F1 = c(scores_l1_mtg[[2]],scores_l2_mtg[[2]]))
rownames(score) <- c('l1_mtg','l2_mtg')
saveRDS(score,'./score_mtg.rds')

# MIDAS label transfer
MIDAS <- as.matrix(res[[1]])
rownames(MIDAS) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),4]))
scores_l1_midas<- evaluate(meta[rownames(meta_q),4],pred_label_l1_midas[,1])
meta$pred_l1_midas <- c(meta_ref[,4],pred_label_l1_midas[,1])

pred_label_l2_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),5]))
scores_l2_midas<- evaluate(meta[rownames(meta_q),5],pred_label_l2_midas[,1])
meta$pred_l2_midas <- c(meta_ref[,5],pred_label_l2_midas[,1])

score <- data.frame(Acc = c(scores_l1_midas[[1]],scores_l2_midas[[1]]),
                    Macro_F1 = c(scores_l1_midas[[2]],scores_l2_midas[[2]]))
rownames(score) <- c('l1_midas','l2_midas')
saveRDS(score,'./score_midas.rds')

######################### RNA ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Reference-based_integration/output/Result/rna')
meta <- readRDS('./meta.rds')
method <- c('MIDAS','Multigrate',
            'Palette_harmony','Palette_harmony_supervised',
            'Palette_fastMNN','Palette_fastMNN_supervised',
            'Palette_Seurat','Palette_Seurat_supervised',
            'Palette_denovo','Palette_denovo_supervised')

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}


out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_modality.rds')


############

out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'RQ','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_RQ.rds')

# Multigrate label transfer
mtg <- as.matrix(res[[2]])
rownames(mtg) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),4]))
scores_l1_mtg<- evaluate(meta[rownames(meta_q),4],pred_label_l1_mtg[,1])
meta$pred_l1_mtg <- c(meta_ref[,4],pred_label_l1_mtg[,1])

pred_label_l2_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),5]))
scores_l2_mtg<- evaluate(meta[rownames(meta_q),5],pred_label_l2_mtg[,1])
meta$pred_l2_mtg <- c(meta_ref[,5],pred_label_l2_mtg[,1])

score <- data.frame(Acc = c(scores_l1_mtg[[1]],scores_l2_mtg[[1]]),
                    Macro_F1 = c(scores_l1_mtg[[2]],scores_l2_mtg[[2]]))
rownames(score) <- c('l1_mtg','l2_mtg')
saveRDS(score,'./score_mtg.rds')

# MIDAS label transfer
MIDAS <- as.matrix(res[[1]])
rownames(MIDAS) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),4]))
scores_l1_midas<- evaluate(meta[rownames(meta_q),4],pred_label_l1_midas[,1])
meta$pred_l1_midas <- c(meta_ref[,4],pred_label_l1_midas[,1])

pred_label_l2_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),5]))
scores_l2_midas<- evaluate(meta[rownames(meta_q),5],pred_label_l2_midas[,1])
meta$pred_l2_midas <- c(meta_ref[,5],pred_label_l2_midas[,1])

score <- data.frame(Acc = c(scores_l1_midas[[1]],scores_l2_midas[[1]]),
                    Macro_F1 = c(scores_l1_midas[[2]],scores_l2_midas[[2]]))
rownames(score) <- c('l1_midas','l2_midas')
saveRDS(score,'./score_midas.rds')

######################### RNA + ADT ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Reference-based_integration/output/Result/rna_adt')
meta <- readRDS('./meta.rds')
method <- c('MIDAS','Multigrate',
            'Palette_harmony','Palette_harmony_supervised',
            'Palette_fastMNN','Palette_fastMNN_supervised',
            'Palette_Seurat','Palette_Seurat_supervised',
            'Palette_denovo','Palette_denovo_supervised')

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}


out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_modality.rds')


############

out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'RQ','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_RQ.rds')

# Multigrate label transfer
mtg <- as.matrix(res[[2]])
rownames(mtg) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),4]))
scores_l1_mtg<- evaluate(meta[rownames(meta_q),4],pred_label_l1_mtg[,1])
meta$pred_l1_mtg <- c(meta_ref[,4],pred_label_l1_mtg[,1])

pred_label_l2_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),5]))
scores_l2_mtg<- evaluate(meta[rownames(meta_q),5],pred_label_l2_mtg[,1])
meta$pred_l2_mtg <- c(meta_ref[,5],pred_label_l2_mtg[,1])

score <- data.frame(Acc = c(scores_l1_mtg[[1]],scores_l2_mtg[[1]]),
                    Macro_F1 = c(scores_l1_mtg[[2]],scores_l2_mtg[[2]]))
rownames(score) <- c('l1_mtg','l2_mtg')
saveRDS(score,'./score_mtg.rds')

# MIDAS label transfer
MIDAS <- as.matrix(res[[1]])
rownames(MIDAS) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),4]))
scores_l1_midas<- evaluate(meta[rownames(meta_q),4],pred_label_l1_midas[,1])
meta$pred_l1_midas <- c(meta_ref[,4],pred_label_l1_midas[,1])

pred_label_l2_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),5]))
scores_l2_midas<- evaluate(meta[rownames(meta_q),5],pred_label_l2_midas[,1])
meta$pred_l2_midas <- c(meta_ref[,5],pred_label_l2_midas[,1])

score <- data.frame(Acc = c(scores_l1_midas[[1]],scores_l2_midas[[1]]),
                    Macro_F1 = c(scores_l1_midas[[2]],scores_l2_midas[[2]]))
rownames(score) <- c('l1_midas','l2_midas')
saveRDS(score,'./score_midas.rds')

######################### RNA + ADT + ATAC ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Reference-based_integration/output/Result/rna_adt_atac')
meta <- readRDS('./meta.rds')
method <- c('MIDAS','Multigrate',
            'Palette_harmony','Palette_harmony_supervised',
            'Palette_fastMNN','Palette_fastMNN_supervised',
            'Palette_Seurat','Palette_Seurat_supervised',
            'Palette_denovo','Palette_denovo_supervised')

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}


out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_modality.rds')


############

out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'RQ','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_RQ.rds')

# Multigrate label transfer
mtg <- as.matrix(res[[2]])
rownames(mtg) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),4]))
scores_l1_mtg<- evaluate(meta[rownames(meta_q),4],pred_label_l1_mtg[,1])
meta$pred_l1_mtg <- c(meta_ref[,4],pred_label_l1_mtg[,1])

pred_label_l2_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),5]))
scores_l2_mtg<- evaluate(meta[rownames(meta_q),5],pred_label_l2_mtg[,1])
meta$pred_l2_mtg <- c(meta_ref[,5],pred_label_l2_mtg[,1])

score <- data.frame(Acc = c(scores_l1_mtg[[1]],scores_l2_mtg[[1]]),
                    Macro_F1 = c(scores_l1_mtg[[2]],scores_l2_mtg[[2]]))
rownames(score) <- c('l1_mtg','l2_mtg')
saveRDS(score,'./score_mtg.rds')

# MIDAS label transfer
MIDAS <- as.matrix(res[[1]])
rownames(MIDAS) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),4]))
scores_l1_midas<- evaluate(meta[rownames(meta_q),4],pred_label_l1_midas[,1])
meta$pred_l1_midas <- c(meta_ref[,4],pred_label_l1_midas[,1])

pred_label_l2_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),5]))
scores_l2_midas<- evaluate(meta[rownames(meta_q),5],pred_label_l2_midas[,1])
meta$pred_l2_midas <- c(meta_ref[,5],pred_label_l2_midas[,1])

score <- data.frame(Acc = c(scores_l1_midas[[1]],scores_l2_midas[[1]]),
                    Macro_F1 = c(scores_l1_midas[[2]],scores_l2_midas[[2]]))
rownames(score) <- c('l1_midas','l2_midas')
saveRDS(score,'./score_midas.rds')


######################### RNA + ATAC ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Reference-based_integration/output/Result/rna_atac')
meta <- readRDS('./meta.rds')
method <- c('MIDAS','Multigrate',
            'Palette_harmony','Palette_harmony_supervised',
            'Palette_fastMNN','Palette_fastMNN_supervised',
            'Palette_Seurat','Palette_Seurat_supervised',
            'Palette_denovo','Palette_denovo_supervised')

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}


out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_modality.rds')


############

out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'RQ','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("MIDAS",'Multigrate',
                   'Palette_harmony','Palette_harmony_supervised',
                   'Palette_fastMNN','Palette_fastMNN_supervised',
                   "Palette_Seurat","Palette_Seurat_supervised",
                   'Palette_denovo','Palette_denovo_supervised')
saveRDS(out,'./lisi_score_RQ.rds')

# Multigrate label transfer
mtg <- as.matrix(res[[2]])
rownames(mtg) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),4]))
scores_l1_mtg<- evaluate(meta[rownames(meta_q),4],pred_label_l1_mtg[,1])
meta$pred_l1_mtg <- c(meta_ref[,4],pred_label_l1_mtg[,1])

pred_label_l2_mtg <- Trans.Knowledge(mtg[rownames(meta_ref),],
                                     mtg[rownames(meta_q),],
                                     data.frame(meta[rownames(meta_ref),5]))
scores_l2_mtg<- evaluate(meta[rownames(meta_q),5],pred_label_l2_mtg[,1])
meta$pred_l2_mtg <- c(meta_ref[,5],pred_label_l2_mtg[,1])

score <- data.frame(Acc = c(scores_l1_mtg[[1]],scores_l2_mtg[[1]]),
                    Macro_F1 = c(scores_l1_mtg[[2]],scores_l2_mtg[[2]]))
rownames(score) <- c('l1_mtg','l2_mtg')
saveRDS(score,'./score_mtg.rds')

# MIDAS label transfer
MIDAS <- as.matrix(res[[1]])
rownames(MIDAS) <- rownames(meta)
meta_ref <- dplyr::filter(meta,RQ == 'reference')
meta_q <- dplyr::filter(meta,RQ == 'query')

pred_label_l1_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),4]))
scores_l1_midas<- evaluate(meta[rownames(meta_q),4],pred_label_l1_midas[,1])
meta$pred_l1_midas <- c(meta_ref[,4],pred_label_l1_midas[,1])

pred_label_l2_midas <- Trans.Knowledge(MIDAS[rownames(meta_ref),],
                                       MIDAS[rownames(meta_q),],
                                       data.frame(meta[rownames(meta_ref),5]))
scores_l2_midas<- evaluate(meta[rownames(meta_q),5],pred_label_l2_midas[,1])
meta$pred_l2_midas <- c(meta_ref[,5],pred_label_l2_midas[,1])

score <- data.frame(Acc = c(scores_l1_midas[[1]],scores_l2_midas[[1]]),
                    Macro_F1 = c(scores_l1_midas[[2]],scores_l2_midas[[2]]))
rownames(score) <- c('l1_midas','l2_midas')
saveRDS(score,'./score_midas.rds')



