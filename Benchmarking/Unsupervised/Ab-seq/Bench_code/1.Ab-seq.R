######################### Random 1 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/Ab-seq/output/random1')

method = c( "Palette_s1",
            'Multigrate_s1',
            'scMoMaT_s1',
            'midas_s1',
            'scVAEIT_s1',
            "stab_s1")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta_s1.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS','scVAEIT',"StabMap")

saveRDS(out,'./lisi_score_s1.rds')

######################### Random 2 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/Ab-seq/output/random2')

method = c( "Palette_s2",
            'Multigrate_s2',
            'scMoMaT_s2',
            'midas_s2',
            'scVAEIT_s2',
            "stab_s2")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta_s2.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS','scVAEIT',"StabMap")

saveRDS(out,'./lisi_score_s2.rds')

######################### Random 3 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/Ab-seq/output/random3')

method = c( "Palette_s3",
            'Multigrate_s3',
            'scMoMaT_s3',
            'midas_s3',
            'scVAEIT_s3',
            "stab_s3")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta_s3.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS','scVAEIT',"StabMap")

saveRDS(out,'./lisi_score_s3.rds')

######################### Random 4 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/Ab-seq/output/random4')

method = c( "Palette",
            'Multigrate',
            'scMoMaT',
            'midas',
            'scVAEIT',
            "stab")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS','scVAEIT',"StabMap")

saveRDS(out,'./lisi_score.rds')

######################### Random 5 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/Ab-seq/output/random5')

method = c( "Palette_s4",
            'Multigrate_s4',
            'scMoMaT_s4',
            'midas_s4',
            'scVAEIT_s4',
            "stab_s4")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta_s4.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS','scVAEIT',"StabMap")

saveRDS(out,'./lisi_score_s4.rds')

