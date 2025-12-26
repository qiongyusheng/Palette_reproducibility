######################### Random 1 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/BMMC_s1/output/random1')

method = c( "Palette_drop5",
            'Multigrate_drop5',
            'scMoMaT_drop5',
            'midas_drop5',
            'scVAEIT_drop5',
            "stab_drop5")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta_drop5.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS','scVAEIT',"StabMap")

saveRDS(out,'./lisi_score_drop5.rds')

######################### Random 2 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/BMMC_s1/output/random2')

method = c( "Palette_drop1",
            'Multigrate_drop1',
            'scMoMaT_drop1',
            'midas_drop1',
            'scVAEIT_drop1',
            "stab_drop1")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta_drop1.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS','scVAEIT',"StabMap")

saveRDS(out,'./lisi_score_drop1.rds')

######################### Random 3 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/BMMC_s1/output/random3')

method = c( "Palette_drop6",
            'Multigrate_drop6',
            'scMoMaT_drop6',
            'midas_drop6',
            'scVAEIT_drop6',
            "stab_drop6")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta_drop6.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS','scVAEIT',"StabMap")

saveRDS(out,'./lisi_score_drop6.rds')

######################### Random 4 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/BMMC_s1/output/random4')

method = c( "Palette_drop4",
            'Multigrate_drop4',
            'scMoMaT_drop4',
            'midas_drop4',
            'scVAEIT_drop4',
            "stab_drop4")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta_drop4.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS','scVAEIT',"StabMap")

saveRDS(out,'./lisi_score_drop4.rds')

######################### Random 5 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/BMMC_s1/output/random5')

method = c( "Palette_drop2",
            'Multigrate_drop2',
            'scMoMaT_drop2',
            'midas_drop2',
            'scVAEIT_drop2',
            "stab_drop2")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta_drop2.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS','scVAEIT',"StabMap")

saveRDS(out,'./lisi_score_drop2.rds')

