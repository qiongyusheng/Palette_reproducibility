######################### Random 1 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/BMMC_s2/output/random1')

method = c( "Palette_34",
            'Multigrate_34',
            'scMoMaT_34',
            'midas_34',
            'scVAEIT_34',
            "stab_34")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta_34.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS','scVAEIT',"StabMap")

saveRDS(out,'./lisi_score_34.rds')

######################### Random 2 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/BMMC_s2/output/random2')

method = c( "Palette_35",
            'Multigrate_35',
            'scMoMaT_35',
            'midas_35',
            "stab_35")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta_35.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS',"StabMap")

saveRDS(out,'./lisi_score_35.rds')

######################### Random 3 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/BMMC_s2/output/random3')

method = c( "Palette_36",
            'Multigrate_36',
            'scMoMaT_36',
            'midas_36',
            "stab_36")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta_36.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS',"StabMap")

saveRDS(out,'./lisi_score_36.rds')

######################### Random 4 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/BMMC_s2/output/random4')

method = c( "Palette_24",
            'Multigrate_24',
            'scMoMaT_24',
            'midas_24',
            "stab_24")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta_24.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS',"StabMap")

saveRDS(out,'./lisi_score_24.rds')

######################### Random 5 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/BMMC_s2/output/random5')

method = c( "Palette_26",
            'Multigrate_26',
            'scMoMaT_26',
            'midas_26',
            "stab_26")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta_26.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype.l2',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS',"StabMap")

saveRDS(out,'./lisi_score_26.rds')

