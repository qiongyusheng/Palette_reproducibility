######################### Random 1 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/TEA_s2/output/random1')

method = c( "Palette_v1",
            'Multigrate_v1',
            'scMoMaT_v1',
            'midas_v1',
            'scVAEIT_v1',
            "stab_v1")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta_v1.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'batch','celltype',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS','scVAEIT',"StabMap")

saveRDS(out,'./lisi_score_v1.rds')

######################### Random 2 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/TEA_s2/output/random3')

method = c( "Palette_v2",
            'Multigrate_v2',
            'scMoMaT_v2',
            'midas_v2',
            'scVAEIT_v2',
            "stab_v2")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta_v2.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'batch','celltype',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS','scVAEIT',"StabMap")

saveRDS(out,'./lisi_score_v2.rds')

######################### Random 3 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/TEA_s2/output/random2')

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
  lisi <- CiLISI_sample(res[[i]],meta,'batch','celltype',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS','scVAEIT',"StabMap")

saveRDS(out,'./lisi_score.rds')


######################### Random 4 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/TEA_s2/output/random4')

method = c( "Palette_v3",
            'Multigrate_v3',
            'scMoMaT_v3',
            'midas_v3',
            'scVAEIT_v3',
            "stab_v3")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta_v3.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'batch','celltype',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS','scVAEIT',"StabMap")

saveRDS(out,'./lisi_score_v3.rds')

######################### Random 5 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Unsupervised/TEA_s2/output/random5')

method = c( "Palette_v4",
            'Multigrate_v4',
            'scMoMaT_v4',
            'midas_v4',
            'scVAEIT_v4',
            "stab_v4")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}

meta <- readRDS('./meta_v4.rds')
out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'batch','celltype',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS','scVAEIT',"StabMap")

saveRDS(out,'./lisi_score_v4.rds')

