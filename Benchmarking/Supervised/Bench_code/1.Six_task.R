######################### TEA_s1 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Supervised/output/TEA_s1')

data_name <- c('Palette_v1','Palette','Palette_v2','Palette_v3','Palette_v4')
meta_name <- c('meta_v1','meta','meta_v2','meta_v3','meta_v4')

out <- data.frame(CiLISI = rep(0,length(data_name)))
rownames(out) <- data_name

for(i in 1:length(data_name)){
  res <- fread(paste0('./',data_name[i],'.csv'),data.table = F,header = T)
  meta <- readRDS(paste0('./',meta_name[i],'.rds'))
  lisi <- CiLISI_sample(res,meta,'batch','celltype',n_repeat = 4)
  out[data_name[i],'CiLISI'] <- lisi
}

saveRDS(out,'./lisi_score.rds')

######################### TEA_s2 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Supervised/output/TEA_s2')

data_name <- c('Palette_v1','Palette_v2','Palette','Palette_v3','Palette_v4')
meta_name <- c('meta_v1','meta_v2','meta','meta_v3','meta_v4')

out <- data.frame(CiLISI = rep(0,length(data_name)))
rownames(out) <- data_name

for(i in 1:length(data_name)){
  res <- fread(paste0('./',data_name[i],'.csv'),data.table = F,header = T)
  meta <- readRDS(paste0('./',meta_name[i],'.rds'))
  lisi <- CiLISI_sample(res,meta,'batch','celltype',n_repeat = 4)
  out[data_name[i],'CiLISI'] <- lisi
}

saveRDS(out,'./lisi_score.rds')

######################### BMMC_s1 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Supervised/output/BMMC_s1')
data_name <- c('Palette_drop5','Palette_drop1','Palette_drop6','Palette_drop4','Palette_drop2')
meta_name <- c('meta_drop5','meta_drop1','meta_drop6','meta_drop4','meta_drop2')

out <- data.frame(CiLISI = rep(0,length(data_name)))
rownames(out) <- data_name

for(i in 1:length(data_name)){
  res <- fread(paste0('./',data_name[i],'.csv'),data.table = F,header = T)
  meta <- readRDS(paste0('./',meta_name[i],'.rds'))
  lisi <- CiLISI_sample(res,meta,'modality','celltype.l2',n_repeat = 4)
  out[data_name[i],'CiLISI'] <- lisi
}
saveRDS(out,'./lisi_score.rds')

######################### BMMC_s2 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Supervised/output/BMMC_s2')

data_name <- c('Palette_34','Palette_35','Palette_36','Palette_24','Palette_26')
meta_name <- c('meta_34','meta_35','meta_36','meta_24','meta_26')

out <- data.frame(CiLISI = rep(0,length(data_name)))
rownames(out) <- data_name

for(i in 1:length(data_name)){
  res <- fread(paste0('./',data_name[i],'.csv'),data.table = F,header = T)
  meta <- readRDS(paste0('./',meta_name[i],'.rds'))
  lisi <- CiLISI_sample(res,meta,'modality','celltype.l2',n_repeat = 4)
  out[data_name[i],'CiLISI'] <- lisi
}

saveRDS(out,'./lisi_score.rds')

######################### Retina ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Supervised/output/Retina')
data_name <- c('Palette_s1','Palette_s2','Palette','Palette_s3','Palette_s4')
meta_name <- c('meta_s1','meta_s2','meta','meta_s3','meta_s4')

out <- data.frame(CiLISI = rep(0,length(data_name)))
rownames(out) <- data_name

for(i in 1:length(data_name)){
  res <- fread(paste0('./',data_name[i],'.csv'),data.table = F,header = T)
  meta <- readRDS(paste0('./',meta_name[i],'.rds'))
  lisi <- CiLISI_sample(res,meta,'modality','celltype',n_repeat = 4)
  out[data_name[i],'CiLISI'] <- lisi
}

saveRDS(out,'./lisi_score.rds')

######################### Ab-seq ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Benchmarking/Supervised/output/Ab-seq')
data_name <- c('Palette_s1','Palette_s2','Palette_s3','Palette','Palette_s4')
meta_name <- c('meta_s1','meta_s2','meta_s3','meta','meta_s4')

out <- data.frame(CiLISI = rep(0,length(data_name)))
rownames(out) <- data_name

for(i in 1:length(data_name)){
  res <- fread(paste0('./',data_name[i],'.csv'),data.table = F,header = T)
  meta <- readRDS(paste0('./',meta_name[i],'.rds'))
  lisi <- CiLISI_sample(res,meta,'modality','celltype.l2',n_repeat = 4)
  out[data_name[i],'CiLISI'] <- lisi
}

saveRDS(out,'./lisi_score.rds')