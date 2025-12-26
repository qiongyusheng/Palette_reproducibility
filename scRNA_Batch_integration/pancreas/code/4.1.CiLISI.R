library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('/home/server/sqy/data/mosaic/result/NC_revise_1/Batch_Bench/human_pancreas')
method = c("Palette_unsupervised",
                 "Palette_supervised",
                 'STACAS',
                 'ssSTACAS',
                 'SIGNAL',
                 'scpoli',
                 'scanvi',
                 'fastMNN',
                 'Harmony',
                'Seurat',
                 'raw',
                  
                 'label_miss/Palette_miss5',
                 'label_miss/Palette_miss10',
                 'label_miss/Palette_miss20',
                 'label_miss/SIGNAL_miss5',
                 'label_miss/SIGNAL_miss10',
                 'label_miss/SIGNAL_miss20',
                 'label_miss/ssSTACAS_miss5',
                 'label_miss/ssSTACAS_miss10',
                 'label_miss/ssSTACAS_miss20',
                 'label_miss/scanvi_miss5',
                 'label_miss/scanvi_miss10',
                 'label_miss/scanvi_miss20',
                 'label_miss/scpoli_miss5',
                 'label_miss/scpoli_miss10',
                 'label_miss/scpoli_miss20',
                  
                 'label_shuffled/Palette_shuffled5',
                 'label_shuffled/Palette_shuffled10',
                 'label_shuffled/Palette_shuffled20',
                 'label_shuffled/ssSTACAS_shuffled5',
                 'label_shuffled/ssSTACAS_shuffled10',
                 'label_shuffled/ssSTACAS_shuffled20',
                 'label_shuffled/SIGNAL_shuffled5',
                 'label_shuffled/SIGNAL_shuffled10',
                 'label_shuffled/SIGNAL_shuffled20',
                 'label_shuffled/scpoli_shuffled5',
                 'label_shuffled/scpoli_shuffled10',
                 'label_shuffled/scpoli_shuffled20',
                 'label_shuffled/scanvi_shuffled5',
                 'label_shuffled/scanvi_shuffled10',
                'label_shuffled/scanvi_shuffled20')

res <- list()
for(i in 1:length(method)){
    res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}
meta <- readRDS('./meta.rds')

out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'Batch','CellType',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}
saveRDS(out,'./lisi_score.rds')
