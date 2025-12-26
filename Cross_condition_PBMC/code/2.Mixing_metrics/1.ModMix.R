library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')

setwd('./Cross_condition_PBMC/output/Embedding')
method = c( "Palette",
                  'Multigrate',
                  'scMoMaT','midas',
           'scVAEIT',
           "stab")
res <- list()
for(i in 1:length(method)){
    res[[i]] <- fread(paste0('./Cross_condition_PBMC/output/Embedding/',method[i],'.csv'),data.table = F,header = T)
}
meta <- readRDS('./meta.rds')

out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'modality','celltype',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

setwd('./Cross_condition_PBMC/output/Embedding')
rownames(out) <- c( "Palette",'Multigrate','scMoMaT','MIDAS','scVAEIT',"StabMap")
saveRDS(out,'./lisi_score.rds')
