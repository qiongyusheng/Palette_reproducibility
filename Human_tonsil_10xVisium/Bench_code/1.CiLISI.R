######################### Random 1 ################################
library(data.table)
library(Matrix)
library(matrixStats)
source('./metrics/metrics_R.R')
setwd('./10xVisium_tonsil')
method = c("Palette", 
           'Multigrate',
           'scMoMaT', 
           'midas','scVAEIT',"stab")

res <- list()
for(i in 1:length(method)){
  res[[i]] <- fread(paste0('./',method[i],'.csv'),data.table = F,header = T)
}
meta <- readRDS('./10xVisium_tonsil/meta.rds')

out <- data.frame(CiLISI = rep(0,length(method)))
rownames(out) <- method
for(i in 1:length(method)){
  lisi <- CiLISI_sample(res[[i]],meta,'batch','anno',n_repeat = 4)
  out[method[i],'CiLISI'] <- lisi
}

rownames(out) <- c("Palette",
                   'Multigrate',
                   'scMoMaT',
                   'MIDAS','scVAEIT',
                   "StabMap")

saveRDS(out,'./10xVisium_tonsil/lisi_score.rds')