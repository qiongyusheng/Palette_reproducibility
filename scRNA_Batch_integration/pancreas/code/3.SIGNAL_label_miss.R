library(scInt)
library(RSpectra)
library(SIGNAL)

setwd('./')
data <- readRDS('./dataset.rds')
varg <- readRDS('./vargs.rds')
data <- lapply(data,function(x) x[varg,])
meta <- readRDS('./meta.rds')
X <- do.call(cbind,data)
rm(data)
out_dir <- './label_miss/'
task <- c('label_miss5s','label_miss10s','label_miss20s')
mod <- c('_miss5','_miss10','_miss20')


for(i in 1:3){
    
  res_SIGNAL = Run.gcPCA(X, meta, g_factor = task[i], b_factor = "Batch", 
                         excluded.cells = which(meta[,task[i]] == 'Unknown'),
                         npcs = 20, lambda = 50,do.scale = F,do.cosine = T)
  write.csv(as.matrix(t(res_SIGNAL)),
            paste0(out_dir,'SIGNAL',mod[i],'.csv'),
            row.names = FALSE)
}
