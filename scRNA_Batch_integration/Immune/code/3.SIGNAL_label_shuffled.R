library(scInt)
library(RSpectra)
library(SIGNAL)

setwd('./')
data <- readRDS('./dataset.rds')
varg <- readRDS('./vargs.rds')
data <- lapply(data,function(x) x[varg,])
meta <- readRDS('./meta.rds')
X = do.call(cbind,data)
rm(data)
out_dir <- './label_shuffled/'
task <- c('label_shuffled5','label_shuffled10','label_shuffled20')
mod <- c('_shuffled5','_shuffled10','_shuffled20')
for(i in 1:3){
  res_SIGNAL = Run.gcPCA(X, meta, g_factor = task[i], b_factor = "Batch", 
                         npcs = 20, lambda = 50,do.scale = F,do.cosine = T)
  write.csv(as.matrix(t(res_SIGNAL)),
            paste0(out_dir,'SIGNAL',mod[i],'.csv'),
            row.names = FALSE)
}


