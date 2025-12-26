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

for(i in 1:3){
    
   query_meta <- meta[is.na(meta[,task[i]]), ]
   ref_meta <- meta[!is.na(meta[,task[i]]), ]
   query <- X[,rownames(query_meta)] 
   ref <- X[,rownames(ref_meta)] 
   zero_m = sparseMatrix(i = integer(0), j = integer(0), dims = c(nrow(ref) - nrow(query), ncol(query)))
   rownames(zero_m) = setdiff(rownames(ref), rownames(query))
   colnames(zero_m) = colnames(query)
   query = rbind(query, zero_m)[rownames(ref), ]
   predicted_labels = Run.LabelTransfer.Single(ref, ref_meta[,task[i]], query)
    query_meta[,task[i]] <- predicted_labels$Prediction
    tmp_meta <- rbind(query_meta,ref_meta)
    tmp_meta <- tmp_meta[rownames(meta),]
    
  res_SIGNAL = Run.gcPCA(X, tmp_meta, g_factor = task[i], b_factor = "Batch", 
                         npcs = 20, lambda = 50,do.scale = F,do.cosine = T)
  write.csv(as.matrix(t(res_SIGNAL)),
            paste0(out_dir,'SIGNAL',mod[i],'.csv'),
            row.names = FALSE)
}
