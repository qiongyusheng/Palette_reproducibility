source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
options(future.globals.maxSize = 500 * 1024^3)

setwd('./')

out_dir <- './label_shuffled'
data <- readRDS('./dataset.rds')
varg <- readRDS('./vargs.rds')
data <- lapply(data,function(x) x[varg,])
meta <- readRDS("./meta.rds")

names(data) <- NULL
musa_obj <- Create.Musa.Object(data.list = data, 
                               samples = c('B1','B2','B3','B4','B5'), 
                               modals = c('rna','rna','rna','rna','rna'),
                               filter=T,
                               min.cells = c(0),  
                               min.features = c(1),
                               modal.order = c('rna'), 
                               sparce = TRUE)
rm(data)
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- Add.HVFs(musa_obj, modal = "rna")
cellname <- rownames(meta)
meta <- meta[musa_obj@Meta$cell.name,]
all.equal(musa_obj@Meta$cell.name,rownames(meta))

task <- c('label_shuffled5','label_shuffled10','label_shuffled20')
mod <- c('_shuffled5','_shuffled10','_shuffled20')
for(i in 1:3){
  obj <- musa_obj
  obj@Meta[['celltype']] <- meta[,task[i]]
  obj <- DownSample(obj,
                    modal = c("rna"),
                    supervised = T,
                    group = 'celltype',
                    dims = list(1:50),# 使用哪些维数
                    method = c("PCA"))
  
  obj <- Find.Subspace3(obj,
                        modal = c("rna"),
                        lambda = list(0.8),
                        alpha = list(0),
                        supervised = T,
                        L2.norm = F,
                        sub.dims = list(1:20),
                        cos.dims = list(1:50),
                        Angle.var = c(15),
                        max.Angle = c(50))
  
  emb <- list()
  for(j in 1:5){
    emb[[j]] <- obj@Assays[[j]]@embedding
  }
  emb <- Reduce(cbind,emb)
  emb <- emb[,cellname]
  write.csv(as.matrix(t(emb)),
            paste0(out_dir,'/Palette',paste0(mod[i]),'.csv'),
            row.names = FALSE)
}