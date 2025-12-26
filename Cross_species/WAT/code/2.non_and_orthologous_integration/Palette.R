source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")
options(future.globals.maxSize = 500 * 1024^3)
library(Matrix)

setwd('./human_mouse_WAT')
meta <- readRDS('./meta.rds')
human_homo <- as(as.matrix(readRDS('./human/rna_homo.rds')),'dgCMatrix')
human_non <- as(as.matrix(readRDS('./human/rna.rds')),'dgCMatrix')
human_gene <- readRDS('./human/gene_non.rds')

mouse_homo <- as(as.matrix(readRDS('./mouse/rna_homo.rds')),'dgCMatrix')
mouse_non <- as(as.matrix(readRDS('./mouse/rna.rds')),'dgCMatrix')
mouse_gene <- readRDS('./mouse/gene_non.rds')

human_meta <- readRDS('./human/meta.rds')
mouse_meta <- readRDS('./mouse/meta.rds')

human_homo_list <- list()
human_non_list <- list()
for(i in unique(human_meta$Batch)){
  sub_meta <- dplyr::filter(human_meta,Batch == i)
  human_homo_list[[i]] <- human_homo[,rownames(sub_meta)]
  human_non_list[[i]] <- human_non[,rownames(sub_meta)]
}

mouse_homo_list <- list()
mouse_non_list <- list()
for(i in unique(mouse_meta$Batch)){
  sub_meta <- dplyr::filter(mouse_meta,Batch == i)
  mouse_homo_list[[i]] <- mouse_homo[,rownames(sub_meta)]
  mouse_non_list[[i]] <- mouse_non[,rownames(sub_meta)]
}
rm(sub_meta,human_homo,human_non,mouse_homo,mouse_non)

names(human_homo_list) <- NULL
names(human_non_list) <- NULL
names(mouse_homo_list) <- NULL
names(mouse_non_list) <- NULL

musa_obj <- Create.Musa.Object(data.list = c(human_homo_list,human_non_list,mouse_homo_list,mouse_non_list),
                               samples =  c(unique(human_meta$Batch),
                                            unique(human_meta$Batch),
                                            unique(mouse_meta$Batch),
                                            unique(mouse_meta$Batch)), 
                               modals = c(rep('homo',length(human_homo_list)),
                                          rep('hnon',length(human_homo_list)),
                                          rep('homo',length(mouse_homo_list)),
                                          rep('mnon',length(mouse_homo_list))),
                               filter = F,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("homo",'hnon','mnon'), 
                               sparce = TRUE)
rm(human_homo_list,human_non_list,mouse_homo_list,mouse_non_list)

musa_obj <- Normalize.Data(musa_obj, 
                           modal = "homo",
                           normal.method = "LogNormalize")
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "hnon",
                           normal.method = "LogNormalize")
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "mnon",
                           normal.method = "LogNormalize")
musa_obj <- Add.HVFs(musa_obj,modal = "homo")
musa_obj <- Add.HVFs(musa_obj,modal = "hnon")
musa_obj <- Add.HVFs(musa_obj,modal = "mnon")

meta <- meta[musa_obj@Meta$cell.name,]
musa_obj@Meta[['celltype']] <- meta$CellType
all.equal(musa_obj@Meta$cell.name,rownames(meta))

musa_obj <- DownSample(musa_obj,
                       modal = c("homo",'hnon','mnon'),
                       supervised = T,
                       group = 'celltype',
                       dims = list(1:50,1:50,1:50),
                       method = c("PCA","PCA","PCA"))

musa_obj <- Add.HVFs(musa_obj,modal = "hnon", features = human_gene)
musa_obj <- Add.HVFs(musa_obj,modal = "mnon", features = mouse_gene)

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("homo",'hnon','mnon'),
                           lambda = list(0.8,0.8,0.8),
                           supervised = T,
                           joint = F,
                           L2.norm = F,
                           proj.scale = F,
                           sub.dims = list(1:40,1:10,1:10),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,15,15),
                           max.Angle = c(50,50,50),
                           pre_dims = list(1:3000))

musa_obj <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,scale = T,
                         modal.norm = F,emb_L2 =T,supervised = T,group = 'celltype')


emb <- musa_obj@Int.result$bind
rm(musa_obj)
gc()
emb <- emb[rownames(meta),]
meta <- readRDS('./Cross_species/WAT/output/meta.rds')
write.csv(as.matrix(emb),
          paste0(out_dir,'/Palette_non.csv'),
          row.names = FALSE)