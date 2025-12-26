######################### Palette DEFs ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

dge_atac <- readRDS("./DGE/peak_new.rds")
adt_dge <- readRDS('./DGE/adt_new.rds')
RNA_DGE_rank <- readRDS('./DGE/rna_new.rds')

adt_drop <- adt_dge[,1][adt_dge$padj < 0.01 & (abs(adt_dge$foldchange) > 0.5)]
rna_drop <- RNA_DGE_rank$gene[RNA_DGE_rank$FDR < 0.01 & (abs(RNA_DGE_rank$logFC) > 0.4)]
atac_drop <- rownames(dge_atac)[dge_atac$FDR < 0.01 & (abs(dge_atac$logFC) > 0.5)]

setwd("./")
path <- c('./PBMC_ASAP/control_asap','./PBMC_ASAP/stim_asap',
          './PBMC_CITE/control_cite','./PBMC_CITE/stim_cite')
atac.list <- list()
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in 1:4){
    meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
    tmp <- readRDS(paste0(path[i],'/adt.rds'))
    idx <- match(adt_drop,rownames(tmp))
    tmp <- tmp[-idx,]
    adt.list[[i]] <- tmp
    rm(tmp,idx)
    if(i %in% c(1,2)){
        tmp <- readRDS(paste0(path[i],'/atac.rds'))
        idx <- match(atac_drop,rownames(tmp))
        tmp <- tmp[-idx,]
        atac.list[[i]] <- tmp
        print(paste0("atac drop feature are ",length(idx)))
        rm(idx,tmp)
    }
    if(i %in% c(3,4)){
        tmp <- readRDS(paste0(path[i],'/rna.rds'))
        co_rna <- intersect(rna_drop,rownames(tmp))
        idx <- match(co_rna,rownames(tmp))
        tmp <- tmp[-idx,]
        rna.list[[i-2]] <- tmp
        print(paste0("rna drop feature are ",length(idx)))
        rm(tmp,idx)
    }
}

musa_obj <- Create.Musa.Object(data.list = c(adt.list,atac.list,rna.list), 
                               samples = c('control_asap','stim_asap',
                                           'control_cite','stim_cite',
                                          'control_asap','stim_asap',
                                           'control_cite','stim_cite'), 
                               modals = c('adt','adt','adt','adt',
                                         'atac','atac','rna','rna'),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna",'atac','adt'), 
                               sparce = TRUE)

musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                          normal.method = "LogNormalize")
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "adt",
                          normal.method = "CLR",
                          margin = 2)
musa_obj <- IDFLog_try(musa_obj,modal = "atac")

musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj, modal = "adt")
musa_obj <- Add.HVFs(musa_obj, modal = "atac")

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'atac','adt'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:50,1:50),# 使用哪些维数
                         method = c("PCA","PCA",'LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

musa_obj <- Find.Subspace3(musa_obj,
                          modal = c("rna",'atac','adt'),
                          lambda = list(0.8,0.5,0.8),
                          alpha = list(0,0,0),
                           joint = T,
                           L2.norm = F,
                         sub.dims = list(1:40,1:30,1:30),
                         cos.dims = list(1:50,1:50,1:50),
                          Angle.var = c(15,20,15),
                         max.Angle = c(50,50,50))

a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)

write.csv(as.matrix(a@Int.result[["bind"]]),
          "./Palette_newDGE_drop.csv",
          row.names = FALSE)

mm <- data.frame(batch = c(meta.list[[1]]$batch,meta.list[[2]]$batch,
                           meta.list[[3]]$batch,meta.list[[4]]$batch),
                 modality = c(meta.list[[1]]$modality,meta.list[[2]]$modality,
                              meta.list[[3]]$modality,meta.list[[4]]$modality),
                 condition = c(meta.list[[1]]$condition,meta.list[[2]]$condition,
                               meta.list[[3]]$condition,meta.list[[4]]$condition),
                 celltype = c(as.character(meta.list[[1]][,1]),
                              as.character(meta.list[[2]][,1]),
                              as.character(meta.list[[3]][,1]),
                              as.character(meta.list[[4]][,1]))
)
rownames(mm) <- colnames(Reduce(cbind,adt.list))
mm_sub <- dplyr::filter(mm,celltype == 'T')
idx <- match(rownames(mm_sub),rownames(mm))
library(data.table)
origion_int <- fread('./Embedding/Palette.csv',
                 data.table = F)
newDGE_int <- a@Int.result[["bind"]][idx,]
origion_int <- origion_int[idx,]
write.csv(as.matrix(newDGE_int),
          "./T_cell_analysis/newDGE.csv",
          row.names = FALSE)
write.csv(as.matrix(origion_int),
          "./T_cell_analysis/origion.csv",
          row.names = FALSE)
saveRDS(mm_sub,'./T_cell_analysis/meta.rds')


######################### Mimitou et al. DEFs ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")


setwd('./Cross_condition_PBMC/output/Mimitou_et_al.')
dge_rna <- readRDS("./DGE_RNA_Tcell_Activation.rds")
adt_dge <- readRDS('./DGE_Protein_Tcell_Activation.rds')
adt_drop <- adt_dge[adt_dge$padj < 0.01,]
bool_adt <- rep(FALSE,160)
bool_adt[order(abs(adt_drop$foldchange),decreasing = TRUE)[1:84]] <- TRUE
rna <- readRDS('/home/server/sqy/data/mosaic/input/PBMC_CITE/control_cite/rna.rds')
co_gene <- intersect(dge_rna$gene,rownames(rna))
RNA_DGE_rank <- dplyr::filter(dge_rna,gene %in% co_gene)
rna_drop <- RNA_DGE_rank[RNA_DGE_rank$FDR < 0.01,]
bool_rna <- rep(FALSE,1087)
bool_rna[order(abs(rna_drop$logFC),decreasing = TRUE)[1:492]] <- TRUE

setwd('./')
dir1 <- './PBMC_ASAP/control_asap'
dir2 <- './PBMC_ASAP/stim_asap'
dir3 <- './PBMC_CITE/control_cite'
dir4 <- './PBMC_CITE/stim_cite'

adt1 <- readRDS(paste0(dir1,'/adt.rds'))
adt2 <- readRDS(paste0(dir2,'/adt.rds'))
adt3 <- readRDS(paste0(dir3,'/adt.rds'))
adt4 <- readRDS(paste0(dir4,'/adt.rds'))


adt_name <- adt_drop[,1]
adt_name_replaced <- gsub("-", ".",adt_name)
adt_name_replaced <- gsub("/", ".",adt_name_replaced)
adt_name_replaced <- gsub("-|_|\\(|\\)|/", ".",adt_name_replaced)
idx <- which(adt_name_replaced == 'CD57.Recombinant')
adt_name_replaced[idx] <- 'CD57_Recombinant'

idx <- which(adt_name_replaced == 'Rat.IgG1k.isotypeCtrl')
adt_name_replaced[idx] <- 'Rat_IgG1k_isotypeCtrl'

idx <- which(adt_name_replaced == 'Rat.IgG2c.IsotypeCtrl')
adt_name_replaced[idx] <- 'Rat_IgG2c_IsotypeCtrl'

idx <- match(adt_name_replaced,rownames(adt1))
drop_adt <- idx[bool_adt]


adt1 <- adt1[-drop_adt,]
adt2 <- adt2[-drop_adt,]
adt3 <- adt3[-drop_adt,]
adt4 <- adt4[-drop_adt,]

dir1 <- './PBMC_CITE/control_cite'
dir2 <- './PBMC_CITE/stim_cite'
rna1 <- readRDS(paste0(dir1,'/rna.rds'))
rna2 <- readRDS(paste0(dir2,'/rna.rds'))

idx <- match(rownames(rna_drop),rownames(rna1))
drop_rna <- idx[bool_rna]
rna1 <- rna1[-drop_rna,]
rna2 <- rna2[-drop_rna,]

setwd('./Cross_condition_PBMC/output/Mimitou_et_al.')
diff_atac <- readRDS("./DGE_ATAC_Tcell_Activation.rds")

setwd("./")
dir1 <- './PBMC_ASAP/control_asap'
dir2 <- './PBMC_ASAP/stim_asap'
atac1 <- readRDS(paste0(dir1,'/atac.rds'))
atac2 <- readRDS(paste0(dir2,'/atac.rds'))
co_peak <- intersect(rownames(diff_atac),rownames(atac1))
diff_atac <- diff_atac[co_peak,]

dge_atac <- diff_atac
atac_drop <- dge_atac[dge_atac$FDR < 0.01,]
bool_atac <- rep(FALSE,nrow(atac_drop))
bool_atac[order(abs(atac_drop$logFC),decreasing = TRUE)[1:9548]] <- TRUE
idx <- match(rownames(atac_drop),rownames(atac1))
drop_atac <- idx[bool_atac]

atac1 <- atac1[-drop_atac,]
atac2 <- atac2[-drop_atac,]

musa_obj <- Create.Musa.Object(data.list = list(adt1,adt2,adt3,adt4,atac1,atac2,rna1,rna2), 
                               samples = c('control_asap','stim_asap',
                                           'control_cite','stim_cite',
                                           'control_asap','stim_asap',
                                           'control_cite','stim_cite'), 
                               modals = c('adt','adt','adt','adt',
                                          'atac','atac','rna','rna'),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna",'atac','adt'), 
                               sparce = TRUE)

musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "adt",
                           normal.method = "CLR",
                           margin = 2)
musa_obj <- IDFLog_try(musa_obj,modal = "atac")

musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj, modal = "adt")
musa_obj <- Add.HVFs(musa_obj, modal = "atac")
musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'atac','adt'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:50,1:50),# 使用哪些维数
                         method = c("PCA","PCA",'LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)
musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'atac','adt'),
                           lambda = list(0.8,0.5,0.8),
                           alpha = list(0,0,0),
                           joint = T,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,20,15),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
setwd("./")
path <- c('./PBMC_ASAP/control_asap','./PBMC_ASAP/stim_asap',
          './PBMC_CITE/control_cite','./PBMC_CITE/stim_cite')
meta.list <- list()
for(i in 1:4){
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
}
mm <- data.frame(batch = c(meta.list[[1]]$batch,meta.list[[2]]$batch,
                           meta.list[[3]]$batch,meta.list[[4]]$batch),
                 modality = c(meta.list[[1]]$modality,meta.list[[2]]$modality,
                              meta.list[[3]]$modality,meta.list[[4]]$modality),
                 condition = c(meta.list[[1]]$condition,meta.list[[2]]$condition,
                               meta.list[[3]]$condition,meta.list[[4]]$condition),
                 celltype = c(as.character(meta.list[[1]][,1]),
                              as.character(meta.list[[2]][,1]),
                              as.character(meta.list[[3]][,1]),
                              as.character(meta.list[[4]][,1]))
)

write.csv(as.matrix(a@Int.result[["bind"]]),
          "./Palette_DGE_drop.csv",
          row.names = FALSE)

adt.list <- list(adt1,adt2,adt3,adt4)
rownames(mm) <- colnames(Reduce(cbind,adt.list))
mm_sub <- dplyr::filter(mm,celltype == 'T')
idx <- match(rownames(mm_sub),rownames(mm))
library(data.table)
DGE_int <- a@Int.result[["bind"]][idx,]
write.csv(as.matrix(DGE_int),
          "./T_cell_analysis/DGE.csv",
          row.names = FALSE)

######################### Random Features ################################
source("./Palette/FUN.R")
folder_path <- "./Palette/utils"
files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
for (file in files) {
  source(file)
}
sourceCpp("./Palette/utils.cpp")
source("./Palette/utils.R")

setwd("./")
path <- c('./PBMC_ASAP/control_asap','./PBMC_ASAP/stim_asap',
          './PBMC_CITE/control_cite','./PBMC_CITE/stim_cite')
atac.list <- list()
adt.list <- list()
rna.list <- list()
meta.list <- list()
for(i in 1:4){
  meta.list[[i]] <- readRDS(paste0(path[i],'/meta.rds'))
  adt.list[[i]] <- readRDS(paste0(path[i],'/adt.rds'))
  if(i %in% c(1,2)){
    atac.list[[i]] <- readRDS(paste0(path[i],'/atac.rds'))
  }
  if(i %in% c(3,4)){
    rna.list[[i-2]] <- readRDS(paste0(path[i],'/rna.rds'))
  }
}

set.seed(123)
adt_random <- sample(1:227,84,replace = FALSE)
atac_random <- sample(1:99942,9548,replace = FALSE)
rna_random <- sample(1:2936,492,replace = FALSE)
for(i in 1:4){
  if(i %in% c(1,2)){
    adt.list[[i]] <- adt.list[[i]][-adt_random,]
    rna.list[[i]] <- rna.list[[i]][-rna_random,]
    atac.list[[i]] <- atac.list[[i]][-atac_random,]
  }else{
    adt.list[[i]] <- adt.list[[i]][-adt_random,]
  }
}

musa_obj <- Create.Musa.Object(data.list = c(adt.list,atac.list,rna.list), 
                               samples = c('control_asap','stim_asap',
                                           'control_cite','stim_cite',
                                           'control_asap','stim_asap',
                                           'control_cite','stim_cite'), 
                               modals = c('adt','adt','adt','adt',
                                          'atac','atac','rna','rna'),
                               filter=T,
                               min.cells = c(0,0,0),  
                               min.features = c(1,1,1),
                               modal.order = c("rna",'atac','adt'), 
                               sparce = TRUE)

musa_obj <- Normalize.Data(musa_obj, 
                           modal = "rna",
                           normal.method = "LogNormalize")
musa_obj <- Normalize.Data(musa_obj, 
                           modal = "adt",
                           normal.method = "CLR",
                           margin = 2)
musa_obj <- IDFLog_try(musa_obj,modal = "atac")

musa_obj <- Add.HVFs(musa_obj,modal = "rna")
musa_obj <- Add.HVFs(musa_obj, modal = "adt")
musa_obj <- Add.HVFs(musa_obj, modal = "atac")

musa_obj <- Find.Cluster(object = musa_obj, 
                         modal = c("rna",'atac','adt'),
                         sample = NULL,
                         dim.reduce = TRUE,# 是否在降维矩阵上计算
                         dims = list(1:50,1:50,1:50),# 使用哪些维数
                         method = c("PCA","PCA",'LSI'),
                         joint = F,
                         resolution = 1,#聚类参数
                         nn.k = 20,# knn snn参数
                         nn.method = "annoy",# knn snn参数
                         annoy.metric = "euclidean",
                         verbose = TRUE)

musa_obj <- Find.Subspace3(musa_obj,
                           modal = c("rna",'atac','adt'),
                           lambda = list(0.8,0.5,0.8),
                           alpha = list(0,0,0),
                           joint = T,
                           L2.norm = F,
                           sub.dims = list(1:40,1:30,1:30),
                           cos.dims = list(1:50,1:50,1:50),
                           Angle.var = c(15,20,15),
                           max.Angle = c(50,50,50))
a <- Run.Palette2(musa_obj,nn=2,nDims = 20,lambda = 0.5,emb_L2 = T)
mm <- data.frame(batch = c(meta.list[[1]]$batch,meta.list[[2]]$batch,
                           meta.list[[3]]$batch,meta.list[[4]]$batch),
                 modality = c(meta.list[[1]]$modality,meta.list[[2]]$modality,
                              meta.list[[3]]$modality,meta.list[[4]]$modality),
                 condition = c(meta.list[[1]]$condition,meta.list[[2]]$condition,
                               meta.list[[3]]$condition,meta.list[[4]]$condition),
                 celltype = c(as.character(meta.list[[1]][,1]),
                              as.character(meta.list[[2]][,1]),
                              as.character(meta.list[[3]][,1]),
                              as.character(meta.list[[4]][,1]))
)

write.csv(as.matrix(a@Int.result[["bind"]]),
          "./Palette_random_drop.csv",
          row.names = FALSE)

rownames(mm) <- colnames(Reduce(cbind,adt.list))
mm_sub <- dplyr::filter(mm,celltype == 'T')
idx <- match(rownames(mm_sub),rownames(mm))

library(data.table)
random_int <- a@Int.result[["bind"]][idx,]
write.csv(as.matrix(random_int),
          "./T_cell_analysis/random.csv",
          row.names = FALSE)
